package sampler

import (
	"fmt"
	"io"
	"math/rand"
	"os"
	"sort"
	"time"

	bamnative "github.com/rainoffallingstar/bamdriver-go/pkg/bamnative"
)

// ============================================================================
// 蓄水池算法实现
// ============================================================================

// NewReservoir 创建新的蓄水池
// capacity: 蓄水池容量 (即期望的采样数量)
func NewReservoir(capacity int64, seed int64) *Reservoir {
	if seed == 0 {
		seed = time.Now().UnixNano()
	}
	r := &Reservoir{
		capacity: capacity,
		records:  make([]interface{}, 0, capacity),
		rand:     rand.New(rand.NewSource(seed)),
	}
	return r
}

// Add 向蓄水池添加记录 (标准蓄水池算法)
//
// 算法伪代码:
// if size < capacity:
//
//	records[size] = item
//
// else:
//
//	j = random(0, size)  // [0, size]范围内随机
//	if j < capacity:
//	    records[j] = item
func (r *Reservoir) Add(record interface{}) {
	r.size++
	if r.size <= r.capacity {
		// 前capacity个直接放入
		r.records = append(r.records, record)
	} else {
		// 以 capacity/size 的概率替换
		j := r.rand.Int63n(r.size)
		if j < r.capacity {
			r.records[j] = record
		}
	}
}

// GetRecords 获取蓄水池中的所有记录
func (r *Reservoir) GetRecords() []interface{} {
	// 返回副本
	result := make([]interface{}, len(r.records))
	copy(result, r.records)
	return result
}

// Size 获取当前处理的记录数（即所有曾调用 Add 的次数，非蓄水池实际大小）
func (r *Reservoir) Size() int64 {
	return r.size
}

// Capacity 获取蓄水池容量
func (r *Reservoir) Capacity() int64 {
	return r.capacity
}

// ============================================================================
// BAM采样器主类型
// ============================================================================

// BAMSampler BAM采样器
type BAMSampler struct {
	config *SamplingConfig
	stats  *SamplingStats
}

// NewBAMSampler 创建新的采样器
func NewBAMSampler(config *SamplingConfig) *BAMSampler {
	return &BAMSampler{
		config: config,
		stats:  &SamplingStats{},
	}
}

// GetStats 获取采样统计信息
func (s *BAMSampler) GetStats() *SamplingStats {
	return s.stats
}

// ============================================================================
// 单BAM文件采样 (保证配对关系)
// ============================================================================

// SampleSingleBAM 单BAM文件采样 - 保证R1/R2配对关系
//
// 采用两遍流式策略以降低内存占用：
//   - 第一遍：仅扫描 read name，构建配对状态映射（内存 ≈ O(n_names × name_size)）
//   - 第二遍：只读取选中的 reads，内存中仅保留 k 条采样记录（O(k)）
func (s *BAMSampler) SampleSingleBAM(inputPath, outputPath string) error {
	// ── 第一遍：仅读取 read names，构建配对状态映射 ──────────────────────────
	infoMap, header, err := s.readPairedNamesOnly(inputPath)
	if err != nil {
		return fmt.Errorf("failed to read record names: %w", err)
	}

	// 构建配对名称列表，同时统计孤儿 reads
	var pairedNames []string
	for name, ni := range infoMap {
		if !ni.hasR1 || !ni.hasR2 {
			s.stats.OrphanReads++
			continue
		}
		if s.config.ProperPairsOnly && !ni.isProper {
			s.stats.FilteredByProperPair++
			continue
		}
		pairedNames = append(pairedNames, name)
	}
	sort.Strings(pairedNames)
	s.stats.TotalPairedInput = int64(len(pairedNames))

	// ── 蓄水池采样（以配对名称为单位）────────────────────────────────────────
	var sampleCount int64
	if s.config.Mode == ModeRatio {
		sampleCount = int64(float64(len(pairedNames)) * s.config.Ratio)
		if sampleCount == 0 && len(pairedNames) > 0 {
			fmt.Fprintf(os.Stderr, "Warning: ratio %.6f with %d pairs yields 0; output will be empty\n",
				s.config.Ratio, len(pairedNames))
		}
	} else {
		sampleCount = s.config.AbsoluteCount
		if sampleCount > int64(len(pairedNames)) {
			fmt.Fprintf(os.Stderr, "Warning: requested -count=%d exceeds available pairs=%d; clamping to %d\n",
				sampleCount, len(pairedNames), len(pairedNames))
			sampleCount = int64(len(pairedNames))
		}
	}

	reservoir := NewReservoir(sampleCount, s.config.Seed)
	for _, name := range pairedNames {
		reservoir.Add(name)
	}

	// 构建选中名称的 hash set
	sampledItems := reservoir.GetRecords()
	selectedNames := make(map[string]bool, len(sampledItems))
	for _, item := range sampledItems {
		selectedNames[item.(string)] = true
	}

	// ── 第二遍：仅读取选中的 reads（O(k) 内存）──────────────────────────────
	outputRecords, pairedOutput, err := s.collectSelectedRecords(inputPath, selectedNames)
	if err != nil {
		return fmt.Errorf("failed to collect selected records: %w", err)
	}
	s.stats.TotalPairedOutput = pairedOutput
	s.stats.TotalOutputRecords = int64(len(outputRecords))

	// ── 排序并更新 Header SortOrder ─────────────────────────────────────────
	switch s.config.SortOrder {
	case SortByCoord:
		sortRecordsByCoord(outputRecords)
		header.SortOrder = "coordinate"
	case SortByName:
		// 按 (name, R1优先) 排序，保证同名配对内 R1 在 R2 之前
		sort.Slice(outputRecords, func(i, j int) bool {
			ri, rj := outputRecords[i], outputRecords[j]
			if ri.Name != rj.Name {
				return ri.Name < rj.Name
			}
			return IsR1(ri) && !IsR1(rj)
		})
		header.SortOrder = "queryname"
	default:
		header.SortOrder = "unknown"
	}

	// ── 写入输出文件 ──────────────────────────────────────────────────────────
	if err := s.writeOutput(outputRecords, outputPath, header); err != nil {
		return fmt.Errorf("failed to write output: %w", err)
	}

	return nil
}

// readPairedNamesOnly 第一遍扫描：仅读取 read name 与标志，构建配对状态映射。
// 内存占用约为 O(n_names × (name_len + overhead))，远低于全量加载所有 records。
func (s *BAMSampler) readPairedNamesOnly(path string) (map[string]*nameInfo, *bamnative.Header, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	defer f.Close()

	reader, err := bamnative.NewReader(f)
	if err != nil {
		return nil, nil, err
	}

	header := reader.Header()
	infoMap := make(map[string]*nameInfo)
	var count int64

	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, nil, err
		}

		count++
		s.stats.TotalInputRecords++
		if count%1_000_000 == 0 {
			fmt.Fprintf(os.Stderr, "  [Pass 1] Processed %dM records...\n", count/1_000_000)
		}

		// 跳过 secondary / supplementary
		if rec.Flags&bamnative.FlagSecondary != 0 || rec.Flags&bamnative.FlagSupplementary != 0 {
			s.stats.SecondaryDropped++
			continue
		}

		ni, ok := infoMap[rec.Name]
		if !ok {
			ni = &nameInfo{}
			infoMap[rec.Name] = ni
		}

		if IsR1(rec) {
			if ni.hasR1 {
				s.stats.MultimapDropped++
			}
			ni.hasR1 = true
		} else if IsR2(rec) {
			if ni.hasR2 {
				s.stats.MultimapDropped++
			}
			ni.hasR2 = true
		} else {
			s.stats.UnmarkedReads++
		}

		if IsProperPair(rec) {
			ni.isProper = true
		}
	}

	return infoMap, header, nil
}

// collectSelectedRecords 第二遍扫描：从 BAM 文件中只收集名称在 selectedNames 中的 reads。
// 每个名称的第一条 R1 和第一条 R2 各保留一条，与第一遍的去多比对逻辑保持一致。
func (s *BAMSampler) collectSelectedRecords(inputPath string, selectedNames map[string]bool) ([]*bamnative.Record, int64, error) {
	f, err := os.Open(inputPath)
	if err != nil {
		return nil, 0, err
	}
	defer f.Close()

	reader, err := bamnative.NewReader(f)
	if err != nil {
		return nil, 0, err
	}

	seenR1 := make(map[string]bool, len(selectedNames))
	seenR2 := make(map[string]bool, len(selectedNames))
	outputRecords := make([]*bamnative.Record, 0, len(selectedNames)*2)
	var count int64

	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, 0, err
		}

		count++
		if count%1_000_000 == 0 {
			fmt.Fprintf(os.Stderr, "  [Pass 2] Processed %dM records...\n", count/1_000_000)
		}

		if !selectedNames[rec.Name] {
			continue
		}
		if rec.Flags&bamnative.FlagSecondary != 0 || rec.Flags&bamnative.FlagSupplementary != 0 {
			continue
		}

		if IsR1(rec) && !seenR1[rec.Name] {
			seenR1[rec.Name] = true
			outputRecords = append(outputRecords, rec)
		} else if IsR2(rec) && !seenR2[rec.Name] {
			seenR2[rec.Name] = true
			outputRecords = append(outputRecords, rec)
		}
	}

	validSelected := make(map[string]bool, len(selectedNames))
	for name := range selectedNames {
		if seenR1[name] && seenR2[name] {
			validSelected[name] = true
		}
	}

	// 极端情况下（输入文件在两遍读取间发生变化）保护一致性
	if len(validSelected) != len(selectedNames) {
		fmt.Fprintf(os.Stderr, "Warning: selected pair count changed between scan and collect passes (%d -> %d)\n",
			len(selectedNames), len(validSelected))
		outputRecords = filterRecordsByNameSet(outputRecords, validSelected)
	}

	return outputRecords, int64(len(validSelected)), nil
}

// ============================================================================
// 排序辅助函数
// ============================================================================

// sortRecordsByCoord 按各 record 自身的参考序列ID和坐标严格排序，保证 BAI 索引有效
func sortRecordsByCoord(records []*bamnative.Record) {
	sort.Slice(records, func(i, j int) bool {
		ri, rj := records[i], records[j]
		if ri.RefID != rj.RefID {
			return ri.RefID < rj.RefID
		}
		return ri.Pos < rj.Pos
	})
}

// ============================================================================
// 内部IO辅助函数
// ============================================================================

// writeOutput 写入输出BAM文件 [DESIGN-1: 参数改为 []*bamnative.Record]
func (s *BAMSampler) writeOutput(records []*bamnative.Record, outputPath string, header *bamnative.Header) error {
	writer, err := bamnative.NewWriter(outputPath, header)
	if err != nil {
		return fmt.Errorf("failed to create BAM writer: %w", err)
	}
	defer writer.Close()

	for i, rec := range records {
		if (i+1)%1_000_000 == 0 {
			fmt.Fprintf(os.Stderr, "  Writing %dM records...\n", (i+1)/1_000_000)
		}
		if err := writer.Write(rec); err != nil {
			return fmt.Errorf("failed to write record: %w", err)
		}
	}

	return nil
}
