package sampler

import (
	"fmt"
	"io"
	"os"
	"sort"

	bamnative "github.com/rainoffallingstar/bamdriver-go/pkg/bamnative"
)

// ============================================================================
// 双BAM文件采样 (排序匹配版本)
// ============================================================================

// SampleDualBAMSorted 双BAM采样 - CLI的默认调用路径
func (s *BAMSampler) SampleDualBAMSorted(r1Path, r2Path, outputPrefix string) error {
	infoMap := make(map[string]*dualNameInfo)

	// Pass 1: 扫描 R1 和 R2，仅记录 name 级统计状态
	r1Header, _, err := s.scanDualFilePass1(r1Path, true, infoMap)
	if err != nil {
		return fmt.Errorf("failed to scan R1: %w", err)
	}
	r2Header, _, err := s.scanDualFilePass1(r2Path, false, infoMap)
	if err != nil {
		return fmt.Errorf("failed to scan R2: %w", err)
	}

	// 统计配对、孤儿与多比对信息（按 pair-name 口径）
	var pairedNames []string
	for name, info := range infoMap {
		if info.r1Count > 0 && info.r2Count > 0 {
			if s.config.ProperPairsOnly && !(info.r1Proper && info.r2Proper) {
				s.stats.FilteredByProperPair++
				continue
			}
			pairedNames = append(pairedNames, name)
			if info.r1Count > 1 {
				s.stats.MultimapDropped += info.r1Count - 1
			}
			if info.r2Count > 1 {
				s.stats.MultimapDropped += info.r2Count - 1
			}
			continue
		}
		// 单侧出现的 read name 记为 1 个 orphan pair
		if info.r1Count > 0 || info.r2Count > 0 {
			s.stats.OrphanReads++
		}
	}
	sort.Strings(pairedNames)
	s.stats.TotalPairedInput = int64(len(pairedNames))

	// 计算采样数量
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

	// 蓄水池采样（按配对 name）
	reservoir := NewReservoir(sampleCount, s.config.Seed)
	for _, name := range pairedNames {
		reservoir.Add(name)
	}

	sampledNamesRaw := reservoir.GetRecords()
	selectedNames := make(map[string]bool, len(sampledNamesRaw))
	for _, item := range sampledNamesRaw {
		selectedNames[item.(string)] = true
	}

	// Pass 2: 分别收集 R1/R2 中被选中的首条主比对记录
	r1Out, r1Seen, err := s.collectDualSelectedPass2(r1Path, selectedNames)
	if err != nil {
		return fmt.Errorf("failed to collect R1 selected records: %w", err)
	}
	r2Out, r2Seen, err := s.collectDualSelectedPass2(r2Path, selectedNames)
	if err != nil {
		return fmt.Errorf("failed to collect R2 selected records: %w", err)
	}

	validSelected := make(map[string]bool, len(selectedNames))
	for name := range selectedNames {
		if r1Seen[name] && r2Seen[name] {
			validSelected[name] = true
		}
	}

	// 极端情况下（输入文件在两遍读取间发生变化）保护一致性
	if len(validSelected) != len(selectedNames) {
		fmt.Fprintf(os.Stderr, "Warning: selected pair count changed between scan and collect passes (%d -> %d)\n",
			len(selectedNames), len(validSelected))
		r1Out = filterRecordsByNameSet(r1Out, validSelected)
		r2Out = filterRecordsByNameSet(r2Out, validSelected)
	}

	s.stats.TotalPairedOutput = int64(len(validSelected))
	s.stats.TotalOutputRecords = s.stats.TotalPairedOutput * 2

	// 写入输出文件
	r1Output := outputPrefix + "_R1.bam"
	r2Output := outputPrefix + "_R2.bam"

	if s.config.SortOrder == SortByCoord {
		sortRecordsByCoord(r1Out)
		sortRecordsByCoord(r2Out)
		r1Header.SortOrder = "coordinate"
		r2Header.SortOrder = "coordinate"
	} else if s.config.SortOrder == SortByName {
		sortRecordsByName(r1Out)
		sortRecordsByName(r2Out)
		r1Header.SortOrder = "queryname"
		r2Header.SortOrder = "queryname"
	} else {
		r1Header.SortOrder = "unknown"
		r2Header.SortOrder = "unknown"
	}

	if err := s.writeOutput(r1Out, r1Output, r1Header); err != nil {
		return fmt.Errorf("failed to write R1 output: %w", err)
	}
	if err := s.writeOutput(r2Out, r2Output, r2Header); err != nil {
		return fmt.Errorf("failed to write R2 output: %w", err)
	}

	return nil
}

type dualNameInfo struct {
	r1Count  int64
	r2Count  int64
	r1Proper bool
	r2Proper bool
}

type dualFileScanStats struct {
	r1LikeReads int64
	r2LikeReads int64
}

func (s *BAMSampler) scanDualFilePass1(path string, isR1 bool, infoMap map[string]*dualNameInfo) (*bamnative.Header, *dualFileScanStats, error) {
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
	roleStats := &dualFileScanStats{}
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
			fmt.Fprintf(os.Stderr, "  [Pass 1] Processed %dM records from %s...\n", count/1_000_000, path)
		}

		if rec.Flags&bamnative.FlagSecondary != 0 || rec.Flags&bamnative.FlagSupplementary != 0 {
			s.stats.SecondaryDropped++
			continue
		}

		if IsR1(rec) && !IsR2(rec) {
			roleStats.r1LikeReads++
		} else if IsR2(rec) && !IsR1(rec) {
			roleStats.r2LikeReads++
		}

		info, ok := infoMap[rec.Name]
		if !ok {
			info = &dualNameInfo{}
			infoMap[rec.Name] = info
		}
		if isR1 {
			info.r1Count++
			if IsProperPair(rec) {
				info.r1Proper = true
			}
		} else {
			info.r2Count++
			if IsProperPair(rec) {
				info.r2Proper = true
			}
		}
	}

	if err := validateDualFileRole(path, isR1, roleStats); err != nil {
		return nil, nil, err
	}

	return header, roleStats, nil
}

func validateDualFileRole(path string, isR1 bool, roleStats *dualFileScanStats) error {
	if roleStats == nil {
		return nil
	}

	if isR1 {
		if roleStats.r2LikeReads > roleStats.r1LikeReads {
			return fmt.Errorf("input file %s appears to contain more R2 than R1 reads (%d vs %d)", path, roleStats.r2LikeReads, roleStats.r1LikeReads)
		}
		if roleStats.r2LikeReads > 0 {
			fmt.Fprintf(os.Stderr, "Warning: detected %d R2-like reads in R1 input file %s\n", roleStats.r2LikeReads, path)
		}
		return nil
	}

	if roleStats.r1LikeReads > roleStats.r2LikeReads {
		return fmt.Errorf("input file %s appears to contain more R1 than R2 reads (%d vs %d)", path, roleStats.r1LikeReads, roleStats.r2LikeReads)
	}
	if roleStats.r1LikeReads > 0 {
		fmt.Fprintf(os.Stderr, "Warning: detected %d R1-like reads in R2 input file %s\n", roleStats.r1LikeReads, path)
	}
	return nil
}

func (s *BAMSampler) collectDualSelectedPass2(path string, selectedNames map[string]bool) ([]*bamnative.Record, map[string]bool, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	defer f.Close()

	reader, err := bamnative.NewReader(f)
	if err != nil {
		return nil, nil, err
	}

	seen := make(map[string]bool, len(selectedNames))
	output := make([]*bamnative.Record, 0, len(selectedNames))
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
		if count%1_000_000 == 0 {
			fmt.Fprintf(os.Stderr, "  [Pass 2] Processed %dM records from %s...\n", count/1_000_000, path)
		}

		if !selectedNames[rec.Name] || seen[rec.Name] {
			continue
		}
		if rec.Flags&bamnative.FlagSecondary != 0 || rec.Flags&bamnative.FlagSupplementary != 0 {
			continue
		}
		if s.config.ProperPairsOnly && !IsProperPair(rec) {
			continue
		}

		seen[rec.Name] = true
		output = append(output, rec)
	}

	return output, seen, nil
}

func filterRecordsByNameSet(records []*bamnative.Record, keep map[string]bool) []*bamnative.Record {
	filtered := make([]*bamnative.Record, 0, len(records))
	for _, rec := range records {
		if keep[rec.Name] {
			filtered = append(filtered, rec)
		}
	}
	return filtered
}

func sortRecordsByName(records []*bamnative.Record) {
	sort.Slice(records, func(i, j int) bool {
		return records[i].Name < records[j].Name
	})
}
