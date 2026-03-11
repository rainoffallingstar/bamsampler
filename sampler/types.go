package sampler

import (
	"math/rand"

	bamnative "github.com/rainoffallingstar/bamdriver-go/pkg/bamnative"
)

// ============================================================================
// 采样配置类型
// ============================================================================

// SamplingMode 采样模式
type SamplingMode int

const (
	ModeRatio    SamplingMode = iota // 比例采样
	ModeAbsolute                     // 绝对数量采样
)

// InputMode 输入模式
type InputMode int

const (
	InputSingleBAM InputMode = iota // 单BAM文件
	InputDualBAM                    // 双BAM (R1+R2)
)

// SortOrder 排序方式常量
const (
	SortByName  = "name"
	SortByCoord = "coord"
	SortNone    = "none"
)

// SamplingConfig 采样配置
type SamplingConfig struct {
	Mode            SamplingMode // 采样模式
	InputMode       InputMode    // 输入模式
	Ratio           float64      // 比例 (0.0-1.0)
	AbsoluteCount   int64        // 绝对数量
	Seed            int64        // 随机种子 (可选)
	ProperPairsOnly bool         // 只保留properly paired
	SortOrder       string       // SortByName 或 SortByCoord
	CreateIndex     bool         // 创建BAI索引
}

// ============================================================================
// 蓄水池类型
// ============================================================================

// Reservoir 蓄水池结构
// 非并发安全：所有调用方均在单 goroutine 内，无需加锁。
type Reservoir struct {
	capacity int64
	size     int64
	records  []interface{} // 使用interface{}来存储采样项（当前为read name）
	rand     *rand.Rand
}

// ============================================================================
// 统计信息类型
// ============================================================================

// SamplingStats 采样统计信息
type SamplingStats struct {
	TotalInputRecords    int64 // 输入总记录数
	TotalPairedInput     int64 // 输入配对数
	TotalOutputRecords   int64 // 输出总记录数
	TotalPairedOutput    int64 // 输出配对数
	FilteredByProperPair int64 // 因未通过 proper-pair 过滤而跳过的 pair 数
	OrphanReads          int64 // 在对侧找不到配对的孤儿 pair 数（single/dual 均统计）
	UnmarkedReads        int64 // 单BAM模式中既非R1也非R2的记录数（静默丢弃）
	MultimapDropped      int64 // 因多比对保留第一条而丢弃的额外记录数
	SecondaryDropped     int64 // secondary/supplementary reads 被过滤的记录数
}

// ============================================================================
// 内部流式采样辅助类型
// ============================================================================

// nameInfo 第一遍扫描时按 read name 记录的最小状态
// 仅用于 SampleSingleBAM 两遍流式采样的内部实现，不对外暴露。
type nameInfo struct {
	hasR1    bool // 该 name 在 BAM 中至少出现过一条 R1
	hasR2    bool // 该 name 在 BAM 中至少出现过一条 R2
	isProper bool // 该 name 的至少一条 read 通过 proper-pair 检查
}

// ============================================================================
// 辅助函数
// ============================================================================

// IsProperPair 检查是否是proper pair
func IsProperPair(record *bamnative.Record) bool {
	return record.Flags&bamnative.FlagPaired != 0 &&
		record.Flags&bamnative.FlagProperPair != 0 &&
		record.Flags&bamnative.FlagUnmapped == 0 &&
		record.Flags&bamnative.FlagMateUnmapped == 0
}

// IsR1 检查是否是R1
func IsR1(record *bamnative.Record) bool {
	return record.Flags&bamnative.FlagFirstInPair != 0
}

// IsR2 检查是否是R2
func IsR2(record *bamnative.Record) bool {
	return record.Flags&bamnative.FlagSecondInPair != 0
}
