package main

import (
	"flag"
	"fmt"
	"os"
	"time"

	"github.com/rainoffallingstar/bamdriver-go/pkg/bamnative"

	"bamsampler/sampler"
)

// CLI参数结构体
type CLIArgs struct {
	// 采样参数 (互斥,二选一)
	Ratio         float64
	AbsoluteCount int64
	Seed          int64

	// 输入模式
	InputMode string

	// 过滤选项
	ProperPairsOnly bool

	// 排序和索引
	SortOrder   string
	CreateIndex bool

	// 位置参数
	InputFiles []string
}

func main() {
	args := parseArgs()

	// 验证参数
	if err := validateArgs(args); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		usage()
		os.Exit(1)
	}

	// 创建采样配置
	config := createConfig(args)

	// 创建采样器
	s := sampler.NewBAMSampler(config)

	// 记录开始时间
	startTime := time.Now()

	// 根据输入模式执行采样
	var err error
	if args.InputMode == "dual" {
		// 双BAM模式：使用 SampleDualBAMSorted（更可靠的配对逻辑）[BUG-3]
		err = s.SampleDualBAMSorted(args.InputFiles[0], args.InputFiles[1], args.InputFiles[2])
	} else {
		// 单BAM模式
		err = s.SampleSingleBAM(args.InputFiles[0], args.InputFiles[1])
	}

	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	// 输出统计信息
	stats := s.GetStats()
	elapsed := time.Since(startTime)

	fmt.Printf("=== Sampling Complete ===\n")
	fmt.Printf("Total input records:  %d\n", stats.TotalInputRecords)
	if stats.TotalPairedInput > 0 {
		fmt.Printf("Total input pairs:    %d\n", stats.TotalPairedInput)
	}
	fmt.Printf("Total output records: %d\n", stats.TotalOutputRecords)
	if stats.TotalPairedOutput > 0 {
		fmt.Printf("Total output pairs:   %d\n", stats.TotalPairedOutput)
	}
	if stats.FilteredByProperPair > 0 {
		fmt.Printf("Filtered pairs (not proper pair): %d\n", stats.FilteredByProperPair)
	}
	if stats.OrphanReads > 0 {
		fmt.Printf("Orphan pairs (no mate found): %d\n", stats.OrphanReads)
		// 超过 5% 时发出 WARNING，提醒用户采样比例可能偏差
		totalEligible := stats.TotalPairedInput + stats.OrphanReads
		if totalEligible > 0 {
			orphanRate := float64(stats.OrphanReads) / float64(totalEligible)
			if orphanRate > 0.05 {
				fmt.Fprintf(os.Stderr, "WARNING: High orphan pair rate: %d pairs (%.1f%%) have no mate. Sampling accuracy may be affected.\n",
					stats.OrphanReads, orphanRate*100)
			}
		}
	}
	if stats.UnmarkedReads > 0 {
		fmt.Printf("Warning: unmarked reads (neither R1 nor R2, dropped): %d\n", stats.UnmarkedReads)
	}
	if stats.MultimapDropped > 0 {
		fmt.Printf("Warning: multimapper records dropped (only first alignment kept): %d\n", stats.MultimapDropped)
	}
	if stats.SecondaryDropped > 0 {
		fmt.Printf("Secondary/supplementary reads filtered: %d\n", stats.SecondaryDropped)
	}
	fmt.Printf("Time elapsed:         %v\n", elapsed)

	// 如果需要创建索引
	indexFailed := false
	if args.CreateIndex {
		if args.InputMode == "dual" {
			outputPrefix := args.InputFiles[2]
			for _, suffix := range []string{"_R1.bam", "_R2.bam"} {
				path := outputPrefix + suffix
				if err := createIndex(path); err != nil {
					fmt.Fprintf(os.Stderr, "Error: failed to create index for %s: %v\n", path, err)
					indexFailed = true
				}
			}
		} else {
			outputPath := args.InputFiles[1]
			if err := createIndex(outputPath); err != nil {
				fmt.Fprintf(os.Stderr, "Error: failed to create index: %v\n", err)
				indexFailed = true
			}
		}
	}
	if indexFailed {
		os.Exit(1)
	}
}

func parseArgs() *CLIArgs {
	args := &CLIArgs{}

	flag.Float64Var(&args.Ratio, "ratio", 0, "Sampling ratio (0.0-1.0)")
	flag.Int64Var(&args.AbsoluteCount, "count", 0, "Absolute number of read pairs to sample")
	flag.Int64Var(&args.Seed, "seed", 0, "Random seed (optional, for reproducible sampling)")
	flag.StringVar(&args.InputMode, "mode", "single", "Input mode: single or dual")
	flag.BoolVar(&args.ProperPairsOnly, "proper-pairs-only", false, "Only keep properly paired reads")
	flag.StringVar(&args.SortOrder, "sort", "none", "Sort order for output: name, coord, or none (default: none)")
	flag.BoolVar(&args.CreateIndex, "index", false, "Create BAI index for output")

	flag.Usage = usage
	flag.Parse()

	args.InputFiles = flag.Args()
	return args
}

func validateArgs(args *CLIArgs) error {
	// 检查采样参数
	if args.Ratio == 0 && args.AbsoluteCount == 0 {
		return fmt.Errorf("either -ratio or -count must be specified")
	}
	if args.Ratio != 0 && args.AbsoluteCount != 0 {
		return fmt.Errorf("only one of -ratio or -count can be specified")
	}

	// 验证ratio范围
	if args.Ratio < 0 || args.Ratio > 1 {
		return fmt.Errorf("ratio must be between 0.0 and 1.0")
	}

	// 验证count范围
	if args.AbsoluteCount < 0 {
		return fmt.Errorf("count must be positive")
	}

	// 验证输入模式
	if args.InputMode != "single" && args.InputMode != "dual" {
		return fmt.Errorf("input mode must be 'single' or 'dual'")
	}

	// 验证位置参数
	if err := validatePositionalArgs(args); err != nil {
		return err
	}

	// 验证排序方式
	if args.SortOrder != sampler.SortByName && args.SortOrder != sampler.SortByCoord && args.SortOrder != sampler.SortNone {
		return fmt.Errorf("sort order must be 'name', 'coord', or 'none'")
	}

	// 若 -index 但未指定 -sort coord，BAI 索引要求坐标排序，直接报错
	if args.CreateIndex && args.SortOrder != sampler.SortByCoord {
		return fmt.Errorf("-index requires -sort coord")
	}

	return nil
}

func validatePositionalArgs(args *CLIArgs) error {
	if args.InputMode == "dual" {
		if len(args.InputFiles) != 3 {
			return fmt.Errorf("dual mode requires 3 positional arguments: <R1.bam> <R2.bam> <output_prefix>")
		}
		return nil
	}
	if len(args.InputFiles) != 2 {
		return fmt.Errorf("single mode requires 2 positional arguments: <input.bam> <output.bam>")
	}
	return nil
}

func createConfig(args *CLIArgs) *sampler.SamplingConfig {
	config := &sampler.SamplingConfig{
		InputMode:       sampler.InputSingleBAM,
		Ratio:           args.Ratio,
		AbsoluteCount:   args.AbsoluteCount,
		Seed:            args.Seed,
		ProperPairsOnly: args.ProperPairsOnly,
		SortOrder:       args.SortOrder,
		CreateIndex:     args.CreateIndex,
	}

	if args.Ratio != 0 {
		config.Mode = sampler.ModeRatio
	} else {
		config.Mode = sampler.ModeAbsolute
	}

	if args.InputMode == "dual" {
		config.InputMode = sampler.InputDualBAM
	}

	return config
}

func createIndex(bamPath string) error {
	fmt.Printf("Creating BAI index for %s...\n", bamPath)
	return bamnative.BuildIndex(bamPath)
}

func usage() {
	fmt.Fprintf(os.Stderr, `BAM Sampler - Random sampling of BAM files using reservoir sampling

Usage:
  bamsampler [options] <input.bam> <output.bam>
  bamsampler [options] -mode dual <R1.bam> <R2.bam> <output_prefix>

Sampling options:
  -ratio FLOAT          Sampling ratio (0.0-1.0), e.g., 0.1 for 10%%
  -count INT            Absolute number of read pairs to sample
  -seed INT             Random seed (optional, for reproducible results)

Input options:
  -mode string          Input mode: single or dual (default: single)

Filtering options:
  -proper-pairs-only    Only keep properly paired reads

Output options:
  -sort string          Sort order for output: name, coord, or none (default: none)
  -index                Create BAI index for output

Examples:
  # Sample 10%% of reads from a single BAM file
  bamsampler -ratio 0.1 input.bam output.bam

  # Sample 1 million read pairs from a single BAM file
  bamsampler -count 1000000 input.bam output.bam

  # Sample 10%% of reads from paired R1/R2 BAM files
  bamsampler -ratio 0.1 -mode dual R1.bam R2.bam output

  # Sample with specific random seed for reproducibility
  bamsampler -ratio 0.1 -seed 42 input.bam output.bam

  # Only keep properly paired reads
  bamsampler -ratio 0.1 -proper-pairs-only input.bam output.bam

`)
	flag.PrintDefaults()
}
