package sampler

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	bamnative "github.com/rainoffallingstar/bamdriver-go/pkg/bamnative"
)

// ============================================================================
// Integration-test helpers
// ============================================================================

const (
	// Minimum flags for R1 / R2 primary reads.
	itFlagR1 = bamnative.FlagPaired | bamnative.FlagFirstInPair
	itFlagR2 = bamnative.FlagPaired | bamnative.FlagSecondInPair
)

// itMakeRecord builds a minimal but round-trip-safe BAM record.
// CIGAR and Seq/Qual are kept consistent so the native reader can re-parse the file.
func itMakeRecord(name string, flags uint16) *bamnative.Record {
	return &bamnative.Record{
		Name:      name,
		Flags:     flags,
		RefID:     -1,
		Pos:       0,
		MapQ:      60,
		Cigar:     []bamnative.CigarOp{{Op: 'M', Len: 4}},
		MateRefID: -1,
		MatePos:   0,
		TLen:      0,
		Seq:       "ACGT",
		Qual:      []byte{30, 30, 30, 30},
		Aux:       []*bamnative.AuxField{},
	}
}

// itWriteBAM writes records to a new temp BAM file and returns its path.
func itWriteBAM(t *testing.T, records []*bamnative.Record) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "test.bam")
	header := &bamnative.Header{
		SortOrder:  "unknown",
		OtherLines: make(map[string]string),
	}
	w, err := bamnative.NewWriter(path, header)
	if err != nil {
		t.Fatalf("itWriteBAM: NewWriter: %v", err)
	}
	for _, rec := range records {
		if err := w.Write(rec); err != nil {
			_ = w.Close()
			t.Fatalf("itWriteBAM: Write: %v", err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("itWriteBAM: Close: %v", err)
	}
	return path
}

// itReadBAM reads all records from a BAM file.
func itReadBAM(t *testing.T, path string) []*bamnative.Record {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("itReadBAM: Open %s: %v", path, err)
	}
	defer f.Close()
	reader, err := bamnative.NewReader(f)
	if err != nil {
		t.Fatalf("itReadBAM: NewReader: %v", err)
	}
	var records []*bamnative.Record
	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("itReadBAM: Read: %v", err)
		}
		records = append(records, rec)
	}
	return records
}

func itRecordSignatures(records []*bamnative.Record) []string {
	signatures := make([]string, 0, len(records))
	for _, rec := range records {
		signatures = append(signatures, fmt.Sprintf("%s|%d|%d|%d", rec.Name, rec.Flags, rec.RefID, rec.Pos))
	}
	return signatures
}

// itConfig returns a default SamplingConfig for tests.
func itConfig() *SamplingConfig {
	return &SamplingConfig{
		Mode:            ModeRatio,
		InputMode:       InputSingleBAM,
		Ratio:           0.5,
		Seed:            42,
		ProperPairsOnly: false,
		SortOrder:       SortNone,
		CreateIndex:     false,
	}
}

// ============================================================================
// SampleSingleBAM integration tests
// ============================================================================

// TestIntegration_SingleBAM_RatioSampling verifies that ratio sampling produces
// the expected number of output pairs.
func TestIntegration_SingleBAM_RatioSampling(t *testing.T) {
	const nPairs = 100
	var records []*bamnative.Record
	for i := 0; i < nPairs; i++ {
		name := fmt.Sprintf("read%04d", i)
		records = append(records, itMakeRecord(name, itFlagR1))
		records = append(records, itMakeRecord(name, itFlagR2))
	}

	inputPath := itWriteBAM(t, records)
	outputPath := filepath.Join(t.TempDir(), "out.bam")

	s := NewBAMSampler(itConfig()) // ratio=0.5, seed=42
	if err := s.SampleSingleBAM(inputPath, outputPath); err != nil {
		t.Fatalf("SampleSingleBAM: %v", err)
	}

	stats := s.GetStats()
	wantPairs := int64(nPairs / 2) // reservoir capacity = floor(100 * 0.5) = 50
	if stats.TotalPairedOutput != wantPairs {
		t.Errorf("TotalPairedOutput = %d; want %d", stats.TotalPairedOutput, wantPairs)
	}
	if stats.TotalOutputRecords != wantPairs*2 {
		t.Errorf("TotalOutputRecords = %d; want %d", stats.TotalOutputRecords, wantPairs*2)
	}
}

// TestIntegration_SingleBAM_PairIntegrity verifies every sampled R1 has a matching R2.
func TestIntegration_SingleBAM_PairIntegrity(t *testing.T) {
	const nPairs = 20
	var records []*bamnative.Record
	for i := 0; i < nPairs; i++ {
		name := fmt.Sprintf("read%04d", i)
		records = append(records, itMakeRecord(name, itFlagR1))
		records = append(records, itMakeRecord(name, itFlagR2))
	}

	inputPath := itWriteBAM(t, records)
	outputPath := filepath.Join(t.TempDir(), "out.bam")

	s := NewBAMSampler(itConfig())
	if err := s.SampleSingleBAM(inputPath, outputPath); err != nil {
		t.Fatalf("SampleSingleBAM: %v", err)
	}

	outRecs := itReadBAM(t, outputPath)
	r1Names := make(map[string]bool)
	r2Names := make(map[string]bool)
	for _, rec := range outRecs {
		if IsR1(rec) {
			r1Names[rec.Name] = true
		} else if IsR2(rec) {
			r2Names[rec.Name] = true
		}
	}
	for name := range r1Names {
		if !r2Names[name] {
			t.Errorf("R1 %q in output has no matching R2", name)
		}
	}
	for name := range r2Names {
		if !r1Names[name] {
			t.Errorf("R2 %q in output has no matching R1", name)
		}
	}
}

// TestIntegration_SingleBAM_MultimapFirstKeep verifies that only the first alignment for a
// read name is retained and extras are counted in MultimapDropped.
func TestIntegration_SingleBAM_MultimapFirstKeep(t *testing.T) {
	records := []*bamnative.Record{
		itMakeRecord("dupread", itFlagR1), // first alignment – kept
		itMakeRecord("dupread", itFlagR1), // duplicate – should be dropped
		itMakeRecord("dupread", itFlagR2),
	}

	inputPath := itWriteBAM(t, records)
	outputPath := filepath.Join(t.TempDir(), "out.bam")

	cfg := itConfig()
	cfg.Ratio = 1.0 // keep all pairs
	s := NewBAMSampler(cfg)
	if err := s.SampleSingleBAM(inputPath, outputPath); err != nil {
		t.Fatalf("SampleSingleBAM: %v", err)
	}

	stats := s.GetStats()
	if stats.MultimapDropped != 1 {
		t.Errorf("MultimapDropped = %d; want 1", stats.MultimapDropped)
	}
	if stats.TotalPairedOutput != 1 {
		t.Errorf("TotalPairedOutput = %d; want 1", stats.TotalPairedOutput)
	}

	outRecs := itReadBAM(t, outputPath)
	if len(outRecs) != 2 { // 1 R1 + 1 R2
		t.Errorf("output record count = %d; want 2", len(outRecs))
	}
}

// TestIntegration_SingleBAM_SecondaryFilter verifies secondary and supplementary reads
// are excluded and counted in SecondaryDropped.
func TestIntegration_SingleBAM_SecondaryFilter(t *testing.T) {
	flagSecondary := uint16(itFlagR1 | bamnative.FlagSecondary)
	flagSupplementary := uint16(itFlagR1 | bamnative.FlagSupplementary)

	// read1: primary pair (both kept)
	// read1 R1 secondary: filtered
	// read2 R1 supplementary: filtered  → read2 R2 becomes orphan, not output
	records := []*bamnative.Record{
		itMakeRecord("read1", itFlagR1),
		itMakeRecord("read1", itFlagR2),
		itMakeRecord("read1", flagSecondary),
		itMakeRecord("read2", flagSupplementary),
		itMakeRecord("read2", itFlagR2),
	}

	inputPath := itWriteBAM(t, records)
	outputPath := filepath.Join(t.TempDir(), "out.bam")

	cfg := itConfig()
	cfg.Ratio = 1.0
	s := NewBAMSampler(cfg)
	if err := s.SampleSingleBAM(inputPath, outputPath); err != nil {
		t.Fatalf("SampleSingleBAM: %v", err)
	}

	stats := s.GetStats()
	if stats.SecondaryDropped != 2 {
		t.Errorf("SecondaryDropped = %d; want 2", stats.SecondaryDropped)
	}
	if stats.TotalPairedOutput != 1 {
		t.Errorf("TotalPairedOutput = %d; want 1 (only read1 pair survives)", stats.TotalPairedOutput)
	}
}

// ============================================================================
// SampleDualBAMSorted integration tests
// ============================================================================

// TestIntegration_DualBAM_PairIntegrity verifies dual-mode output contains complete pairs.
func TestIntegration_DualBAM_PairIntegrity(t *testing.T) {
	const nPairs = 20
	var r1Recs, r2Recs []*bamnative.Record
	for i := 0; i < nPairs; i++ {
		name := fmt.Sprintf("read%04d", i)
		r1Recs = append(r1Recs, itMakeRecord(name, itFlagR1))
		r2Recs = append(r2Recs, itMakeRecord(name, itFlagR2))
	}

	r1Path := itWriteBAM(t, r1Recs)
	r2Path := itWriteBAM(t, r2Recs)
	outDir := t.TempDir()
	prefix := filepath.Join(outDir, "out")

	s := NewBAMSampler(itConfig()) // ratio=0.5
	if err := s.SampleDualBAMSorted(r1Path, r2Path, prefix); err != nil {
		t.Fatalf("SampleDualBAMSorted: %v", err)
	}

	r1Out := itReadBAM(t, prefix+"_R1.bam")
	r2Out := itReadBAM(t, prefix+"_R2.bam")

	if len(r1Out) != len(r2Out) {
		t.Errorf("output R1 count %d != R2 count %d", len(r1Out), len(r2Out))
	}

	r1Names := make(map[string]bool)
	for _, rec := range r1Out {
		r1Names[rec.Name] = true
	}
	for _, rec := range r2Out {
		if !r1Names[rec.Name] {
			t.Errorf("R2 %q in output has no matching R1", rec.Name)
		}
	}
}

func TestIntegration_SingleBAM_DeterministicSeed(t *testing.T) {
	const nPairs = 120
	var records []*bamnative.Record
	for i := 0; i < nPairs; i++ {
		name := fmt.Sprintf("read%04d", i)
		records = append(records, itMakeRecord(name, itFlagR1))
		records = append(records, itMakeRecord(name, itFlagR2))
	}

	inputPath := itWriteBAM(t, records)
	out1 := filepath.Join(t.TempDir(), "out1.bam")
	out2 := filepath.Join(t.TempDir(), "out2.bam")

	cfg := itConfig()
	cfg.Ratio = 0.37
	cfg.Seed = 20260306

	s1 := NewBAMSampler(cfg)
	if err := s1.SampleSingleBAM(inputPath, out1); err != nil {
		t.Fatalf("SampleSingleBAM run1: %v", err)
	}

	cfg2 := itConfig()
	cfg2.Ratio = cfg.Ratio
	cfg2.Seed = cfg.Seed
	s2 := NewBAMSampler(cfg2)
	if err := s2.SampleSingleBAM(inputPath, out2); err != nil {
		t.Fatalf("SampleSingleBAM run2: %v", err)
	}

	sig1 := itRecordSignatures(itReadBAM(t, out1))
	sig2 := itRecordSignatures(itReadBAM(t, out2))
	if !reflect.DeepEqual(sig1, sig2) {
		t.Fatalf("same seed should produce identical output signatures")
	}
}

func TestIntegration_DualBAM_DeterministicSeed(t *testing.T) {
	const nPairs = 120
	var r1Recs, r2Recs []*bamnative.Record
	for i := 0; i < nPairs; i++ {
		name := fmt.Sprintf("read%04d", i)
		r1Recs = append(r1Recs, itMakeRecord(name, itFlagR1))
		r2Recs = append(r2Recs, itMakeRecord(name, itFlagR2))
	}

	r1Path := itWriteBAM(t, r1Recs)
	r2Path := itWriteBAM(t, r2Recs)
	prefix1 := filepath.Join(t.TempDir(), "run1")
	prefix2 := filepath.Join(t.TempDir(), "run2")

	cfg := itConfig()
	cfg.Ratio = 0.29
	cfg.Seed = 20260306

	s1 := NewBAMSampler(cfg)
	if err := s1.SampleDualBAMSorted(r1Path, r2Path, prefix1); err != nil {
		t.Fatalf("SampleDualBAMSorted run1: %v", err)
	}

	cfg2 := itConfig()
	cfg2.Ratio = cfg.Ratio
	cfg2.Seed = cfg.Seed
	s2 := NewBAMSampler(cfg2)
	if err := s2.SampleDualBAMSorted(r1Path, r2Path, prefix2); err != nil {
		t.Fatalf("SampleDualBAMSorted run2: %v", err)
	}

	r1Sig1 := itRecordSignatures(itReadBAM(t, prefix1+"_R1.bam"))
	r1Sig2 := itRecordSignatures(itReadBAM(t, prefix2+"_R1.bam"))
	r2Sig1 := itRecordSignatures(itReadBAM(t, prefix1+"_R2.bam"))
	r2Sig2 := itRecordSignatures(itReadBAM(t, prefix2+"_R2.bam"))

	if !reflect.DeepEqual(r1Sig1, r1Sig2) || !reflect.DeepEqual(r2Sig1, r2Sig2) {
		t.Fatalf("same seed should produce identical dual-mode outputs")
	}
}

// itMakeRecordAt builds a BAM record at a specific genomic position.
func itMakeRecordAt(name string, flags uint16, refID int32, pos int32) *bamnative.Record {
	rec := itMakeRecord(name, flags)
	rec.RefID = refID
	rec.Pos = pos
	return rec
}

// ============================================================================
// Sort-path integration tests
// ============================================================================

// TestIntegration_SingleBAM_CoordSort verifies that -sort coord produces output
// where every individual record (R1 and R2) is in strictly ascending coordinate order.
func TestIntegration_SingleBAM_CoordSort(t *testing.T) {
	// Deliberately create pairs where R2 positions are NOT in the same order as R1 positions.
	// pair "read1": R1 at pos 300, R2 at pos 100
	// pair "read2": R1 at pos 100, R2 at pos 300
	// pair "read3": R1 at pos 200, R2 at pos 200
	// A pair-level sort by R1 coord would yield: read2(100),read1(300),read3(200) — wrong.
	// Per-record sort should yield global positions: 100,100,200,200,300,300.
	records := []*bamnative.Record{
		itMakeRecordAt("read1", itFlagR1, 0, 300),
		itMakeRecordAt("read1", itFlagR2, 0, 100),
		itMakeRecordAt("read2", itFlagR1, 0, 100),
		itMakeRecordAt("read2", itFlagR2, 0, 300),
		itMakeRecordAt("read3", itFlagR1, 0, 200),
		itMakeRecordAt("read3", itFlagR2, 0, 200),
	}

	inputPath := itWriteBAM(t, records)
	outputPath := filepath.Join(t.TempDir(), "out.bam")

	cfg := itConfig()
	cfg.Ratio = 1.0
	cfg.SortOrder = SortByCoord
	s := NewBAMSampler(cfg)
	if err := s.SampleSingleBAM(inputPath, outputPath); err != nil {
		t.Fatalf("SampleSingleBAM: %v", err)
	}

	outRecs := itReadBAM(t, outputPath)
	if len(outRecs) != 6 {
		t.Fatalf("expected 6 output records, got %d", len(outRecs))
	}

	// Every consecutive pair of records must be in non-decreasing coordinate order.
	for i := 1; i < len(outRecs); i++ {
		prev, cur := outRecs[i-1], outRecs[i]
		if cur.RefID < prev.RefID || (cur.RefID == prev.RefID && cur.Pos < prev.Pos) {
			t.Errorf("record[%d] (RefID=%d Pos=%d) is before record[%d] (RefID=%d Pos=%d); not in coord order",
				i, cur.RefID, cur.Pos, i-1, prev.RefID, prev.Pos)
		}
	}
}

// TestIntegration_SingleBAM_NameSort verifies that -sort name produces output where
// records with the same name are grouped and groups appear in lexicographic order.
func TestIntegration_SingleBAM_NameSort(t *testing.T) {
	// Insert records out of alphabetical order.
	records := []*bamnative.Record{
		itMakeRecord("charlie", itFlagR1),
		itMakeRecord("charlie", itFlagR2),
		itMakeRecord("alpha", itFlagR1),
		itMakeRecord("alpha", itFlagR2),
		itMakeRecord("bravo", itFlagR1),
		itMakeRecord("bravo", itFlagR2),
	}

	inputPath := itWriteBAM(t, records)
	outputPath := filepath.Join(t.TempDir(), "out.bam")

	cfg := itConfig()
	cfg.Ratio = 1.0
	cfg.SortOrder = SortByName
	s := NewBAMSampler(cfg)
	if err := s.SampleSingleBAM(inputPath, outputPath); err != nil {
		t.Fatalf("SampleSingleBAM: %v", err)
	}

	outRecs := itReadBAM(t, outputPath)
	if len(outRecs) != 6 {
		t.Fatalf("expected 6 output records, got %d", len(outRecs))
	}

	// Consecutive records with the same name must be adjacent, and group names must be sorted.
	var prevGroupName string
	for i := 0; i < len(outRecs); i += 2 {
		r1, r2 := outRecs[i], outRecs[i+1]
		if r1.Name != r2.Name {
			t.Errorf("record pair at [%d,%d]: names %q and %q should match", i, i+1, r1.Name, r2.Name)
		}
		if prevGroupName != "" && r1.Name < prevGroupName {
			t.Errorf("group %q appears before %q; not in name order", r1.Name, prevGroupName)
		}
		prevGroupName = r1.Name
	}
}

// TestIntegration_DualBAM_CoordSort verifies that dual-mode -sort coord independently sorts
// each output file by its own records' coordinates (not by the mate's coordinates).
func TestIntegration_DualBAM_CoordSort(t *testing.T) {
	// R1 positions: read1=300, read2=100, read3=200
	// R2 positions: read1=100, read2=300, read3=200
	// Sorting R1 file by R1 coord: read2(100), read3(200), read1(300) ✓
	// Sorting R2 file by R2 coord: read1(100), read3(200), read2(300) ✓
	// If R2 were naively ordered by R1 coord it would be: read1(100→R2pos=100)... actually
	// let's use more distinct values to make the test unambiguous.
	//
	// R1 positions: read1=300, read2=100
	// R2 positions: read1=50,  read2=400
	// Pair-level sort by R1: [read2(100), read1(300)]
	//   → R1 file: 100, 300 (correct)
	//   → R2 file: 400, 50  (WRONG – not sorted by R2's own position)
	// Per-record sort for R2: [50, 400] (correct)
	r1Recs := []*bamnative.Record{
		itMakeRecordAt("read1", itFlagR1, 0, 300),
		itMakeRecordAt("read2", itFlagR1, 0, 100),
	}
	r2Recs := []*bamnative.Record{
		itMakeRecordAt("read1", itFlagR2, 0, 50),
		itMakeRecordAt("read2", itFlagR2, 0, 400),
	}

	r1Path := itWriteBAM(t, r1Recs)
	r2Path := itWriteBAM(t, r2Recs)
	outDir := t.TempDir()
	prefix := filepath.Join(outDir, "out")

	cfg := itConfig()
	cfg.Ratio = 1.0
	cfg.SortOrder = SortByCoord
	s := NewBAMSampler(cfg)
	if err := s.SampleDualBAMSorted(r1Path, r2Path, prefix); err != nil {
		t.Fatalf("SampleDualBAMSorted: %v", err)
	}

	r1Out := itReadBAM(t, prefix+"_R1.bam")
	r2Out := itReadBAM(t, prefix+"_R2.bam")

	checkCoordOrder := func(label string, recs []*bamnative.Record) {
		for i := 1; i < len(recs); i++ {
			prev, cur := recs[i-1], recs[i]
			if cur.RefID < prev.RefID || (cur.RefID == prev.RefID && cur.Pos < prev.Pos) {
				t.Errorf("%s: record[%d] (Pos=%d) is before record[%d] (Pos=%d); not in coord order",
					label, i, cur.Pos, i-1, prev.Pos)
			}
		}
	}
	checkCoordOrder("R1 file", r1Out)
	checkCoordOrder("R2 file", r2Out)

	// Verify R2 file is sorted by R2's own positions, not R1's.
	if len(r2Out) == 2 && r2Out[0].Pos != 50 {
		t.Errorf("R2 file first record Pos = %d; want 50 (sorted by R2's own coord)", r2Out[0].Pos)
	}
}

// TestIntegration_DualBAM_FilteredByProperPair verifies that dual-mode -proper-pairs-only
// correctly counts pairs filtered by the proper-pair flag.
func TestIntegration_DualBAM_FilteredByProperPair(t *testing.T) {
	// Flags for a proper pair: Paired + ProperPair + neither unmapped nor mate-unmapped.
	properFlagR1 := uint16(bamnative.FlagPaired | bamnative.FlagProperPair | bamnative.FlagFirstInPair)
	properFlagR2 := uint16(bamnative.FlagPaired | bamnative.FlagProperPair | bamnative.FlagSecondInPair)
	// Non-proper pair: both sides are paired but without ProperPair.
	nonProperFlagR1 := uint16(bamnative.FlagPaired | bamnative.FlagFirstInPair)
	nonProperFlagR2 := uint16(bamnative.FlagPaired | bamnative.FlagSecondInPair)

	r1Recs := []*bamnative.Record{
		itMakeRecord("proper", properFlagR1),       // kept
		itMakeRecord("nonproper", nonProperFlagR1), // filtered out as a pair
	}
	r2Recs := []*bamnative.Record{
		itMakeRecord("proper", properFlagR2),
		itMakeRecord("nonproper", nonProperFlagR2),
	}

	r1Path := itWriteBAM(t, r1Recs)
	r2Path := itWriteBAM(t, r2Recs)
	outDir := t.TempDir()
	prefix := filepath.Join(outDir, "out")

	cfg := itConfig()
	cfg.Ratio = 1.0
	cfg.ProperPairsOnly = true
	s := NewBAMSampler(cfg)
	if err := s.SampleDualBAMSorted(r1Path, r2Path, prefix); err != nil {
		t.Fatalf("SampleDualBAMSorted: %v", err)
	}

	stats := s.GetStats()
	if stats.FilteredByProperPair != 1 {
		t.Errorf("FilteredByProperPair = %d; want 1", stats.FilteredByProperPair)
	}
	if stats.TotalPairedOutput != 1 {
		t.Errorf("TotalPairedOutput = %d; want 1 (only 'proper' pair survives)", stats.TotalPairedOutput)
	}
}
func TestIntegration_DualBAM_OrphanReads(t *testing.T) {
	r1Recs := []*bamnative.Record{
		itMakeRecord("paired", itFlagR1),
		itMakeRecord("orphan", itFlagR1), // no R2 counterpart
	}
	r2Recs := []*bamnative.Record{
		itMakeRecord("paired", itFlagR2),
	}

	r1Path := itWriteBAM(t, r1Recs)
	r2Path := itWriteBAM(t, r2Recs)
	outDir := t.TempDir()
	prefix := filepath.Join(outDir, "out")

	cfg := itConfig()
	cfg.Ratio = 1.0
	s := NewBAMSampler(cfg)
	if err := s.SampleDualBAMSorted(r1Path, r2Path, prefix); err != nil {
		t.Fatalf("SampleDualBAMSorted: %v", err)
	}

	stats := s.GetStats()
	if stats.OrphanReads != 1 {
		t.Errorf("OrphanReads = %d; want 1", stats.OrphanReads)
	}
	if stats.TotalPairedOutput != 1 {
		t.Errorf("TotalPairedOutput = %d; want 1", stats.TotalPairedOutput)
	}
}

func TestIntegration_DualBAM_InputRoleValidation(t *testing.T) {
	r1Wrong := []*bamnative.Record{
		itMakeRecord("read1", itFlagR2),
		itMakeRecord("read2", itFlagR2),
		itMakeRecord("read3", itFlagR2),
	}
	r2Wrong := []*bamnative.Record{
		itMakeRecord("read1", itFlagR1),
		itMakeRecord("read2", itFlagR1),
		itMakeRecord("read3", itFlagR1),
	}

	r1Path := itWriteBAM(t, r1Wrong)
	r2Path := itWriteBAM(t, r2Wrong)
	prefix := filepath.Join(t.TempDir(), "out")

	cfg := itConfig()
	cfg.Ratio = 1.0
	s := NewBAMSampler(cfg)
	err := s.SampleDualBAMSorted(r1Path, r2Path, prefix)
	if err == nil {
		t.Fatalf("expected validation error for swapped R1/R2 inputs")
	}
	if !strings.Contains(err.Error(), "appears to contain more R2 than R1") {
		t.Fatalf("unexpected error: %v", err)
	}
}

// ============================================================================
// Boundary condition tests [L-1]
// ============================================================================

// itReadBAMHeader opens a BAM file and returns only its header.
func itReadBAMHeader(t *testing.T, path string) *bamnative.Header {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("itReadBAMHeader: Open %s: %v", path, err)
	}
	defer f.Close()
	reader, err := bamnative.NewReader(f)
	if err != nil {
		t.Fatalf("itReadBAMHeader: NewReader: %v", err)
	}
	return reader.Header()
}

// TestIntegration_SingleBAM_EmptyInput verifies that an empty BAM produces an
// empty output file without error.
func TestIntegration_SingleBAM_EmptyInput(t *testing.T) {
	inputPath := itWriteBAM(t, nil)
	outputPath := filepath.Join(t.TempDir(), "out.bam")

	cfg := itConfig()
	cfg.Ratio = 0.5
	s := NewBAMSampler(cfg)
	if err := s.SampleSingleBAM(inputPath, outputPath); err != nil {
		t.Fatalf("SampleSingleBAM on empty input: %v", err)
	}

	outRecs := itReadBAM(t, outputPath)
	if len(outRecs) != 0 {
		t.Errorf("expected 0 output records for empty input, got %d", len(outRecs))
	}
	stats := s.GetStats()
	if stats.TotalOutputRecords != 0 {
		t.Errorf("TotalOutputRecords = %d; want 0", stats.TotalOutputRecords)
	}
}

// TestIntegration_SingleBAM_AbsoluteCountExceedsPairs verifies that when
// AbsoluteCount is greater than the actual number of pairs, all pairs are
// returned without panic.
func TestIntegration_SingleBAM_AbsoluteCountExceedsPairs(t *testing.T) {
	const nPairs = 5
	var records []*bamnative.Record
	for i := 0; i < nPairs; i++ {
		name := fmt.Sprintf("read%04d", i)
		records = append(records, itMakeRecord(name, itFlagR1))
		records = append(records, itMakeRecord(name, itFlagR2))
	}

	inputPath := itWriteBAM(t, records)
	outputPath := filepath.Join(t.TempDir(), "out.bam")

	cfg := itConfig()
	cfg.Mode = ModeAbsolute
	cfg.AbsoluteCount = 1000 // far more than nPairs
	s := NewBAMSampler(cfg)
	if err := s.SampleSingleBAM(inputPath, outputPath); err != nil {
		t.Fatalf("SampleSingleBAM: %v", err)
	}

	stats := s.GetStats()
	if stats.TotalPairedOutput != nPairs {
		t.Errorf("TotalPairedOutput = %d; want %d (all pairs)", stats.TotalPairedOutput, nPairs)
	}
	outRecs := itReadBAM(t, outputPath)
	if len(outRecs) != nPairs*2 {
		t.Errorf("output record count = %d; want %d", len(outRecs), nPairs*2)
	}
}

// TestIntegration_SingleBAM_SortNoneHeaderUnknown verifies that when SortOrder
// is SortNone, the output BAM header has SortOrder = "unknown".
func TestIntegration_SingleBAM_SortNoneHeaderUnknown(t *testing.T) {
	records := []*bamnative.Record{
		itMakeRecord("r1", itFlagR1),
		itMakeRecord("r1", itFlagR2),
	}

	inputPath := itWriteBAM(t, records)
	outputPath := filepath.Join(t.TempDir(), "out.bam")

	cfg := itConfig()
	cfg.Ratio = 1.0
	cfg.SortOrder = SortNone
	s := NewBAMSampler(cfg)
	if err := s.SampleSingleBAM(inputPath, outputPath); err != nil {
		t.Fatalf("SampleSingleBAM: %v", err)
	}

	hdr := itReadBAMHeader(t, outputPath)
	if hdr.SortOrder != "unknown" {
		t.Errorf("SortOrder in output header = %q; want %q", hdr.SortOrder, "unknown")
	}
}

// TestIntegration_SingleBAM_UnmarkedReads verifies that reads with neither the
// R1 nor the R2 flag are silently dropped and counted in UnmarkedReads.
func TestIntegration_SingleBAM_UnmarkedReads(t *testing.T) {
	// "unmarked" has FlagPaired but neither FlagFirstInPair nor FlagSecondInPair.
	flagUnmarked := uint16(bamnative.FlagPaired)
	records := []*bamnative.Record{
		itMakeRecord("paired", itFlagR1),
		itMakeRecord("paired", itFlagR2),
		itMakeRecord("unmarked", flagUnmarked), // should be dropped
	}

	inputPath := itWriteBAM(t, records)
	outputPath := filepath.Join(t.TempDir(), "out.bam")

	cfg := itConfig()
	cfg.Ratio = 1.0
	s := NewBAMSampler(cfg)
	if err := s.SampleSingleBAM(inputPath, outputPath); err != nil {
		t.Fatalf("SampleSingleBAM: %v", err)
	}

	stats := s.GetStats()
	if stats.UnmarkedReads != 1 {
		t.Errorf("UnmarkedReads = %d; want 1", stats.UnmarkedReads)
	}
	if stats.TotalPairedOutput != 1 {
		t.Errorf("TotalPairedOutput = %d; want 1", stats.TotalPairedOutput)
	}
	outRecs := itReadBAM(t, outputPath)
	if len(outRecs) != 2 {
		t.Errorf("output record count = %d; want 2 (only the paired reads)", len(outRecs))
	}
}
