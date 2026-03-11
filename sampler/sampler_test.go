package sampler

import (
	"testing"
)

// TestReservoirSamplesAll verifies that when total items ≤ capacity, every item is kept.
func TestReservoirSamplesAll(t *testing.T) {
	r := NewReservoir(10, 42)
	for i := 0; i < 5; i++ {
		r.Add(i)
	}
	got := r.GetRecords()
	if len(got) != 5 {
		t.Errorf("expected 5 records, got %d", len(got))
	}
}

// TestReservoirExactCapacity verifies that when items == capacity, all are kept.
func TestReservoirExactCapacity(t *testing.T) {
	r := NewReservoir(5, 42)
	for i := 0; i < 5; i++ {
		r.Add(i)
	}
	got := r.GetRecords()
	if len(got) != 5 {
		t.Errorf("expected 5 records, got %d", len(got))
	}
}

// TestReservoirRespectsCapacity verifies that the reservoir never exceeds its capacity.
func TestReservoirRespectsCapacity(t *testing.T) {
	const cap int64 = 10
	r := NewReservoir(cap, 42)
	for i := 0; i < 1000; i++ {
		r.Add(i)
	}
	got := r.GetRecords()
	if int64(len(got)) != cap {
		t.Errorf("expected %d records, got %d", cap, len(got))
	}
}

// TestReservoirSizeTracks verifies that Size() reflects all items ever added, not just kept ones.
func TestReservoirSizeTracks(t *testing.T) {
	r := NewReservoir(5, 42)
	for i := 0; i < 100; i++ {
		r.Add(i)
	}
	if r.Size() != 100 {
		t.Errorf("expected Size()=100, got %d", r.Size())
	}
}

// TestReservoirDeterministic verifies that the same seed produces the same sample.
func TestReservoirDeterministic(t *testing.T) {
	const n = 1000
	const cap int64 = 50

	run := func() []interface{} {
		r := NewReservoir(cap, 12345)
		for i := 0; i < n; i++ {
			r.Add(i)
		}
		return r.GetRecords()
	}

	a := run()
	b := run()

	if len(a) != len(b) {
		t.Fatalf("length mismatch: %d vs %d", len(a), len(b))
	}
	for i := range a {
		if a[i] != b[i] {
			t.Errorf("index %d: %v != %v", i, a[i], b[i])
		}
	}
}

// TestReservoirCapacityAccessor verifies Capacity() returns the configured value.
func TestReservoirCapacityAccessor(t *testing.T) {
	r := NewReservoir(42, 1)
	if r.Capacity() != 42 {
		t.Errorf("expected capacity 42, got %d", r.Capacity())
	}
}

// TestReservoirGetRecordsIsCopy verifies that GetRecords returns a copy, not the internal slice.
func TestReservoirGetRecordsIsCopy(t *testing.T) {
	r := NewReservoir(5, 1)
	r.Add("a")
	r.Add("b")

	got := r.GetRecords()
	got[0] = "mutated"

	got2 := r.GetRecords()
	if got2[0] == "mutated" {
		t.Error("GetRecords should return a copy, not a reference to internal slice")
	}
}

// TestNewReservoirAutoSeed verifies that seed=0 does not panic (falls back to time-based seed).
func TestNewReservoirAutoSeed(t *testing.T) {
	defer func() {
		if r := recover(); r != nil {
			t.Errorf("NewReservoir panicked with seed=0: %v", r)
		}
	}()
	r := NewReservoir(10, 0)
	r.Add("x")
	if r.Size() != 1 {
		t.Errorf("expected Size()=1, got %d", r.Size())
	}
}
