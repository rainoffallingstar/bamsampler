# bamsampler

`bamsampler` is a Go CLI for random sampling of mapped BAM data with pair integrity guarantees.

## Features

- Single-BAM mode: sample paired reads from one BAM while preserving R1/R2 pairing.
- Dual-BAM mode: sample matching read names from separate `R1.bam` and `R2.bam`.
- Two sampling modes:
  - `-ratio` (fraction of read pairs)
  - `-count` (absolute number of read pairs)
- Optional filters/sorting:
  - `-proper-pairs-only`
  - `-sort name|coord|none`
  - `-index` (BAI generation; requires `-sort coord`)
- Seeded reproducibility:
  - with fixed input and same `-seed`, outputs are deterministic.

## Usage

```bash
# Single-BAM ratio sampling
bamsampler -ratio 0.1 -seed 42 input.bam output.bam

# Single-BAM absolute pair count
bamsampler -count 1000000 input.bam output.bam

# Dual-BAM ratio sampling
bamsampler -mode dual -ratio 0.1 -seed 42 R1.bam R2.bam out_prefix
```

## Production notes

- In dual mode, `R1.bam` should primarily contain first-in-pair reads and `R2.bam` should primarily contain second-in-pair reads. Obvious swaps are rejected.
- `-index` requires coordinate-sorted output (`-sort coord`); otherwise the CLI exits with an error.
- Sampling is pair-name based; multi-mappers keep the first primary alignment per side and count dropped extras in stats.
- Orphan pairs are reported as read names without a usable mate after filtering.
