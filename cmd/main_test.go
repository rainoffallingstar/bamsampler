package main

import "testing"

func TestValidateArgsSamplingModeMutualExclusion(t *testing.T) {
	tests := []struct {
		name    string
		args    *CLIArgs
		wantErr bool
	}{
		{
			name: "ratio only",
			args: &CLIArgs{
				Ratio:      0.1,
				InputMode:  "single",
				SortOrder:  "none",
				InputFiles: []string{"in.bam", "out.bam"},
			},
			wantErr: false,
		},
		{
			name: "count only",
			args: &CLIArgs{
				AbsoluteCount: 100,
				InputMode:     "single",
				SortOrder:     "none",
				InputFiles:    []string{"in.bam", "out.bam"},
			},
			wantErr: false,
		},
		{
			name: "missing ratio and count",
			args: &CLIArgs{
				InputMode:  "single",
				SortOrder:  "none",
				InputFiles: []string{"in.bam", "out.bam"},
			},
			wantErr: true,
		},
		{
			name: "ratio and count both set",
			args: &CLIArgs{
				Ratio:         0.1,
				AbsoluteCount: 100,
				InputMode:     "single",
				SortOrder:     "none",
				InputFiles:    []string{"in.bam", "out.bam"},
			},
			wantErr: true,
		},
	}

	for _, tc := range tests {
		err := validateArgs(tc.args)
		if tc.wantErr && err == nil {
			t.Fatalf("%s: expected error, got nil", tc.name)
		}
		if !tc.wantErr && err != nil {
			t.Fatalf("%s: unexpected error: %v", tc.name, err)
		}
	}
}

func TestValidatePositionalArgs(t *testing.T) {
	tests := []struct {
		name    string
		args    *CLIArgs
		wantErr bool
	}{
		{
			name: "single mode needs two args",
			args: &CLIArgs{
				InputMode:  "single",
				InputFiles: []string{"in.bam"},
			},
			wantErr: true,
		},
		{
			name: "single mode ok",
			args: &CLIArgs{
				InputMode:  "single",
				InputFiles: []string{"in.bam", "out.bam"},
			},
			wantErr: false,
		},
		{
			name: "dual mode needs three args",
			args: &CLIArgs{
				InputMode:  "dual",
				InputFiles: []string{"r1.bam", "r2.bam"},
			},
			wantErr: true,
		},
		{
			name: "dual mode ok",
			args: &CLIArgs{
				InputMode:  "dual",
				InputFiles: []string{"r1.bam", "r2.bam", "prefix"},
			},
			wantErr: false,
		},
	}

	for _, tc := range tests {
		err := validatePositionalArgs(tc.args)
		if tc.wantErr && err == nil {
			t.Fatalf("%s: expected error, got nil", tc.name)
		}
		if !tc.wantErr && err != nil {
			t.Fatalf("%s: unexpected error: %v", tc.name, err)
		}
	}
}

func TestValidateArgsIndexRequiresCoordSort(t *testing.T) {
	tests := []struct {
		name    string
		args    *CLIArgs
		wantErr bool
	}{
		{
			name: "index without coord sort errors",
			args: &CLIArgs{
				Ratio:       0.1,
				InputMode:   "single",
				SortOrder:   "none",
				CreateIndex: true,
				InputFiles:  []string{"in.bam", "out.bam"},
			},
			wantErr: true,
		},
		{
			name: "index with coord sort passes",
			args: &CLIArgs{
				Ratio:       0.1,
				InputMode:   "single",
				SortOrder:   "coord",
				CreateIndex: true,
				InputFiles:  []string{"in.bam", "out.bam"},
			},
			wantErr: false,
		},
	}

	for _, tc := range tests {
		err := validateArgs(tc.args)
		if tc.wantErr && err == nil {
			t.Fatalf("%s: expected error, got nil", tc.name)
		}
		if !tc.wantErr && err != nil {
			t.Fatalf("%s: unexpected error: %v", tc.name, err)
		}
	}
}
