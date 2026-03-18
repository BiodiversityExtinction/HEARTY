# HEARTY

HEARTY (Heterozygosity Error and Ancient-damage Robust Thresholding sYstem) is a lightweight toolkit for profiling heterozygosity from BAM files with ANGSD, assessing the impact of sequencing error and ancient DNA damage, choosing empirical heterozygosity thresholds, and carrying those choices into downstream ROH and PSMC workflows.

The main idea is:
- build one reusable basecall table from ANGSD counts
- always inspect the exploratory `0.05` threshold
- add stricter thresholds when needed
- regenerate downstream summaries with different depth filters without rerunning ANGSD

## Features

- Heterozygosity calling with a fixed exploratory `0.05` threshold plus optional extra thresholds
- One reusable basecall table per sample with threshold-specific call columns
- 100 kb window heterozygosity summaries for all sites and transversions only
- Consolidated minor allele frequency summaries for threshold exploration
- Standalone plotting tools for minor allele frequency and ROH/window summaries
- Standalone PSMC helper that uses HEARTY threshold-specific het calls for masking
- Batch processing for both BAM inputs and minorfreq/ROH plotting

## Installation

```bash
git clone https://github.com/BiodiversityExtinction/HEARTY.git
cd HEARTY
chmod +x hearty.py
chmod +x psmc_tool.sh
```

## Requirements

Core workflow:
- Python 3.x
- [ANGSD](http://www.popgen.dk/angsd/) (tested with v0.935)
- Standard Unix tools: `bash`, `awk`, `grep`, `paste`
- `pigz` optional but recommended for multithreaded compression/decompression

Plotting:
- R
- Packages used by `roh_plot_tool.R`: `tidyverse`, `ggplot2`

PSMC helper:
- `bcftools`
- `vcfutils.pl`
- `seqtk`
- `bedtools`
- `fq2psmcfa`
- `splitfa`
- `psmc`

## Main Tool

`hearty.py` is the main caller.

### Quick start

Single BAM:

```bash
python3 hearty.py -b sample.bam -o sample01
```

Single BAM with an extra stricter threshold:

```bash
python3 hearty.py -b sample.bam -o sample01 -t 0.25
```

Region-restricted run:

```bash
python3 hearty.py -b sample.bam -o sample01 -r "chr1:1-50000000"
```

Reuse an existing basecall file and only rebuild downstream summaries:

```bash
python3 hearty.py -b sample.bam -o sample01 -d . --reuse-existing -m 20 -M 60 -t 0.25
```

### Usage

```text
python3 hearty.py [OPTIONS]

Required arguments:
  -b, --bam FILE             BAM file for analysis
  -l, --bam-list FILE        Two-column file: BAM path and output prefix
  -o, --out-prefix STR       Output prefix for all files (required for single BAM)

Optional arguments:
  -d, --outdir DIR           Output directory [default: results]
  -m, --min-depth INT        Minimum depth for downstream summaries [default: 10]
  -M, --max-depth INT        Maximum depth for downstream summaries
  -t, --threshold FLOAT      Extra het calling threshold for window outputs (can be repeated; 0.05 always runs internally)
  -c, --cores INT            Number of ANGSD/compression threads [default: 8]
  --angsd PATH               ANGSD executable [default: auto]
  -r, --regions STR          Region string for ANGSD -r
  -f, --rf FILE              Regions file for ANGSD -rf
  -n, --dry-run              Show commands without executing them
  -F, --force                Overwrite existing outputs
  --reuse-existing           Reuse existing ANGSD/basecall outputs when available
  -h, --help                 Show help
```

### BAM list format

Each line must contain two columns:

```text
/path/to/sample.bam sample_prefix
```

### Main outputs

For each sample, HEARTY generates:

- `{prefix}.pos.gz` - genomic positions from ANGSD
- `{prefix}.counts.gz` - A/C/G/T counts from ANGSD
- `{prefix}.basecall.txt.gz` - unfiltered basecall table with shared frequencies plus per-threshold call/status columns
- `{prefix}_het{threshold}.mincov{depth}[ _maxcov{depth} ].windows.txt.gz` - 100 kb window heterozygosity summary for explicitly requested `-t` thresholds
- `{prefix}_het{threshold}.mincov{depth}[ _maxcov{depth} ].windows_trv.txt.gz` - 100 kb transversion-only window heterozygosity summary for explicitly requested `-t` thresholds
- `{prefix}.mincov{depth}[ _maxcov{depth} ].minorfreq.txt` - consolidated minor allele frequency table after downstream depth filtering

Important behavior:
- `0.05` always runs internally for the exploratory minor allele frequency output, even if you only specify stricter thresholds
- `--min-depth` and `--max-depth` are applied downstream from the raw basecall table, not during ANGSD counting
- `--reuse-existing` lets you rebuild downstream summaries from an existing basecall file without rerunning ANGSD
- if you rerun with `--reuse-existing` and request a new threshold, HEARTY will rebuild the basecall table from the existing `.pos.gz` and `.counts.gz` files and keep the old threshold columns

### Basecall format

The basecall table is gzip-compressed and tab-delimited. It contains:

- position columns copied from `pos.gz` such as `chr`, `pos`, `totDepth`
- frequency columns: `A`, `C`, `G`, `T`
- for each threshold:
  - `Base_<threshold>`
  - `Status_<threshold>`

Example threshold-specific columns:

```text
Base_0.05   Status_0.05   Base_0.25   Status_0.25
```

This is the format expected by `psmc_tool.sh`.

## Minor Frequency Plotting

`plot_minorfreq.R` plots the consolidated minorfreq table produced by HEARTY.

### Single-sample usage

```bash
Rscript plot_minorfreq.R \
  --minorfreq-file results/sample01.mincov10.minorfreq.txt \
  --out results/plots/sample01 \
  --format pdf
```

With capped y-axis based on the tail:

```bash
Rscript plot_minorfreq.R \
  --minorfreq-file results/sample01.mincov10.minorfreq.txt \
  --out results/plots/sample01 \
  --format pdf \
  --y-cap-after-x 0.25
```

### Batch usage

```bash
Rscript plot_minorfreq.R \
  --sample-list minorfreq_batch.txt \
  --format pdf \
  --y-cap-after-x 0.25
```

Batch list format:

```text
/path/sample1.mincov10.minorfreq.txt   results/plots/sample1
/path/sample2.mincov10.minorfreq.txt   results/plots/sample2
```

Outputs:
- `<output_prefix>_minorfreq_overlay.<format>`
- `<output_prefix>_minorfreq_summary.txt`

The plotting input is the single consolidated HEARTY `*.minorfreq.txt` table with columns:
- `freq`
- `AC`, `AG`, `AT`, `CG`, `CT`, `GT`
- `total`
- `transition`
- `transversion`

## ROH / Window Plotting

`roh_plot_tool.R` takes HEARTY `*.windows.txt.gz` files and plots windowed heterozygosity / derived ROH patterns.

### Single-sample usage

```bash
Rscript roh_plot_tool.R \
  --mode single \
  --input results/sample01_het0.25.mincov10.windows.txt.gz \
  --output-pdf results/plots/sample01_roh.pdf
```

### Batch usage

```bash
Rscript roh_plot_tool.R \
  --mode batch \
  --sample-list roh_samples.tsv \
  --summary-out roh_summary.tsv \
  --combined-pdf roh_combined.pdf
```

Batch list format:

```text
sampleA<TAB>/path/sampleA_het0.25.mincov10.windows.txt.gz
sampleB<TAB>/path/sampleB_het0.25.mincov10.windows.txt.gz
```

Common tuning arguments:
- `--min-length-mb`
- `--max-scaffolds`
- `--fill-cap`
- `--bin-size-bp`
- `--fallback-resolution`
- `--roh-threshold`

## PSMC Helper

`psmc_tool.sh` prepares masked PSMC input using a HEARTY basecall table and a chosen threshold column.

### Usage

```bash
bash psmc_tool.sh \
  --prefix SAMPLE \
  --out-dir /path/to/output_dir \
  --bam /path/to/sample.bam \
  --threads 8 \
  --bed /path/to/regions.bed \
  --reference /path/to/reference.fa \
  --min-coverage-diploid 10 \
  --min-depth-goodhets 10 \
  --het-file /path/to/sample.basecall.txt.gz \
  --threshold 0.25
```

What it does:
- builds a diploid fastq with `bcftools`
- extracts HEARTY-supported het sites from `Status_<threshold>`
- compares those against `seqtk listhet`
- masks unsupported candidate hets
- runs `psmc`
- optionally runs bootstrap replicates

Useful options:
- `--threshold 0.25`
- `--min-coverage-diploid 10`
- `--min-depth-goodhets 10`
- `--max-depth-goodhets 80`
- `--bootstrap-rounds 50`
- `--skip-bootstrap`
- `--skip-diploid`
- `--skip-mask`

## Notes

- `0.05` is always included internally so you can inspect the exploratory minor allele frequency spectrum before choosing stricter thresholds.
- Window outputs are only written for thresholds you explicitly request with `-t`.
- The raw basecall table is intentionally left unfiltered for coverage so you can rerun downstream summaries with different `--min-depth` and `--max-depth` values.
- `plot_minorfreq.R`, `roh_plot_tool.R`, and `psmc_tool.sh` are standalone helpers and can be run independently of the main HEARTY workflow once the relevant HEARTY outputs exist.

## Repository Layout

- `hearty.py` - main HEARTY caller
- `plot_minorfreq.R` - minor allele frequency overlay plotter
- `roh_plot_tool.R` - ROH/window plotting utility
- `psmc_tool.sh` - PSMC preparation and run helper

## Contact

micwe@dtu.dk
