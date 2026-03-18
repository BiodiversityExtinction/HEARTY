# HEARTY

HEARTY (Heterozygosity Error and Ancient-damage Robust Thresholding sYstem) is a lightweight toolkit for profiling heterozygosity from BAM files with ANGSD, assessing the impact of sequencing error and ancient DNA damage, choosing empirical heterozygosity thresholds, and carrying those choices into downstream ROH and PSMC workflows.

The main idea is:
- build one reusable basecall table from ANGSD counts
- always inspect the exploratory `0.05` threshold
- add stricter thresholds when needed
- regenerate downstream summaries with different depth filters without rerunning ANGSD

## Workflow Overview

1. Run `HEARTY.py` to generate ANGSD counts, a reusable basecall table, the exploratory minor allele frequency table, and threshold-specific 100 kb window summaries.
2. Use `Plot_MinorFreq.R` to inspect the exploratory `0.05` minor allele frequency spectrum and choose a stricter heterozygosity threshold.
3. Use `ROH_Plot_Tool.R` to visualize the threshold-specific window summaries and inspect ROH-like patterns.
4. Use `PSMC_Tool.sh` to prepare masked PSMC input from HEARTY outputs and run PSMC, including optional bootstrap replicates.

## Features

- Heterozygosity calling with a fixed exploratory `0.05` threshold plus optional extra thresholds
- One reusable basecall table per sample with threshold-specific call columns
- 100 kb window heterozygosity summaries for all sites and transversions only
- Consolidated minor allele frequency summaries for threshold exploration
- Standalone plotting tools for minor allele frequency and ROH/window summaries
- Standalone PSMC helper that uses HEARTY threshold-specific het calls for masking and runs PSMC
- Batch processing for both BAM inputs and minorfreq/ROH plotting

## Installation

```bash
git clone https://github.com/BiodiversityExtinction/HEARTY.git
cd HEARTY
chmod +x HEARTY.py
chmod +x PSMC_Tool.sh
```

## Requirements

Core workflow:
- Python 3.x
- [ANGSD](http://www.popgen.dk/angsd/) (tested with v0.935)
- Standard Unix tools: `bash`, `awk`, `grep`, `paste`
- `pigz` optional but recommended for multithreaded compression/decompression

Plotting:
- R
- Packages used by `ROH_Plot_Tool.R`: `tidyverse`, `ggplot2`

PSMC helper:
- `bcftools`
- `vcfutils.pl`
- `seqtk`
- `bedtools`
- `fq2psmcfa`
- `splitfa`
- `psmc`

## Main Tool

`HEARTY.py` is the main caller.

### Quick start

Single BAM:

```bash
python3 HEARTY.py -b sample.bam -o sample01
```

Single BAM with an extra stricter threshold:

```bash
python3 HEARTY.py -b sample.bam -o sample01 -t 0.25
```

Region-restricted run:

```bash
python3 HEARTY.py -b sample.bam -o sample01 -r "chr1:1-50000000"
```

Reuse an existing basecall file and only rebuild downstream summaries:

```bash
python3 HEARTY.py -b sample.bam -o sample01 -d . --reuse-existing -m 20 -M 60 -t 0.25
```

### Usage

```text
python3 HEARTY.py [OPTIONS]

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

This is the format expected by `PSMC_Tool.sh`.

## Minor Frequency Plotting

`Plot_MinorFreq.R` plots the consolidated minorfreq table produced by HEARTY.

### Single-sample usage

```bash
Rscript Plot_MinorFreq.R \
  --minorfreq-file results/sample01.mincov10.minorfreq.txt \
  --out results/plots/sample01 \
  --format pdf
```

With capped y-axis based on the tail:

```bash
Rscript Plot_MinorFreq.R \
  --minorfreq-file results/sample01.mincov10.minorfreq.txt \
  --out results/plots/sample01 \
  --format pdf \
  --y-cap-after-x 0.25
```

### Batch usage

```bash
Rscript Plot_MinorFreq.R \
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

`ROH_Plot_Tool.R` takes HEARTY `*.windows.txt.gz` files and turns them into visual summaries of low-heterozygosity regions.

The input files are the 100 kb window summaries written by `HEARTY.py`, for example:

```text
scaffold_1 0-99999: 28741 sites, 3 sites with 'HET', Proportion: 0.0001044
```

For each window, the script reads:
- the scaffold name
- the window coordinates
- the number of callable sites in that window
- the number of HET sites in that window
- the window heterozygosity proportion

The plotting script then:
- colors windows by heterozygosity proportion
- flags windows at or below a chosen ROH-like threshold
- optionally “bridges” single-bin gaps, so an isolated non-ROH bin between two ROH bins can still be treated as part of a continuous ROH block
- calculates segment-based FROH summaries from the raw and bridged ROH calls

So this script is not calling ROH directly from reads or genotype likelihoods. It is interpreting the window-level heterozygosity profiles produced by HEARTY.

### Basic usage

Show help:

```bash
Rscript ROH_Plot_Tool.R --help
```

### Single-sample usage

```bash
Rscript ROH_Plot_Tool.R \
  --mode single \
  --input results/sample01_het0.25.mincov10.windows.txt.gz \
  --output-pdf results/plots/sample01_roh.pdf
```

Full example with explicit tuning parameters:

```bash
Rscript ROH_Plot_Tool.R \
  --mode single \
  --input sample.ROH.txt \
  --output-pdf sample.ROH.pdf \
  --min-length-mb 20 \
  --max-scaffolds 40 \
  --fill-cap 0.0005 \
  --bin-size-bp 500000 \
  --fallback-resolution 100 \
  --roh-threshold 0.0001
```

Single-sample mode writes one PDF for one sample and prints summary statistics to the terminal, including:
- number of windows parsed
- number of scaffolds plotted
- raw and bridged ROH bin counts
- FROH values for segments longer than 100 kb, 1 Mb, and 5 Mb

### Batch usage

```bash
Rscript ROH_Plot_Tool.R \
  --mode batch \
  --sample-list roh_samples.tsv \
  --summary-out roh_summary.tsv \
  --combined-pdf roh_combined.pdf
```

Full batch example:

```bash
Rscript ROH_Plot_Tool.R \
  --mode batch \
  --sample-list Sample_table.txt \
  --summary-out ROH_batch_summary.tsv \
  --combined-pdf ROH_batch_combined.pdf \
  --min-length-mb 20 \
  --max-scaffolds 40 \
  --fill-cap 0.0005 \
  --bin-size-bp 500000 \
  --fallback-resolution 100 \
  --roh-threshold 0.0001
```

Batch list format:

```text
sampleA<TAB>/path/sampleA_het0.25.mincov10.windows.txt.gz
sampleB<TAB>/path/sampleB_het0.25.mincov10.windows.txt.gz
```

Headered sample list files are also accepted, for example:

```text
sample	file
SampleA	/path/SampleA.ROH.txt
SampleB	/path/SampleB.ROH.txt
```

Batch mode writes:
- a summary TSV with FROH and plotting metadata for each sample
- a combined PDF with multiple samples arranged together for comparison

### Required inputs

- `--mode`
  Must be either `single` or `batch`.
- Single mode:
  Requires `--input` for one window-summary file and `--output-pdf` for one output PDF.
- Batch mode:
  Requires `--sample-list`, `--summary-out`, and `--combined-pdf`.

### Interpreting the plot

- Window fill color:
  Lower heterozygosity windows are shown toward the low end of the color scale and higher heterozygosity windows toward the high end.
- Black ticks:
  Windows whose heterozygosity proportion is less than or equal to `--roh-threshold`.
- Orange ticks:
  “Bridged” ROH calls, where a single non-ROH window between two ROH windows is filled in to create a more continuous segment.

This gives you two related views:
- a strict raw ROH-like signal
- a slightly smoothed version that is less sensitive to one noisy window interrupting a longer low-heterozygosity tract

### Parameter meanings

- `--roh-threshold`
  The heterozygosity proportion cutoff used to mark a window as ROH-like. Lower values are stricter. This is the most important biological setting.
  A window is treated as raw ROH if heterozygosity is less than or equal to this threshold.
- `--fill-cap`
  Caps the displayed color scale at a maximum heterozygosity value so very high windows do not dominate the gradient. This changes the visualization, not the ROH calling itself.
- `--min-length-mb`
  Minimum scaffold length to include in the plot. Short scaffolds are excluded from the visual output.
- `--max-scaffolds`
  Maximum number of scaffolds to plot per sample, ranked by scaffold length.
- `--bin-size-bp`
  Bin size used for plotting aggregation. This controls visual smoothing in the PDF and does not change the original HEARTY windows file.
- `--fallback-resolution`
  PDF fallback resolution used by `ggsave`; mainly relevant for rendering quality and file export.

### FROH in the outputs

The script reports FROH for both raw and bridged ROH using segment cutoffs:
- `>100 kb`
- `>1 Mb`
- `>5 Mb`

The genome size denominator is:
- `number_of_windows * 100,000 bp`

### Recommended defaults

- Start by matching the threshold used to make the windows file, for example a `0.25` HEARTY run with a conservative `--roh-threshold`.
- Use the transversion-only windows file if you want a damage-robust view.
- For most runs, start with:
  `--min-length-mb 20`
  `--max-scaffolds 40`
  `--fill-cap 0.0005`
  `--bin-size-bp 500000`
  `--fallback-resolution 100`
  `--roh-threshold 0.0001`
- If the PDF is too large:
  increase `--bin-size-bp` such as `1000000`
  reduce `--max-scaffolds` such as `20`
  lower `--fallback-resolution` such as `80`

## PSMC Helper

`PSMC_Tool.sh` prepares masked PSMC input using a HEARTY basecall table and a chosen threshold column, and then runs PSMC on the masked result.

### Usage

```bash
bash PSMC_Tool.sh \
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
- `Plot_MinorFreq.R`, `ROH_Plot_Tool.R`, and `PSMC_Tool.sh` are standalone helpers and can be run independently of the main HEARTY workflow once the relevant HEARTY outputs exist.

## Repository Layout

- `HEARTY.py` - main HEARTY caller
- `Plot_MinorFreq.R` - minor allele frequency overlay plotter
- `ROH_Plot_Tool.R` - ROH/window plotting utility
- `PSMC_Tool.sh` - PSMC preparation and run helper

## Contact

micwe@dtu.dk
