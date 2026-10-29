# HEARTY

HEARTY (Heterozygosity Error and Ancient-damage Robust Thresholding sYstem) is a lightweight toolkit for profiling heterozygosity from BAM files with ANGSD, assessing the impact of sequencing error and ancient DNA damage, choosing empirical heterozygosity thresholds, and carrying those choices into downstream ROH and PSMC workflows.

The main idea is:
- build one reusable basecall table from ANGSD counts
- always inspect the exploratory `0.05` threshold
- add stricter thresholds when needed
- regenerate downstream summaries with different depth filters without rerunning ANGSD
- optionally downsample sites to a fixed depth for cross-coverage comparisons

## Workflow Overview

1. Run `HEARTY.py` to generate ANGSD counts, a reusable basecall table, the exploratory minor allele frequency table, and threshold-specific 100 kb window summaries.
2. Use `Plot_MinorFreq.R` to inspect the exploratory `0.05` minor allele frequency spectrum and choose a stricter heterozygosity threshold.
3. Use `ROH_Window_Summary.R` to define window-based ROH-like regions, calculate FROH summaries, and visualize the threshold-specific window summaries.
4. Use `PSMC_Tool.sh` to prepare masked PSMC input from HEARTY outputs and run PSMC, including optional bootstrap replicates.

## Features

- Heterozygosity calling with a fixed exploratory `0.05` threshold plus optional extra thresholds
- One reusable basecall table per sample with raw counts, frequencies, and threshold-specific call columns
- 100 kb window heterozygosity summaries for all sites and transversions only
- Consolidated minor allele frequency summaries for threshold exploration
- Optional replicate-averaged downsampling for minorfreq and heterozygosity summaries
- Standalone plotting and summary tools for minor allele frequency and ROH/window summaries
- Standalone PSMC helper that uses HEARTY threshold-specific het calls for masking and runs PSMC
- Batch processing for both BAM inputs and minorfreq/ROH window summaries

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
- Packages used by `ROH_Window_Summary.R`: `tidyverse`, `ggplot2`

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

Depth-standardized summaries by downsampling each site to 10 reads, 10 times:

```bash
python3 HEARTY.py -b sample.bam -o sample01 -t 0.25 --downsample-depth 10 --downsample-reps 10
```

### Usage

```text
python3 HEARTY.py [OPTIONS]

Required arguments:
  -b, --bam FILE             BAM file for analysis
  -l, --bam-list FILE        Two-column file: BAM path and output prefix
  --basecall-file FILE       Existing basecall .txt or .txt.gz; downstream summaries only
  -o, --out-prefix STR       Output prefix (required for single BAM or basecall input)

Optional arguments:
  -d, --outdir DIR           Output directory [default: results]
  -m, --min-depth INT        Minimum depth for downstream summaries [default: 10]
  -M, --max-depth INT        Maximum depth for downstream summaries
  -t, --threshold FLOAT      Extra het calling threshold for window outputs (can be repeated; 0.05 always runs internally)
  --downsample-depth INT     Downsample each site to this depth for downstream summaries
  --downsample-reps INT      Number of downsampling replicates [default: 10]
  --downsample-seed INT      Seed used to make per-site downsampling deterministic [default: 1]
  --downstream-workers INT   Worker processes for single-pass summaries [default: 1]
  --downstream-chunk-size INT Input rows per worker task [default: 10000]
  --progress-interval FLOAT  Seconds between progress reports; 0 disables [default: 30]
  -c, --cores INT            Number of ANGSD/compression threads [default: 8]
  --angsd PATH               ANGSD executable [default: auto]
  --ref FILE                Reference FASTA with .fai index; add Ref to the basecall table
  --angsd-params STR         Extra ANGSD parameters, quoted as one string
  -r, --regions STR          Region string for ANGSD -r
  -f, --rf FILE              Regions file for ANGSD -rf
  -n, --dry-run              Show commands without executing them
  -F, --force                Overwrite existing outputs
  --reuse-existing           Reuse existing ANGSD/basecall outputs when available
  -h, --help                 Show help
```

See [Downsampling parameters](#downsampling-parameters) for what the five downsampling and worker options do and how to choose them.

### Reusing a shared basecall file

Use `--basecall-file` as an alternative to `--bam` or `--bam-list` to generate minor-frequency and window summaries directly from an existing table. No BAM, ANGSD executable, or `.pos.gz`/`.counts.gz` files are needed.

```bash
python3 HEARTY.py \
  --basecall-file /path/to/SW_07_newtool.basecall.txt.gz \
  -o SW_07_depth20_60 \
  -d results/depth_comparison \
  -m 20 -M 60 \
  -t 0.20 -t 0.25
```

The source table is never rewritten, even with `--force`. Plain `.txt` and gzip-compressed `.txt.gz` tables are supported. Set `--out-prefix` and `--outdir` independently of the source filename so multiple runs can share the same input. Choose distinct output prefixes or directories for combinations that would otherwise produce the same filenames, such as different downsampling seeds.

Existing threshold-specific calls are reused. When a requested threshold is absent, HEARTY calculates its calls from `count_A`, `count_C`, `count_G`, and `count_T` during summarization without adding columns to the source file. Older tables without counts must already contain the requested `Base_<threshold>` and `Status_<threshold>` columns, including `0.05` for minor-frequency output. Downsampling always requires raw counts.

The usual `--min-depth`, `--max-depth`, and downsampling options apply. `--reuse-existing` is unnecessary to select this input mode; existing downstream outputs require `--force` or `--reuse-existing` to refresh them. `--dry-run` validates the input header and previews summaries without changing files. Region and custom ANGSD settings cannot be applied to this input mode: regenerate the original counts if those filters need changing.

### Single-pass downstream processing and parallelism

Minor-frequency and all requested window summaries are generated together in **one full scan** of the basecall table. Each passing site is downsampled once (including all requested replicates), its mean allele frequencies are calculated once, and those frequencies are reused for the internal `0.05` threshold and every requested `-t` threshold. Both all-sites and transversion-only window counts are updated together. No per-threshold or per-replicate file scans are performed.

For example, reuse an existing whole-genome basecall table with 15-read downsampling, 10 replicates, and two window thresholds:

```bash
python3 HEARTY.py \
  --basecall-file /path/to/sample.basecall.txt.gz \
  -o sample_ds15 -d results/downsampled \
  -m 15 -M 80 \
  -t 0.25 -t 0.30 \
  --downsample-depth 15 --downsample-reps 10 --downsample-seed 1 \
  --downstream-workers 8 \
  --downstream-chunk-size 10000 \
  --progress-interval 30 \
  -c 2
```

Use your original depth filters and seed when comparing against previous results. Neither ANGSD nor basecall generation runs in this input mode, and the input table stays unchanged. Use a new output directory/prefix for comparisons; add `--force` only when you intend to replace existing summaries.

- `--downstream-workers` controls **processes**, not threads. The default is `1`, which still uses the optimized single-pass implementation. With multiple workers, substantial chunks of sites are processed independently; each worker returns aggregated minor-frequency and window counts, not per-site results. There is no thread or process per replicate.
- `--downstream-chunk-size` is the number of input rows in each task (default `10000`); rows failing filters are removed before submission. Larger chunks reduce scheduling overhead but use more memory and delay progress updates, especially with one worker. Smaller chunks can improve load balancing. At most twice the worker count is submitted at once, so the entire genome is never queued in memory. Memory is used for these bounded chunks and accumulated window-level counts, not all sites.
- `--progress-interval` sets the reporting interval in seconds (default `30`); `0` disables progress. Reports go to stderr and show completed input rows, passing sites, approximate input-byte percentage, elapsed time, processing rate, and estimated time remaining. For `.gz` input, the byte position refers to the compressed file. Read-ahead, variable compression, and uneven site costs make percentage and ETA approximate, especially early in a run. Counts advance as chunks are merged; single-worker reports occur between chunks. Completion is reported only after outputs have been written.
- `-c`/`--cores` continues to control ANGSD/compression threads, **not** downstream worker processes. Request enough CPUs from your scheduler for the workers plus the parent reader and compression processes; avoid assigning all allocated CPUs to both settings. Samples in a BAM list are still processed sequentially, with the chosen worker count per sample.

The coordinate-derived RNG seed and downsampling algorithm are unchanged, so worker count, chunk size, and scheduling do not alter the sampled frequencies. Output naming, frequency binning, coverage rules, window boundaries, and all-sites/transversion counting retain the previous behavior. Without downsampling, existing call columns and their stored frequencies are still respected; missing thresholds are calculated from raw counts. Summary files are staged after successful computation before replacing their targets.

Regression tests compare the new summaries against the previous serial routines, including scaffold changes, coordinate gaps, chunk boundaries, and different worker counts. Run them with `python3 -m unittest discover -s tests -v`.

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

If downsampling is enabled, the downstream outputs are tagged, for example:
- `{prefix}.mincov{depth}.downsampled_depth10_reps10.minorfreq.txt`
- `{prefix}_het0.25.mincov{depth}.downsampled_depth10_reps10.windows.txt.gz`
- `{prefix}_het0.25.mincov{depth}.downsampled_depth10_reps10.windows_trv.txt.gz`

Important behavior:
- `0.05` always runs internally for the exploratory minor allele frequency output, even if you only specify stricter thresholds
- `--min-depth` and `--max-depth` are applied downstream from the raw basecall table, not during ANGSD counting
- `--reuse-existing` lets you rebuild downstream summaries from an existing basecall file without rerunning ANGSD
- if you rerun with `--reuse-existing` and request a new threshold, HEARTY will rebuild the basecall table from the existing `.pos.gz` and `.counts.gz` files and keep the old threshold columns
- if `--downsample-depth` is used, only sites with total depth greater than or equal to that depth contribute to the downsampled summaries

### Custom ANGSD parameters

HEARTY uses a fixed default ANGSD command for generating nucleotide counts, but you can add extra ANGSD options with `--angsd-params`.

Example:

```bash
python3 HEARTY.py \
  -b sample.bam \
  -o sample01 \
  --angsd-params "-minmapq 30 -baq 1"
```

The value must be quoted as one string. HEARTY appends these options to the ANGSD command after its defaults, so the final command is visible with `--dry-run`. Explicit `-minq` and `-minmapq` values replace HEARTY's defaults rather than being passed twice.

This option is intended for ANGSD filtering and model settings such as mapping-quality, base-alignment-quality, or BAQ-related choices. HEARTY still manages the required counting/output options, including `-i`, `-out`, `-docounts`, `-dumpCounts`, and `-nThreads`.

### Reference allele annotation

Supply `--ref` to add the reference nucleotide at each retained position to the basecall table:

```bash
# Once per reference, if the index does not already exist:
samtools faidx /path/to/reference.fa

python3 HEARTY.py \
  -b sample.bam \
  --ref /path/to/reference.fa \
  -o sample01 -d results -t 0.25
```

The FASTA must be **uncompressed**, with its matching `.fai` index alongside it. It must use the same assembly and contig names as the BAM. HEARTY looks up the one-based genomic coordinates using the index without loading the whole genome into memory. Missing contigs and out-of-range positions are errors, not silently assigned `N`.

The optional `Ref` column is appended to the table. Bases are uppercase; `N` and IUPAC ambiguity codes are preserved as annotations. `Ref` is the nucleotide in the supplied assembly, **not necessarily the sample's major allele**. It does not change HEARTY's count-based calls, minor-frequency calculations, depth filters, or downsampling. The reference is also passed to ANGSD as `-ref` for options that need it, such as BAQ; supply it through `--ref`, not additionally through `--angsd-params`.

To annotate an existing basecall table without rerunning ANGSD, repeat its original BAM/output-prefix command with `--ref` and `--reuse-existing`. HEARTY adds or updates `Ref` atomically; it preserves existing calls and counts. If you also request new threshold columns, the normal rebuild from existing `.pos.gz`/`.counts.gz` applies. Subsequent threshold rebuilds preserve `Ref` even when `--ref` is omitted. Adding a reference during reuse does not retroactively change ANGSD filters such as BAQ: regenerate the counts if those settings need changing.

`--ref` is not accepted with the read-only `--basecall-file` input mode, which never edits or copies the source table. Annotate the original table first if needed; summaries can read tables with or without `Ref`. Without `--ref`, newly generated tables retain the existing format without a reference column.

### Basecall format

The basecall table is gzip-compressed and tab-delimited. It contains:

- position columns copied from `pos.gz` such as `chr`, `pos`, `totDepth`
- raw count columns:
  - `count_A`
  - `count_C`
  - `count_G`
  - `count_T`
- frequency columns: `A`, `C`, `G`, `T`
- optional `Ref` column (when `--ref` was supplied), containing the reference nucleotide
- for each threshold:
  - `Base_<threshold>`
  - `Status_<threshold>`

Example threshold-specific columns:

```text
Base_0.05   Status_0.05   Base_0.25   Status_0.25
```

This is the format expected by `PSMC_Tool.sh`.

### Downsampling mode

When `--downsample-depth` is provided, HEARTY keeps the basecall table unchanged but builds the downstream summaries from replicate-averaged downsampled counts.

Per-site logic:
- the site must first pass the normal `--min-depth` / `--max-depth` filters
- the sum of its raw A/C/G/T counts must also be at least `--downsample-depth`
- the observed nucleotide counts are downsampled without replacement to the requested depth
- this is repeated `--downsample-reps` times
- mean allele frequencies across replicates are calculated
- thresholding is applied to those mean frequencies

That means:
- the `minorfreq` table is built from the mean downsampled `Base_0.05` / `Status_0.05` calls
- the window heterozygosity summaries are built from the mean downsampled calls at each requested threshold

This mode is intended for comparing samples across different coverages without letting deeper samples dominate simply because their allele proportions are estimated more precisely.

### Downsampling parameters

These five options serve two different purposes: the first three control **which reads are sampled and how the calls are made**; the last two control **how the work is divided among processes**. You can use the worker options without enabling downsampling.

| Option | What it does | Practical choice |
| --- | --- | --- |
| `--downsample-depth X` | At each site, randomly selects `X` bases **without replacement** from its observed `count_A`, `count_C`, `count_G`, and `count_T`. Sites with fewer than `X` counted bases are excluded from minor-frequency and window summaries, including the callable-site denominator. The original basecall file is unchanged. | Pick a common depth that your samples can support. For example, `15` makes each retained site contribute a sample of 15 bases. This option enables downsampling; without it, the other two downsampling settings have no effect. |
| `--downsample-reps N` | Draws `N` separate samples of `X` bases **from the full observed counts at each site**. HEARTY averages the A/C/G/T proportions from those draws, then applies `0.05` and all requested `-t` thresholds to that one average. It does not make `N` files or treat `N × X` as extra sequencing depth. Default: `10`. | Use `10` to match the default; `1` uses a single random draw. More replicates reduce random variation from the downsampling procedure but cannot recover alleles absent from the observed reads. More replicates take more CPU time. |
| `--downsample-seed S` | Sets the deterministic random sequence for each genomic coordinate. The same counts, coordinate, depth, replicate count, and seed give the same result regardless of worker count or chunk size. Default: `1`. | Keep the seed fixed for comparisons and reruns. Change it only to assess sampling variability. The seed is **not included in output filenames**, so use another output prefix/directory to retain results from multiple seeds. |
| `--downstream-workers N` | Uses `N` worker **processes** to process chunks of sites. Each worker returns aggregated minor-frequency and window counts. It does not create a process or thread per replicate. Default: `1`; the single pass optimization still applies. | Try a value that fits your allocated CPUs, such as `4` or `8` for large genomes. This is separate from `-c`/`--cores`, which controls ANGSD and compression threads. Leave CPU capacity for reading and compression. |
| `--downstream-chunk-size N` | Sends up to `N` input rows per worker task **before depth filtering**. At most twice the worker count is queued. Default: `10000`. It changes scheduling and memory use, not the counts or random draws. | Start with `10000`. Larger chunks can reduce process overhead but use more memory and delay progress updates; smaller chunks improve balancing when site costs vary. |

For example, if a site has 12 A reads and 8 C reads, `--downsample-depth 15 --downsample-reps 10` draws 15 of those 20 reads ten times. HEARTY averages the ten observed A and C proportions and classifies the site once per threshold from that average. The reads in the basecall table are never replaced or expanded.

`--min-depth` and `--max-depth` are separate filters on the original site depth. They are applied before downsampling; a site must pass both the coverage filters and the counted-base requirement. Use the same coverage filters, downsample depth, replicate count, and seed when comparing samples or reproducing an earlier run.

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

## ROH Window Summary

`ROH_Window_Summary.R` takes HEARTY `*.windows.txt.gz` files, defines window-based ROH-like regions, calculates FROH summaries, and plots low-heterozygosity patterns.

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

The summary script then:
- colors windows by heterozygosity proportion
- flags windows at or below a chosen ROH-like threshold
- optionally “bridges” single-bin gaps, so an isolated non-ROH bin between two ROH bins can still be treated as part of a continuous ROH block
- calculates segment-based FROH summaries from the raw and bridged ROH calls

So this script is not calling ROH directly from reads or genotype likelihoods. It is interpreting the window-level heterozygosity profiles produced by HEARTY.

### Basic usage

Show help:

```bash
Rscript ROH_Window_Summary.R --help
```

### Single-sample usage

```bash
Rscript ROH_Window_Summary.R \
  --mode single \
  --input results/sample01_het0.25.mincov10.windows.txt.gz \
  --output-pdf results/plots/sample01_roh.pdf
```

Full example with explicit tuning parameters:

```bash
Rscript ROH_Window_Summary.R \
  --mode single \
  --input sample.ROH.txt \
  --output-pdf sample.ROH.pdf \
  --min-length-mb 20 \
  --max-scaffolds 40 \
  --fill-cap 0.0005 \
  --bin-size-bp 500000 \
  --fallback-resolution 100 \
  --roh-threshold 0.0001 \
  --roh-threshold-mode absolute \
  --froh-cutoffs "100kb,1mb,5mb"
```

Single-sample mode writes one PDF and automatically saves the full terminal report as a plain-text file beside it. For example, `results/plots/sample01_roh.pdf` produces `results/plots/sample01_roh_summary.txt`. The results are still printed to the terminal, including:
- number of windows parsed
- number of scaffolds plotted
- mean and median heterozygosity across all windows and across non-ROH windows
- raw and bridged ROH bin counts
- FROH values for segments longer than the requested cutoffs, by default 100 kb, 1 Mb, and 5 Mb

The report also records the input filename, dataset type (all-sites or transversion-only), effective ROH threshold, plotting settings, and genome-size denominator used for FROH. To choose a different report filename, add:

```bash
--summary-out results/summaries/sample01_roh.txt
```

In single mode this is a human-readable report containing the same statistics as the terminal, not a TSV. Batch mode continues to write its existing multi-sample TSV using `--summary-out`. A subsequent single run with the same report filename replaces the previous report, so use distinct filenames for parameter comparisons. The report path must differ from both the input file and the PDF output.

### Batch usage

```bash
Rscript ROH_Window_Summary.R \
  --mode batch \
  --sample-list roh_samples.tsv \
  --summary-out roh_summary.tsv \
  --combined-pdf roh_combined.pdf
```

Full batch example:

```bash
Rscript ROH_Window_Summary.R \
  --mode batch \
  --sample-list Sample_table.txt \
  --summary-out ROH_batch_summary.tsv \
  --combined-pdf ROH_batch_combined.pdf \
  --min-length-mb 20 \
  --max-scaffolds 40 \
  --fill-cap 0.0005 \
  --bin-size-bp 500000 \
  --fallback-resolution 100 \
  --roh-threshold 0.0001 \
  --roh-threshold-mode absolute \
  --froh-cutoffs "100kb,1mb,5mb"
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

The batch summary TSV now includes:
- `sample`
  Sample name from the sample list.
- `roh_file`
  Input HEARTY windows file used for that sample.
- `n_windows`
  Number of parsed 100 kb windows in the input file.
- `mean_heterozygosity_all_windows`
  Mean heterozygosity proportion across all parsed windows.
- `median_heterozygosity_all_windows`
  Median heterozygosity proportion across all parsed windows.
- `n_non_roh_raw_windows`
  Number of windows not classified as ROH using the raw window threshold.
- `mean_heterozygosity_non_roh_raw_windows`
  Mean heterozygosity across windows not classified as raw ROH.
- `median_heterozygosity_non_roh_raw_windows`
  Median heterozygosity across windows not classified as raw ROH.
- `n_non_roh_bridged_windows`
  Number of windows outside ROH after isolated non-ROH gaps between two ROH windows have been bridged.
- `mean_heterozygosity_non_roh_bridged_windows`
  Mean heterozygosity across windows outside the bridged ROH definition.
- `median_heterozygosity_non_roh_bridged_windows`
  Median heterozygosity across windows outside the bridged ROH definition.
- `mean_heterozygosity_plotted_windows`
  Mean heterozygosity proportion across the subset of windows that were actually plotted after scaffold-length filtering.
- `median_heterozygosity_plotted_windows`
  Median heterozygosity proportion across the subset of plotted windows.
- `genome_bp_for_froh`
  Genome-size denominator used for FROH, calculated as `n_windows * 100000`.
- `roh_threshold`
  Effective ROH cutoff used for raw ROH calling in the summary. In relative modes, this is sample-specific.
- `roh_threshold_mode`
  How the ROH cutoff was interpreted: `absolute` or `median`.
- `roh_threshold_input`
  The value supplied to `--roh-threshold`. In `median` mode this is multiplied by the sample's median window heterozygosity.
- `n_raw_segments`
  Number of raw ROH segments before bridging.
- `n_bridged_segments`
  Number of ROH segments after single-bin bridging is applied.
- `raw_froh_gt<cutoff>` and `bridged_froh_gt<cutoff>`
  FROH values from raw and bridged ROH segments longer than each requested cutoff. With the defaults, these are `raw_froh_gt100kb`, `raw_froh_gt1Mb`, `raw_froh_gt5Mb`, and the corresponding `bridged_froh_*` columns.
- `raw_roh_bp_gt<cutoff>` and `bridged_roh_bp_gt<cutoff>`
  Total ROH base pairs contributing to each requested raw and bridged FROH class.
- `scaffolds_plotted`
  Number of scaffolds retained for plotting after the minimum-length filter.
- `rows_used_for_plot`
  Number of original windows contributing to the plotted scaffolds.
- `bin_size_bp_used`
  Plotting bin size used when aggregating windows for the PDF.

Non-ROH heterozygosity is calculated from the original 100 kb windows before plotting aggregation or scaffold-length filtering. In relative threshold modes, the sample-specific ROH cutoff is first calculated from all parsed windows, and the non-ROH statistics are then calculated after removing windows classified as ROH.

### Required inputs

- `--mode`
  Must be either `single` or `batch`.
- Single mode:
  Requires `--input` for one window-summary file and `--output-pdf` for one output PDF.
  Optionally use `--summary-out` to set the text report path; otherwise `<output-pdf without .pdf>_summary.txt` is written automatically.
- Batch mode:
  Requires `--sample-list`, `--summary-out`, and `--combined-pdf`.

### Interpreting the plot

- Window fill color:
  Lower heterozygosity windows are shown toward the low end of the color scale and higher heterozygosity windows toward the high end.
- Black ticks:
  Windows whose heterozygosity proportion is less than or equal to the effective ROH threshold. In `absolute` mode this is `--roh-threshold`; in `median` mode it is the requested proportion of that sample's median window heterozygosity.
- Orange ticks:
  “Bridged” ROH calls, where a single non-ROH window between two ROH windows is filled in to create a more continuous segment.

This gives you two related views:
- a strict raw ROH-like signal
- a slightly smoothed version that is less sensitive to one noisy window interrupting a longer low-heterozygosity tract

### Parameter meanings

- `--roh-threshold`
  The value used to mark a window as ROH-like. Lower values are stricter. How this value is interpreted depends on `--roh-threshold-mode`.
- `--roh-threshold-mode`
  Accepts `absolute` (default) or `median`, selecting an absolute heterozygosity cutoff or a relative cutoff based on each sample.
  `absolute` means a window is ROH-like if its heterozygosity is less than or equal to `--roh-threshold`.
  `median` means the effective cutoff is `--roh-threshold * median window heterozygosity`.
  For example, `--roh-threshold 0.2 --roh-threshold-mode median` calls windows ROH-like if they are at or below 20% of that sample's median window heterozygosity. The median is calculated over all parsed input windows before ROH classification and independently for each sample in batch mode, not just the windows retained for plotting. If the median is zero, the effective cutoff is zero and only zero-heterozygosity windows are classified as ROH-like.
  Mean-based ROH thresholds are no longer supported; mean heterozygosity summaries and the mean-based `--fill-cap auto` visualization option are unchanged.
- `--froh-cutoffs`
  Comma-separated ROH segment-length cutoffs used for FROH reporting. Values may use `bp`, `kb`, `mb`, or `gb`; values without a unit are interpreted as base pairs. The default is `100kb,1mb,5mb`. For example, `--froh-cutoffs "100kb,500kb,2mb,10mb"` reports FROH and ROH base pairs for all four requested length classes.
- `--fill-cap`
  Caps the displayed color scale at a maximum heterozygosity value so very high windows do not dominate the gradient. This changes the visualization, not the ROH calling itself.
  You can also use `--fill-cap auto` to derive the cap from the mean heterozygosity of the input windows.
- `--fill-cap-scale`
  Multiplier applied when `--fill-cap auto` is used.
  For example, `--fill-cap auto --fill-cap-scale 0.8` uses 80% of the mean heterozygosity as the display cap.
- `--min-length-mb`
  Minimum scaffold length to include in the plot. Short scaffolds are excluded from the visual output.
- `--max-scaffolds`
  Maximum number of scaffolds to plot per sample, ranked by scaffold length.
- `--bin-size-bp`
  Bin size used for plotting aggregation. This controls visual smoothing in the PDF and does not change the original HEARTY windows file.
- `--fallback-resolution`
  PDF fallback resolution used by `ggsave`; mainly relevant for rendering quality and file export.

### FROH in the outputs

The script reports FROH for both raw and bridged ROH using these default segment cutoffs:
- `>100 kb`
- `>1 Mb`
- `>5 Mb`

Custom cutoffs can be supplied with, for example:

```bash
--froh-cutoffs "100kb,500kb,2mb,10mb"
```

The requested cutoffs are sorted, duplicate values are removed, and the summary TSV columns are generated from the resulting labels. Decimal labels are made column-safe, so a cutoff of `2.5mb` produces columns ending in `gt2p5Mb`.

The genome size denominator is:
- `number_of_windows * 100,000 bp`

In single-sample mode the script also prints:
- mean heterozygosity across all parsed windows
- median heterozygosity across all parsed windows
- mean heterozygosity across plotted scaffolds only
- median heterozygosity across plotted scaffolds only

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
- For a visualization-driven adaptive scale, try:
  `--fill-cap auto --fill-cap-scale 0.8`
- For a sample-aware ROH cutoff, try:
  `--roh-threshold 0.2 --roh-threshold-mode median`
- If the PDF is too large:
  increase `--bin-size-bp` such as `1000000`
  reduce `--max-scaffolds` such as `20`
  lower `--fallback-resolution` such as `80`

When `--fill-cap auto` is used:
- in single mode, the cap is based on that sample
- in batch mode, each sample gets its own cap for visualisation

## PSMC Helper

`PSMC_Tool.sh` prepares masked PSMC input using a HEARTY basecall table and a chosen threshold column, and then runs PSMC on the masked result.

### Basic usage

Show help:

```bash
bash PSMC_Tool.sh --help
```

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

Full example with explicit masking and bootstrap settings:

```bash
bash PSMC_Tool.sh \
  --prefix SampleA \
  --out-dir PSMC_out \
  --bam /path/to/SampleA.bam \
  --threads 8 \
  --bed /path/to/regions.bed \
  --reference /path/to/reference.fa \
  --min-coverage-diploid 10 \
  --min-depth-goodhets 10 \
  --max-depth-goodhets 80 \
  --het-file /path/to/SampleA.basecall.txt.gz \
  --threshold 0.25 \
  --bootstrap-rounds 50
```

What it does:
- builds a diploid fastq with `bcftools`
- extracts HEARTY-supported het sites from `Status_<threshold>`
- compares those against `seqtk listhet`
- masks unsupported candidate hets
- runs `psmc`
- optionally runs bootstrap replicates

So the workflow is:
1. build a diploid fastq backbone from the BAM
2. keep only HEARTY-supported heterozygous sites from the basecall table
3. mask candidate hets seen by `seqtk` that are not supported by HEARTY
4. convert the masked fastq into `psmcfa`
5. run the main PSMC analysis
6. optionally run bootstrap replicates and combine them

### Required inputs

- `--prefix`
  Sample name used as the output stem.
- `--out-dir`
  Output directory for all intermediate and final files.
- `--bam`
  Input BAM file used to build the diploid fastq backbone.
- `--threads`
  Number of threads for `bcftools` and number of parallel bootstrap jobs.
- `--bed`
  BED file of regions passed to `bcftools mpileup -R`.
- `--reference`
  Reference fasta used by `bcftools mpileup`.
- `--min-coverage-diploid`
  Minimum coverage used by `vcfutils.pl vcf2fq -d` when constructing the diploid fastq.
- `--min-depth-goodhets`
  Minimum `totDepth` retained from the HEARTY basecall file when extracting trusted heterozygous sites.
- `--het-file`
  HEARTY basecall table, either `.txt` or `.txt.gz`.
- `--threshold`
  Threshold column to use from the HEARTY basecall table, for example `0.25` for `Status_0.25`.

### Parameter meanings

- `--min-coverage-diploid`
  Coverage threshold for the diploid fastq generation step.
  This affects the underlying PSMC sequence backbone.
- `--min-depth-goodhets`
  Minimum depth required for a HEARTY HET site to be retained in `goodhets.bed`.
  This affects which heterozygous sites are preserved rather than masked.
- `--max-depth-goodhets`
  Optional upper depth cutoff for the HEARTY-supported het list.
  Useful for excluding suspicious high-depth sites.
- `--threshold`
  Selects which `Status_<threshold>` column is used from the HEARTY basecall file.
  This is the most important biological choice for how strict the heterozygous masking is.
- `--bootstrap-rounds`
  Number of bootstrap PSMC replicates to run after the main analysis.
- `--psmc-args`
  Raw argument string passed to `psmc`.
  Use this if you want to change the demographic patterning parameters.
- `--skip-diploid`
  Reuse an existing diploid fastq instead of rebuilding it.
- `--skip-mask`
  Reuse an existing masked fastq and existing `psmcfa` intermediates.
- `--skip-bootstrap`
  Skip the bootstrap replicate stage.

### Main outputs

The script writes:
- `<prefix>_diploid.fq.gz`
  Diploid fastq generated from the BAM.
- `<prefix>_goodhets.bed`
  HEARTY-supported heterozygous sites retained after the `goodhets` depth filter.
- `<prefix>_bcfhets.bed`
  Candidate heterozygous sites found from the diploid fastq with `seqtk listhet`.
- `<prefix>_Tomask.bed`
  Candidate hets to mask because they were seen in the diploid fastq but not supported by HEARTY.
- `<prefix>_diploid_masked.fq.gz`
  Masked diploid fastq used for PSMC conversion.
- `<prefix>.psmcfa`
  Main PSMC input file.
- `<prefix>.split.psmcfa`
  Split PSMC input used for bootstraps.
- `<prefix>.psmc`
  Primary PSMC result.
- `<prefix>.round-<n>.psmc`
  Bootstrap replicate outputs.
- `<prefix>.combined.psmc`
  Combined file containing the main run plus bootstrap replicates.

### Recommended defaults

For most runs, start with:
- `--threshold 0.25`
- `--min-coverage-diploid 10`
- `--min-depth-goodhets 10`
- `--bootstrap-rounds 50`

Add:
- `--max-depth-goodhets 80`

if you want to exclude very high-depth sites from the retained HET list.

### Practical guidance

- Keep `--threshold` matched to the HEARTY threshold you trust from the minor allele frequency plots.
- If you rerun PSMC repeatedly while tuning only the masking logic, `--skip-diploid` can save time.
- If you already have the masked fastq and `psmcfa` files you want, `--skip-mask` can skip straight to PSMC.
- If you only want the main run first, use `--skip-bootstrap` and add bootstraps later.

## Notes

- `0.05` is always included internally so you can inspect the exploratory minor allele frequency spectrum before choosing stricter thresholds.
- Window outputs are only written for thresholds you explicitly request with `-t`.
- The raw basecall table is intentionally left unfiltered for coverage so you can rerun downstream summaries with different `--min-depth` and `--max-depth` values.
- `Plot_MinorFreq.R`, `ROH_Window_Summary.R`, and `PSMC_Tool.sh` are standalone helpers and can be run independently of the main HEARTY workflow once the relevant HEARTY outputs exist.

## Repository Layout

- `HEARTY.py` - main HEARTY caller
- `Plot_MinorFreq.R` - minor allele frequency overlay plotter
- `ROH_Window_Summary.R` - ROH window summary and plotting utility
- `PSMC_Tool.sh` - PSMC preparation and run helper

## Contact

micwe@dtu.dk
