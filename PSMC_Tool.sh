#!/usr/bin/env bash

set -euo pipefail

show_help() {
  cat <<'EOF'
PSMC_Tool.sh - Prepare masked PSMC input using HEARTY heterozygous sites

Usage:
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
    --threshold 0.25 \
    --plot --generation-time 10 --mutation-rate 2.5e-8 --sample-age 15000

Replot an existing <out-dir>/<prefix>.combined.psmc without rerunning PSMC:
  bash PSMC_Tool.sh --prefix SAMPLE --out-dir /path/to/output_dir \
    --plot-only --generation-time 10 --mutation-rate 2.5e-8 --sample-age 15000

Required arguments for a full PSMC run (plot-only needs prefix and out-dir):
  --prefix            Sample stem used for output filenames
  --out-dir           Directory for all outputs
  --bam               Input BAM file
  --threads           Number of threads / parallel bootstrap jobs
  --bed               BED file of regions passed to bcftools mpileup -R
  --reference         Reference fasta
  --min-coverage-diploid
                      Minimum coverage used by vcfutils.pl vcf2fq -d
  --min-depth-goodhets
                      Minimum totDepth retained from HEARTY basecall file
  --het-file          HEARTY basecall file (.txt or .txt.gz)

Optional arguments:
  --max-depth-goodhets
                      Maximum totDepth retained from HEARTY basecall file
  --threshold         Threshold column to use from HEARTY basecall file [default: 0.25]
  --bootstrap-rounds  Number of bootstrap replicates [default: 50]
  --psmc-args         Quoted PSMC arguments [default: -N25 -t10 -r5 -p 4+25*2+4+6]
  --skip-bootstrap    Skip bootstrap replicates
  --skip-diploid      Reuse existing diploid fastq instead of regenerating it
  --skip-mask         Reuse existing masked fastq and psmcfa intermediates
  --plot              Plot the combined PSMC result after the run
  --plot-only         Plot an existing <prefix>.combined.psmc; skip all analysis
  --plot-input        Existing .psmc input for --plot-only [default: <out-dir>/<prefix>.combined.psmc]
  --generation-time   Years per generation; required with --plot or --plot-only
  --mutation-rate     Per-base mutation rate per generation; required for plotting
  --sample-age        Sample age in years before present, added to plotted times [default: 0]
  --psmc-plot         Path to psmc_plot.pl [default: search PATH]
  --help              Show this help message

Plotting requires psmc_plot.pl, gnuplot and epstopdf. Sample age shifts every
curve (including bootstraps) to years before present; it does not modify PSMC.

Outputs:
  <out-dir>/<prefix>_diploid.fq.gz
  <out-dir>/<prefix>_goodhets.bed
  <out-dir>/<prefix>_bcfhets.bed
  <out-dir>/<prefix>_Tomask.bed
  <out-dir>/<prefix>_diploid_masked.fq.gz
  <out-dir>/<prefix>.psmcfa
  <out-dir>/<prefix>.split.psmcfa
  <out-dir>/<prefix>.psmc
  <out-dir>/<prefix>.round-<n>.psmc
  <out-dir>/<prefix>.combined.psmc
  <out-dir>/<prefix>.psmc.pdf                  (with --plot or --plot-only)
  <out-dir>/<prefix>.psmc_plot.txt             (curve and shifted time table)
EOF
}

die() {
  echo "Error: $*" >&2
  exit 1
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found in PATH: $1"
}

log_step() {
  echo
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

if [[ -f "${HOME}/.bashrc" ]]; then
  # Preserve older environments where tools are loaded via bashrc.
  # shellcheck disable=SC1090
  source "${HOME}/.bashrc"
fi

prefix=""
out_dir=""
bam=""
threads=""
bed=""
reference=""
min_coverage_diploid=""
min_depth_goodhets=""
max_depth_goodhets=""
het_file=""
threshold="0.25"
bootstrap_rounds=50
psmc_args=(-N25 -t10 -r5 -p "4+25*2+4+6")
skip_bootstrap=0
skip_diploid=0
skip_mask=0
plot=0
plot_only=0
generation_time=""
mutation_rate=""
sample_age="0"
psmc_plot=""
plot_input=""
sample_age_set=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix) prefix="${2:-}"; shift 2 ;;
    --out-dir) out_dir="${2:-}"; shift 2 ;;
    --bam) bam="${2:-}"; shift 2 ;;
    --threads) threads="${2:-}"; shift 2 ;;
    --bed) bed="${2:-}"; shift 2 ;;
    --reference) reference="${2:-}"; shift 2 ;;
    --min-coverage-diploid) min_coverage_diploid="${2:-}"; shift 2 ;;
    --min-depth-goodhets) min_depth_goodhets="${2:-}"; shift 2 ;;
    --max-depth-goodhets) max_depth_goodhets="${2:-}"; shift 2 ;;
    --het-file) het_file="${2:-}"; shift 2 ;;
    --threshold) threshold="${2:-}"; shift 2 ;;
    --bootstrap-rounds) bootstrap_rounds="${2:-}"; shift 2 ;;
    --psmc-args)
      [[ $# -ge 2 ]] || die "--psmc-args requires a value"
      # shellcheck disable=SC2206
      psmc_args=($2)
      shift 2
      ;;
    --skip-bootstrap) skip_bootstrap=1; shift ;;
    --skip-diploid) skip_diploid=1; shift ;;
    --skip-mask) skip_mask=1; shift ;;
    --plot) plot=1; shift ;;
    --plot-only) plot=1; plot_only=1; shift ;;
    --plot-input) plot_input="${2:-}"; shift 2 ;;
    --generation-time) generation_time="${2:-}"; shift 2 ;;
    --mutation-rate) mutation_rate="${2:-}"; shift 2 ;;
    --sample-age) sample_age="${2:-}"; sample_age_set=1; shift 2 ;;
    --psmc-plot) psmc_plot="${2:-}"; shift 2 ;;
    --help|-h) show_help; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$prefix" ]] || die "--prefix is required"
[[ -n "$out_dir" ]] || die "--out-dir is required"
[[ "$prefix" != */* && "$prefix" != "." && "$prefix" != ".." ]] || die "--prefix must be a filename stem"

if [[ "$plot" -eq 1 ]]; then
  number_re='^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$'
  [[ "$generation_time" =~ $number_re ]] || die "--generation-time must be a positive number of years"
  [[ "$mutation_rate" =~ $number_re ]] || die "--mutation-rate must be a positive number per base per generation"
  [[ "$sample_age" =~ $number_re ]] || die "--sample-age must be a non-negative number of years"
  awk -v n="$generation_time" 'BEGIN { exit !(n > 0) }' || die "--generation-time must be greater than zero"
  awk -v n="$mutation_rate" 'BEGIN { exit !(n > 0) }' || die "--mutation-rate must be greater than zero"
  if [[ -z "$psmc_plot" ]]; then
    psmc_plot="$(command -v psmc_plot.pl || true)"
    [[ -n "$psmc_plot" ]] || die "psmc_plot.pl not found in PATH; provide --psmc-plot PATH"
  fi
  [[ -f "$psmc_plot" ]] || die "psmc_plot.pl not found: $psmc_plot"
  psmc_plot="$(cd "$(dirname -- "$psmc_plot")" && pwd -P)/$(basename -- "$psmc_plot")"
  require_cmd perl
  require_cmd gnuplot
  require_cmd epstopdf
elif [[ "$sample_age_set" -eq 1 || -n "$generation_time" || -n "$mutation_rate" || -n "$psmc_plot" ]]; then
  die "Plot settings require --plot or --plot-only"
fi
if [[ -n "$plot_input" && "$plot_only" -eq 0 ]]; then
  die "--plot-input is only supported with --plot-only"
fi

combined_psmc="${out_dir}/${prefix}.combined.psmc"
plot_input="${plot_input:-$combined_psmc}"
plot_pdf="${out_dir}/${prefix}.psmc.pdf"
plot_table="${out_dir}/${prefix}.psmc_plot.txt"

make_plot() {
  local plot_tmp table_file curve_file curve_index
  [[ -f "$plot_input" ]] || die "PSMC plot input not found: $plot_input"
  plot_tmp="$(mktemp -d)" || die "Could not create plotting temporary directory"
  # psmc_plot.pl uses unquoted shell paths internally, so run it on safe local names.
  cp -- "$plot_input" "$plot_tmp/input.psmc"
  log_step "Plotting PSMC (generation time=$generation_time years, mutation rate=$mutation_rate, sample age=$sample_age years)"
  if ! (cd "$plot_tmp" && perl "$psmc_plot" -g "$generation_time" -u "$mutation_rate" -R plot input.psmc); then
    rm -rf -- "$plot_tmp"
    die "psmc_plot.pl failed"
  fi
  [[ -s "$plot_tmp/plot.0.txt" && -s "$plot_tmp/plot.gp" ]] || {
    rm -rf -- "$plot_tmp"
    die "psmc_plot.pl did not produce curve data and a gnuplot script"
  }
  table_file="$plot_tmp/plot_data.txt"
  printf 'Curve\tTime_years_before_present\tEffective_population_size_1e4\n' > "$table_file"
  curve_index=0
  while [[ -f "$plot_tmp/plot.${curve_index}.txt" ]]; do
    curve_file="$plot_tmp/plot.${curve_index}.txt"
    awk -v age="$sample_age" -v curve="$curve_index" -v table="$table_file" '
      BEGIN { FS = OFS = "\t" }
      NF >= 2 {
        $1 = sprintf("%.15g", $1 + age)
        printf "%d\t%s\t%s\n", curve, $1, $2 >> table
        print
      }
    ' "$curve_file" > "$curve_file.shifted"
    mv -- "$curve_file.shifted" "$curve_file"
    ((curve_index += 1))
  done
  if ! (cd "$plot_tmp" && gnuplot plot.gp && epstopdf plot.eps); then
    rm -rf -- "$plot_tmp"
    die "Could not render PSMC PDF"
  fi
  [[ -s "$plot_tmp/plot.pdf" ]] || {
    rm -rf -- "$plot_tmp"
    die "PSMC plot PDF was not generated"
  }
  cp -- "$plot_tmp/plot.pdf" "$plot_pdf"
  cp -- "$table_file" "$plot_table"
  rm -rf -- "$plot_tmp"
  echo "PSMC plot: $plot_pdf"
  echo "Shifted plot data: $plot_table"
}

if [[ "$plot_only" -eq 1 ]]; then
  mkdir -p "$out_dir"
  make_plot
  exit 0
fi

[[ -n "$bam" ]] || die "--bam is required"
[[ -n "$threads" ]] || die "--threads is required"
[[ -n "$bed" ]] || die "--bed is required"
[[ -n "$reference" ]] || die "--reference is required"
[[ -n "$min_coverage_diploid" ]] || die "--min-coverage-diploid is required"
[[ -n "$min_depth_goodhets" ]] || die "--min-depth-goodhets is required"
[[ -n "$het_file" ]] || die "--het-file is required"
[[ -n "$threshold" ]] || die "--threshold is required"

[[ "$threads" =~ ^[0-9]+$ ]] || die "--threads must be an integer"
[[ "$bootstrap_rounds" =~ ^[0-9]+$ ]] || die "--bootstrap-rounds must be an integer"
[[ "$min_coverage_diploid" =~ ^[0-9]+$ ]] || die "--min-coverage-diploid must be an integer"
[[ "$min_depth_goodhets" =~ ^[0-9]+$ ]] || die "--min-depth-goodhets must be an integer"
if [[ -n "$max_depth_goodhets" ]]; then
  [[ "$max_depth_goodhets" =~ ^[0-9]+$ ]] || die "--max-depth-goodhets must be an integer"
  (( max_depth_goodhets >= min_depth_goodhets )) || die "--max-depth-goodhets must be >= --min-depth-goodhets"
fi

[[ -f "$bam" ]] || die "BAM file not found: $bam"
[[ -f "$bed" ]] || die "BED file not found: $bed"
[[ -f "$reference" ]] || die "Reference file not found: $reference"
[[ -f "$het_file" ]] || die "Het file not found: $het_file"

mkdir -p "$out_dir"

require_cmd bcftools
require_cmd vcfutils.pl
require_cmd seqtk
require_cmd bedtools
require_cmd fq2psmcfa
require_cmd splitfa
require_cmd psmc
require_cmd xargs
require_cmd awk
require_cmd grep
require_cmd gzip

if [[ "$het_file" == *.gz ]]; then
  require_cmd zcat
  het_cat=(zcat -- "$het_file")
else
  het_cat=(cat -- "$het_file")
fi

status_col="Status_${threshold}"

diploid_fq="${out_dir}/${prefix}_diploid.fq.gz"
goodhets_bed="${out_dir}/${prefix}_goodhets.bed"
bcfhets_bed="${out_dir}/${prefix}_bcfhets.bed"
to_mask_bed="${out_dir}/${prefix}_Tomask.bed"
masked_fq="${out_dir}/${prefix}_diploid_masked.fq.gz"
psmcfa="${out_dir}/${prefix}.psmcfa"
split_psmcfa="${out_dir}/${prefix}.split.psmcfa"
psmc_main="${out_dir}/${prefix}.psmc"

log_step "Preparing diploid fastq"
if [[ "$skip_diploid" -eq 1 ]]; then
  [[ -f "$diploid_fq" ]] || die "--skip-diploid was set but diploid fastq is missing: $diploid_fq"
else
  bcftools mpileup -q 20 -Q 20 --threads "$threads" -Ou -R "$bed" -f "$reference" "$bam" \
    | bcftools call --threads "$threads" -c - \
    | vcfutils.pl vcf2fq -d "$min_coverage_diploid" \
    | gzip > "$diploid_fq"
fi

log_step "Building het masks from bcftools and HEARTY outputs"
if [[ "$skip_mask" -eq 1 ]]; then
  [[ -f "$masked_fq" ]] || die "--skip-mask was set but masked fastq is missing: $masked_fq"
  [[ -f "$psmcfa" ]] || die "--skip-mask was set but psmcfa is missing: $psmcfa"
  [[ -f "$split_psmcfa" ]] || die "--skip-mask was set but split psmcfa is missing: $split_psmcfa"
else
  "${het_cat[@]}" \
    | awk -v status_col="$status_col" -v min_depth="$min_depth_goodhets" -v max_depth="${max_depth_goodhets:-}" '
        BEGIN { FS = OFS = "\t" }
        NR == 1 {
          for (i = 1; i <= NF; i++) {
            if ($i == status_col) {
              status_idx = i
            }
            if ($i == "totDepth") {
              depth_idx = i
            }
          }
          if (status_idx == 0) {
            printf("Error: Could not find column %s in basecall file\n", status_col) > "/dev/stderr"
            exit 1
          }
          if (depth_idx == 0) {
            printf("Error: Could not find totDepth column in basecall file\n") > "/dev/stderr"
            exit 1
          }
          next
        }
        {
          depth = $depth_idx + 0
          if (depth < min_depth) {
            next
          }
          if (max_depth != "" && depth > max_depth) {
            next
          }
        }
        $status_idx ~ /^HET/ { print $1, $2, $2 }
      ' > "$goodhets_bed"

  seqtk listhet "$diploid_fq" \
    | awk '{print $1"\t"$2"\t"$2}' > "$bcfhets_bed"

  bedtools intersect -loj -a "$bcfhets_bed" -b "$goodhets_bed" \
    | awk '$4=="."{print $1"\t"$2}' > "$to_mask_bed"

  seqtk seq -n N -M "$to_mask_bed" "$diploid_fq" | gzip > "$masked_fq"

  fq2psmcfa -q20 "$masked_fq" > "$psmcfa"
  splitfa "$psmcfa" > "$split_psmcfa"
fi

log_step "Running primary PSMC"
psmc "${psmc_args[@]}" -o "$psmc_main" "$psmcfa"

if [[ "$skip_bootstrap" -eq 0 && "$bootstrap_rounds" -gt 0 ]]; then
  log_step "Running ${bootstrap_rounds} PSMC bootstrap replicates"
  seq 1 "$bootstrap_rounds" \
    | xargs -P "$threads" -I{} bash -lc '
        set -euo pipefail
        out="$1"
        split="$2"
        shift 2
        psmc "$@" -b -o "$out" "$split"
      ' _ "${out_dir}/${prefix}.round-{}.psmc" "$split_psmcfa" "${psmc_args[@]}"
else
  log_step "Skipping bootstrap replicates"
fi

log_step "Combining PSMC outputs"
if compgen -G "${out_dir}/${prefix}.round-*.psmc" > /dev/null; then
  cat "$psmc_main" "${out_dir}/${prefix}.round-"*.psmc > "$combined_psmc"
else
  cp "$psmc_main" "$combined_psmc"
fi

if [[ "$plot" -eq 1 ]]; then
  make_plot
fi

log_step "Done"
echo "PSMC threshold column: $status_col"
echo "Diploid coverage cutoff: $min_coverage_diploid"
echo "Goodhets depth filter: min=$min_depth_goodhets${max_depth_goodhets:+, max=$max_depth_goodhets}"
echo "Generated files:"
echo "  $diploid_fq"
echo "  $goodhets_bed"
echo "  $bcfhets_bed"
echo "  $to_mask_bed"
echo "  $masked_fq"
echo "  $psmcfa"
echo "  $split_psmcfa"
echo "  $psmc_main"
echo "  $combined_psmc"
if [[ "$plot" -eq 1 ]]; then
  echo "  $plot_pdf"
  echo "  $plot_table"
fi
