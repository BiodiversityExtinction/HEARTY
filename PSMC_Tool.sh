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
    --threshold 0.25

Required arguments:
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
  --help              Show this help message

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
    --help|-h) show_help; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$prefix" ]] || die "--prefix is required"
[[ -n "$out_dir" ]] || die "--out-dir is required"
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
combined_psmc="${out_dir}/${prefix}.combined.psmc"

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

log_step "Done"
echo "Main output: $psmc_main"
echo "Combined output: $combined_psmc"
