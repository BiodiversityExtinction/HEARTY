#!/usr/bin/env bash

set -euo pipefail

for cmd in perl gnuplot epstopdf psmc_plot.pl; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Skipping PSMC plot test: $cmd is unavailable"
    exit 0
  fi
done

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT

cat > "$tmp/input.psmc" <<'EOF'
MM	pattern:4+2, n:6
MM	n_iterations:25, skip:1
MM	is_decoding:0
RD	20
RI	0.1
TR	0.02	0.01
RS	0	0.0	1.0	100	0.1	0.0
RS	1	0.1	0.8	100	0.1	0.0
RS	2	0.2	1.2	100	0.1	0.0
PA	fixture
//
MM	is_decoding:1
RD	20
RI	0.1
TR	0.03	0.01
RS	0	0.0	1.0	100	0.1	0.0
RS	1	0.1	0.8	100	0.1	0.0
RS	2	0.2	1.2	100	0.1	0.0
PA	fixture
//
EOF

for age in 0 15000; do
  bash "$repo_dir/PSMC_Tool.sh" \
    --prefix "age$age" --out-dir "$tmp/out" \
    --plot-only --plot-input "$tmp/input.psmc" \
    --generation-time 10 --mutation-rate 2.5e-8 --sample-age "$age"
  test -s "$tmp/out/age$age.psmc.pdf"
  test -s "$tmp/out/age$age.psmc_plot.txt"
done

awk -F '\t' '
  NR == FNR {
    if (FNR > 1) {
      count[$1]++
      key = $1 SUBSEP count[$1]
      time[key] = $2
      popsize[key] = $3
    }
    next
  }
  FNR > 1 {
    seen[$1]++
    key = $1 SUBSEP seen[$1]
    difference = $2 - time[key] - 15000
    if (!(key in time) || difference > 0.001 || difference < -0.001 || $3 != popsize[key]) {
      bad++
    }
  }
  END {
    if (length(seen) != 2 || bad) {
      print "PSMC age offset test failed" > "/dev/stderr"
      exit 1
    }
  }
' "$tmp/out/age0.psmc_plot.txt" "$tmp/out/age15000.psmc_plot.txt"

echo "PSMC plot age-offset test passed"
