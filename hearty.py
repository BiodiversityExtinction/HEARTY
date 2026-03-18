#!/usr/bin/env python3

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


DEFAULT_ANGSD = "/home/ctools/angsd-0.935/angsd"
WINDOW_SIZE_BP = 100000
WINDOW_MIN_SITES = 25000
DEFAULT_THRESHOLD = 0.05


def show_help():
    help_text = """
HEARTY - Heterozygosity Error and Ancient-damage Robust Thresholding sYstem

USAGE:
    hearty [OPTIONS]

REQUIRED ARGUMENTS:
    -b, --bam FILE             BAM file for analysis
    -l, --bam-list FILE        Two-column file: BAM path and output prefix
    -o, --out-prefix STR       Output prefix for all files (required for single BAM)

OPTIONAL ARGUMENTS:
    -d, --outdir DIR           Output directory [default: results]
    -m, --min-depth INT        Minimum depth for downstream summaries [default: 10]
    -M, --max-depth INT        Maximum depth for downstream summaries
    -t, --threshold FLOAT      Extra het calling threshold for window outputs (can be repeated; 0.05 always runs internally)
    -T, --threshold-list FILE  File with thresholds, one per line
    -R, --roh-min FLOAT        Threshold used in output filenames [default: 0.2]
    -c, --cores INT            Number of ANGSD/compression threads [default: 8]
    --angsd PATH               ANGSD executable [default: auto]
REGION SPECIFICATION (choose one):
    -r, --regions STR          Region string for ANGSD -r
    -f, --rf FILE              Regions file for ANGSD -rf

RUN OPTIONS:
    -n, --dry-run              Show commands without executing them
    -F, --force                Overwrite existing outputs
    --reuse-existing          Reuse existing ANGSD/basecall outputs when available
    -h, --help                 Show this help message

EXAMPLES:
    hearty -b sample.bam -o sample01
    hearty -l bam_list.txt
    hearty -b sample.bam -o sample01 -t 0.10 -t 0.15
    hearty -b sample.bam -o sample01 -r "chr1:1-50000000"
"""
    print(help_text)


def read_file_list(filepath, list_type="values"):
    path = Path(filepath)
    if not path.exists():
        print(f"Error: {list_type} file {filepath} not found", file=sys.stderr)
        sys.exit(1)

    items = []
    with path.open() as handle:
        for line_num, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if list_type == "bams":
                parts = line.split()
                if len(parts) < 2:
                    print(
                        f"Error: BAM list line {line_num} in {filepath} must contain both BAM path and output prefix",
                        file=sys.stderr,
                    )
                    sys.exit(1)
                bam_path = parts[0]
                prefix = parts[1]
                items.append((bam_path, prefix))
            else:
                try:
                    items.append(float(line))
                except ValueError:
                    print(
                        f"Error: Invalid {list_type} value '{line}' on line {line_num} in {filepath}",
                        file=sys.stderr,
                    )
                    sys.exit(1)

    if not items:
        print(f"Error: No valid entries found in {filepath}", file=sys.stderr)
        sys.exit(1)

    return items


def quote_cmd(parts):
    return " ".join(shlex.quote(str(part)) for part in parts)


def run_command(parts, dry_run=False, stdout=None, stderr=None):
    print(f"Running: {quote_cmd(parts)}")
    if dry_run:
        return 0
    return subprocess.run(parts, check=True, stdout=stdout, stderr=stderr).returncode


def ensure_tool(tool_name, configured_path):
    if os.path.isabs(configured_path) or os.sep in configured_path:
        if os.path.exists(configured_path):
            return configured_path
        print(f"Error: {tool_name} not found at {configured_path}", file=sys.stderr)
        sys.exit(1)

    resolved = shutil.which(configured_path)
    if resolved:
        return resolved

    print(f"Error: {tool_name} executable not found: {configured_path}", file=sys.stderr)
    sys.exit(1)


def compression_tools(threads):
    pigz = shutil.which("pigz")
    if pigz:
        thread_count = str(max(1, threads))
        return {
            "compress_cmd": [pigz, "-c", "-p", thread_count],
            "decompress_cmd": [pigz, "-d", "-c", "-p", thread_count],
            "compress_name": "pigz",
            "decompress_name": "pigz",
        }

    gzip = ensure_tool("gzip", "gzip")
    return {
        "compress_cmd": [gzip, "-c"],
        "decompress_cmd": [gzip, "-d", "-c"],
        "compress_name": "gzip",
        "decompress_name": "gzip",
    }


def coerce_rf_reg(rf, reg):
    rf = (rf or "").strip().strip('"').strip("'")
    reg = (reg or "").strip().strip('"').strip("'")
    if rf and any(token in rf for token in (":", "-", ",", ";", "|")) and not reg:
        reg = rf
        rf = ""
    return rf, reg


def maybe_remove(path, force):
    if path.exists() and force:
        path.unlink()


def check_overwrite(paths, force):
    existing = [str(path) for path in paths if path.exists()]
    if existing and not force:
        print("Error: Output files already exist. Re-run with --force to overwrite:", file=sys.stderr)
        for path in existing:
            print(f"  {path}", file=sys.stderr)
        sys.exit(1)


def coverage_tag(min_depth, max_depth):
    tag = f"mincov{min_depth}"
    if max_depth is not None:
        tag += f"_maxcov{max_depth}"
    return tag


def build_outputs(prefix, outdir, min_depth, max_depth, output_thresholds, roh_min):
    output_dir = Path(outdir)
    plots_dir = output_dir / "plots"
    cov_tag = coverage_tag(min_depth, max_depth)
    pos_gz = output_dir / f"{prefix}.pos.gz"
    counts_gz = output_dir / f"{prefix}.counts.gz"

    basecall = output_dir / f"{prefix}.basecall.txt.gz"
    windows = {
        threshold: output_dir / f"{prefix}_het{threshold}.{cov_tag}.windows.txt.gz"
        for threshold in output_thresholds
    }
    windows_trv = {
        threshold: output_dir / f"{prefix}_het{threshold}.{cov_tag}.windows_trv.txt.gz"
        for threshold in output_thresholds
    }
    minorfreq = {
        "table": output_dir / f"{prefix}.{coverage_tag(min_depth, max_depth)}.minorfreq.txt",
    }
    return {
        "outdir": output_dir,
        "plots_dir": plots_dir,
        "pos_gz": pos_gz,
        "counts_gz": counts_gz,
        "basecall": basecall,
        "windows": windows,
        "windows_trv": windows_trv,
        "minorfreq": minorfreq,
    }


def iter_compressed_lines(path, compression):
    proc = subprocess.Popen(
        compression["decompress_cmd"] + [str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        for line in proc.stdout:
            yield line
        stderr = proc.stderr.read()
        returncode = proc.wait()
        if returncode != 0:
            raise subprocess.CalledProcessError(returncode, compression["decompress_cmd"] + [str(path)], stderr=stderr)
    finally:
        if proc.stdout:
            proc.stdout.close()
        if proc.stderr:
            proc.stderr.close()


def write_compressed_lines(output_path, compression, write_fn):
    with tempfile.NamedTemporaryFile("wb", delete=False, dir=output_path.parent) as handle:
        tmp_path = Path(handle.name)

    try:
        with tmp_path.open("wb") as raw_out:
            proc = subprocess.Popen(
                compression["compress_cmd"],
                stdin=subprocess.PIPE,
                stdout=raw_out,
                stderr=subprocess.PIPE,
                text=True,
            )
            try:
                write_fn(proc.stdin)
                proc.stdin.close()
                stderr = proc.stderr.read()
                returncode = proc.wait()
                if returncode != 0:
                    raise subprocess.CalledProcessError(returncode, compression["compress_cmd"], stderr=stderr)
            finally:
                if proc.stderr:
                    proc.stderr.close()
        tmp_path.replace(output_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def threshold_label(threshold):
    return str(threshold)


def threshold_base_col(threshold):
    return f"Base_{threshold_label(threshold)}"


def threshold_status_col(threshold):
    return f"Status_{threshold_label(threshold)}"


def parse_existing_thresholds(header_names):
    thresholds = []
    for name in header_names:
        if not name.startswith("Status_"):
            continue
        raw = name[len("Status_") :]
        try:
            value = float(raw)
        except ValueError:
            continue
        if value not in thresholds:
            thresholds.append(value)
    return thresholds


def merge_thresholds(existing_thresholds, requested_thresholds):
    merged = []
    for threshold in existing_thresholds + requested_thresholds:
        if threshold not in merged:
            merged.append(threshold)
    return merged


def classify_counts_for_threshold(freqs, threshold):
    homo = 1.0 - threshold
    base = "UNKNOWN"
    status = "UNKNOWN"
    if freqs[0] > homo:
        base, status = "AA", "HOM"
    elif freqs[1] > homo:
        base, status = "CC", "HOM"
    elif freqs[2] > homo:
        base, status = "GG", "HOM"
    elif freqs[3] > homo:
        base, status = "TT", "HOM"
    elif freqs[0] > threshold and freqs[1] > threshold:
        base, status = "AC", "HET_trv"
    elif freqs[0] > threshold and freqs[2] > threshold:
        base, status = "AG", "HET_trn"
    elif freqs[0] > threshold and freqs[3] > threshold:
        base, status = "AT", "HET_trv"
    elif freqs[1] > threshold and freqs[2] > threshold:
        base, status = "CG", "HET_trv"
    elif freqs[1] > threshold and freqs[3] > threshold:
        base, status = "CT", "HET_trn"
    elif freqs[2] > threshold and freqs[3] > threshold:
        base, status = "GT", "HET_trv"
    return base, status


def parse_count_freqs(count_line):
    fields = count_line.strip().split()
    counts = [int(value) for value in fields[:4]]
    positive_sum = sum(value for value in counts if value > 0)
    freqs = []
    for value in counts:
        if positive_sum > 0 and value > 0:
            freqs.append(value / positive_sum)
        else:
            freqs.append(0.0)
    return freqs


def basecall_table(pos_gz, counts_gz, output_path, thresholds, compression, dry_run=False):
    preview = (
        f"{quote_cmd(compression['decompress_cmd'])} {shlex.quote(str(pos_gz))} + "
        f"{quote_cmd(compression['decompress_cmd'])} {shlex.quote(str(counts_gz))} "
        f"-> classify thresholds {','.join(map(str, thresholds))} -> {quote_cmd(compression['compress_cmd'])} "
        f"> {shlex.quote(str(output_path))}"
    )
    print(f"Running: {preview}")
    if dry_run:
        return

    def _writer(out_handle):
        pos_iter = iter_compressed_lines(pos_gz, compression)
        counts_iter = iter_compressed_lines(counts_gz, compression)

        first_pos = next(pos_iter, None)
        if first_pos is None:
            raise RuntimeError(f"No positions found in {pos_gz}")
        next(counts_iter, None)

        header_cols = ["A", "C", "G", "T"]
        for threshold in thresholds:
            header_cols.extend([threshold_base_col(threshold), threshold_status_col(threshold)])
        out_handle.write(first_pos.rstrip("\n") + "\t" + "\t".join(header_cols) + "\n")
        for pos_line, count_line in zip(pos_iter, counts_iter):
            freqs = parse_count_freqs(count_line)
            rendered = [f"{value:.2f}" if value > 0 else "0" for value in freqs]
            for threshold in thresholds:
                base, status = classify_counts_for_threshold(freqs, threshold)
                rendered.extend([base, status])
            out_handle.write(pos_line.rstrip("\n"))
            out_handle.write("\t")
            out_handle.write("\t".join(rendered))
            out_handle.write("\n")

    write_compressed_lines(output_path, compression, _writer)


def get_basecall_header_indices(basecall_gz, compression):
    header = next(iter_compressed_lines(basecall_gz, compression)).rstrip("\n").split("\t")
    return {name: idx for idx, name in enumerate(header)}


def depth_passes(fields, depth_idx, min_depth, max_depth):
    if depth_idx is None:
        return True
    try:
        depth = int(fields[depth_idx])
    except (ValueError, IndexError):
        return False
    if depth < min_depth:
        return False
    if max_depth is not None and depth > max_depth:
        return False
    return True


def minorfreq_from_basecall(basecall_gz, threshold, outputs, min_depth, max_depth, compression, dry_run=False):
    if dry_run:
        print(
            f"Running: {quote_cmd(compression['decompress_cmd'])} {shlex.quote(str(basecall_gz))} "
            f"-> minorfreq tables for threshold {threshold} with depth filter min={min_depth}"
            f"{'' if max_depth is None else f',max={max_depth}'}"
        )
        return

    header_index = get_basecall_header_indices(basecall_gz, compression)
    depth_idx = header_index.get("totDepth")
    base_idx = header_index[threshold_base_col(threshold)]
    status_idx = header_index[threshold_status_col(threshold)]
    a_idx = header_index["A"]
    c_idx = header_index["C"]
    g_idx = header_index["G"]
    t_idx = header_index["T"]
    pair_names = ("AC", "AG", "AT", "CG", "CT", "GT")
    counts = {}
    line_iter = iter_compressed_lines(basecall_gz, compression)
    next(line_iter)
    for line in line_iter:
        fields = line.rstrip("\n").split("\t")
        if not depth_passes(fields, depth_idx, min_depth, max_depth):
            continue
        status = fields[status_idx]
        if not status.startswith("HET"):
            continue
        trimmed = [fields[a_idx], fields[c_idx], fields[g_idx], fields[t_idx]]
        trimmed.sort()
        minor = trimmed[-3]
        pair = fields[base_idx]
        counts[(minor, pair)] = counts.get((minor, pair), 0) + 1

    by_minor = {}
    for (minor, pair), count in counts.items():
        by_minor.setdefault(minor, {name: 0 for name in pair_names})
        by_minor[minor][pair] = count

    with outputs["table"].open("w") as handle:
        handle.write("freq\tAC\tAG\tAT\tCG\tCT\tGT\ttotal\ttransition\ttransversion\n")
        for minor in sorted(by_minor, key=lambda x: float(x)):
            row = by_minor[minor]
            total = sum(row.values())
            transition = row["AG"] + row["CT"]
            transversion = total - transition
            handle.write(
                "\t".join(
                    [
                        minor,
                        str(row["AC"]),
                        str(row["AG"]),
                        str(row["AT"]),
                        str(row["CG"]),
                        str(row["CT"]),
                        str(row["GT"]),
                        str(total),
                        str(transition),
                        str(transversion),
                    ]
                )
                + "\n"
            )


def parse_basecall_site(line):
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 2:
        return None
    try:
        site = int(parts[1])
    except ValueError:
        return None
    return parts[0], site


def build_window_summary(basecall_gz, threshold, output_path, min_depth, max_depth, compression, transversions_only=False, dry_run=False):
    print(
        f"Running: {quote_cmd(compression['decompress_cmd'])} {shlex.quote(str(basecall_gz))} "
        f"-> 100kb windows for threshold {threshold} with depth filter min={min_depth}"
        f"{'' if max_depth is None else f',max={max_depth}'}"
        f" -> {quote_cmd(compression['compress_cmd'])} > {shlex.quote(str(output_path))}"
    )
    if dry_run:
        return

    def _writer(out_handle):
        header_index = get_basecall_header_indices(basecall_gz, compression)
        depth_idx = header_index.get("totDepth")
        status_idx = header_index[threshold_status_col(threshold)]
        site_count = 0
        het_count = 0
        current_window_start = 0
        current_window_end = WINDOW_SIZE_BP - 1
        current_scaffold = None

        def flush_window(force=False):
            if current_scaffold is None:
                return
            if force or site_count >= WINDOW_MIN_SITES:
                proportion = het_count / site_count if site_count > 0 else 0.0
                out_handle.write(
                    f"{current_scaffold} {current_window_start}-{current_window_end}: "
                    f"{site_count} sites, {het_count} sites with 'HET', Proportion: {proportion:.7f}\n"
                )

        nonlocal_vars = {
            "site_count": site_count,
            "het_count": het_count,
            "current_window_start": current_window_start,
            "current_window_end": current_window_end,
            "current_scaffold": current_scaffold,
        }

        line_iter = iter_compressed_lines(basecall_gz, compression)
        next(line_iter)
        for line in line_iter:
            parsed = parse_basecall_site(line)
            if parsed is None:
                continue
            scaffold, site = parsed
            fields = line.rstrip("\n").split("\t")
            if not depth_passes(fields, depth_idx, min_depth, max_depth):
                continue
            status = fields[status_idx]

            site_count = nonlocal_vars["site_count"]
            het_count = nonlocal_vars["het_count"]
            current_window_start = nonlocal_vars["current_window_start"]
            current_window_end = nonlocal_vars["current_window_end"]
            current_scaffold = nonlocal_vars["current_scaffold"]

            if scaffold != current_scaffold:
                site_count = 0
                het_count = 0
                current_window_start = 0
                current_window_end = WINDOW_SIZE_BP - 1
                current_scaffold = scaffold

            if site > current_window_end:
                if site_count >= WINDOW_MIN_SITES:
                    proportion = het_count / site_count if site_count > 0 else 0.0
                    out_handle.write(
                        f"{current_scaffold} {current_window_start}-{current_window_end}: "
                        f"{site_count} sites, {het_count} sites with 'HET', Proportion: {proportion:.7f}\n"
                    )
                site_count = 0
                het_count = 0
                current_window_start += WINDOW_SIZE_BP
                current_window_end += WINDOW_SIZE_BP

            site_count += 1
            if transversions_only:
                if status == "HET_trv":
                    het_count += 1
            elif status.startswith("HET"):
                het_count += 1

            nonlocal_vars["site_count"] = site_count
            nonlocal_vars["het_count"] = het_count
            nonlocal_vars["current_window_start"] = current_window_start
            nonlocal_vars["current_window_end"] = current_window_end
            nonlocal_vars["current_scaffold"] = current_scaffold

        flush_window(force=True)

    write_compressed_lines(output_path, compression, _writer)


def run_single_bam(bam_path, prefix, args):
    print(f"\n{'=' * 60}")
    print(f"Processing: {bam_path} -> {prefix}")
    print(f"{'=' * 60}")

    bam = Path(bam_path)
    if not bam.exists():
        print(f"Error: BAM file not found: {bam_path}", file=sys.stderr)
        return False

    angsd = ensure_tool("ANGSD", args.angsd)
    compression = compression_tools(args.cores)
    rf, reg = coerce_rf_reg(args.rf, args.regions)

    requested_thresholds = args.threshold if args.threshold else [0.05]
    output_thresholds = [threshold for threshold in requested_thresholds if threshold != DEFAULT_THRESHOLD]
    outputs = build_outputs(prefix, args.outdir, args.min_depth, args.max_depth, output_thresholds, args.roh_min)

    downstream_outputs = [
        *outputs["windows"].values(),
        *outputs["windows_trv"].values(),
        *outputs["minorfreq"].values(),
    ]
    output_paths = [
        outputs["pos_gz"],
        outputs["counts_gz"],
        outputs["basecall"],
        *downstream_outputs,
    ]
    if not args.reuse_existing:
        check_overwrite(output_paths, args.force)

    if not args.dry_run:
        outputs["outdir"].mkdir(parents=True, exist_ok=True)

    angsd_cmd = [
        angsd,
        "-minq",
        "20",
        "-minmapq",
        "20",
        "-uniqueOnly",
        "1",
        "-remove_bads",
        "1",
        "-docounts",
        "1",
        "-dumpCounts",
        "4",
        "-i",
        str(bam),
        "-out",
        str(outputs["outdir"] / prefix),
        "-nThreads",
        str(max(1, min(args.cores, 64))),
    ]
    if rf:
        angsd_cmd.extend(["-rf", rf])
    elif reg:
        angsd_cmd.extend(["-r", reg])

    log_path = outputs["outdir"] / f"{prefix}.angsd.log"
    print(f"ANGSD log: {log_path}")
    print(
        f"Compression: {compression['compress_name']} for writing, "
        f"{compression['decompress_name']} for reading"
    )

    active_thresholds = list(requested_thresholds)
    if args.reuse_existing and outputs["basecall"].exists() and not args.force:
        header_index = get_basecall_header_indices(outputs["basecall"], compression)
        existing_thresholds = parse_existing_thresholds(list(header_index.keys()))
        missing_thresholds = [threshold for threshold in requested_thresholds if threshold not in existing_thresholds]
        if missing_thresholds:
            active_thresholds = merge_thresholds(existing_thresholds, requested_thresholds)
            print(
                f"Extending existing basecall table with new thresholds: "
                f"{','.join(map(str, missing_thresholds))}"
            )
            if outputs["pos_gz"].exists() and outputs["counts_gz"].exists():
                if not args.dry_run:
                    outputs["basecall"].unlink(missing_ok=True)
                basecall_table(
                    outputs["pos_gz"],
                    outputs["counts_gz"],
                    outputs["basecall"],
                    active_thresholds,
                    compression,
                    dry_run=args.dry_run,
                )
            else:
                print(
                    "Error: Existing basecall file is missing requested thresholds and "
                    "the matching .pos.gz/.counts.gz files are not available to rebuild it.",
                    file=sys.stderr,
                )
                return False
        else:
            active_thresholds = existing_thresholds
            print(f"Reusing existing basecall table: {outputs['basecall']}")
    else:
        if args.reuse_existing and outputs["pos_gz"].exists() and outputs["counts_gz"].exists() and not args.force:
            print(f"Reusing existing ANGSD counts: {outputs['pos_gz']} and {outputs['counts_gz']}")
        else:
            if args.dry_run:
                run_command(angsd_cmd, dry_run=True)
            else:
                with log_path.open("w") as log_handle:
                    run_command(angsd_cmd, stdout=log_handle, stderr=subprocess.STDOUT)
        maybe_remove(outputs["basecall"], args.force)
        basecall_table(
            outputs["pos_gz"],
            outputs["counts_gz"],
            outputs["basecall"],
            requested_thresholds,
            compression,
            dry_run=args.dry_run,
        )
        active_thresholds = list(requested_thresholds)

    first_threshold = DEFAULT_THRESHOLD
    for minorfreq_output in outputs["minorfreq"].values():
        minorfreq_output.unlink(missing_ok=True)
    minorfreq_from_basecall(
        outputs["basecall"],
        first_threshold,
        outputs["minorfreq"],
        args.min_depth,
        args.max_depth,
        compression,
        dry_run=args.dry_run,
    )

    for threshold, output_path in outputs["windows"].items():
        output_path.unlink(missing_ok=True)
        build_window_summary(
            outputs["basecall"],
            threshold,
            output_path,
            args.min_depth,
            args.max_depth,
            compression,
            dry_run=args.dry_run,
        )

    for threshold, output_path in outputs["windows_trv"].items():
        output_path.unlink(missing_ok=True)
        build_window_summary(
            outputs["basecall"],
            threshold,
            output_path,
            args.min_depth,
            args.max_depth,
            compression,
            transversions_only=True,
            dry_run=args.dry_run,
        )

    print(f"Successfully completed: {prefix}")
    return True


def main():
    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1] in ["-h", "--help"]):
        show_help()
        sys.exit(0)

    parser = argparse.ArgumentParser(
        prog="hearty",
        description="HEARTY - Heterozygosity Error and Ancient-damage Robust Thresholding sYstem",
        add_help=False,
    )

    required = parser.add_argument_group("Required arguments (choose one)")
    bam_group = required.add_mutually_exclusive_group(required=True)
    bam_group.add_argument("-b", "--bam", metavar="FILE", help="BAM file for analysis")
    bam_group.add_argument("-l", "--bam-list", metavar="FILE", help="Two-column file: BAM path and output prefix")
    required.add_argument(
        "-o",
        "--out-prefix",
        metavar="STR",
        help="Output prefix (required for single BAM, ignored for BAM list)",
    )

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument("-d", "--outdir", default="results", metavar="DIR", help="Output directory")
    optional.add_argument("-m", "--min-depth", type=int, default=10, metavar="INT", help="Minimum depth for downstream summaries")
    optional.add_argument("-M", "--max-depth", type=int, metavar="INT", help="Maximum depth for downstream summaries")
    optional.add_argument("-c", "--cores", type=int, default=8, metavar="INT", help="ANGSD/compression threads")
    optional.add_argument("--angsd", default=DEFAULT_ANGSD, metavar="PATH", help="ANGSD executable")

    thresh_group = optional.add_mutually_exclusive_group()
    thresh_group.add_argument(
        "-t",
        "--threshold",
        type=float,
        action="append",
        default=[],
        metavar="FLOAT",
        help="Extra het calling threshold for window outputs (can be repeated; 0.05 always runs internally)",
    )
    thresh_group.add_argument("-T", "--threshold-list", metavar="FILE", help="File with thresholds")

    optional.add_argument(
        "-R",
        "--roh-min",
        type=float,
        default=0.2,
        metavar="FLOAT",
        help="Threshold used in output filenames",
    )

    regions = parser.add_argument_group("Region specification (choose one)")
    region_group = regions.add_mutually_exclusive_group()
    region_group.add_argument("-r", "--regions", metavar="STR", help="Region string for ANGSD -r")
    region_group.add_argument("-f", "--rf", metavar="FILE", help="Regions file for ANGSD -rf")

    run_group = parser.add_argument_group("Run options")
    run_group.add_argument("-n", "--dry-run", action="store_true", help="Show commands without executing")
    run_group.add_argument("-F", "--force", action="store_true", help="Overwrite existing outputs")
    run_group.add_argument("--reuse-existing", action="store_true", help="Reuse existing ANGSD/basecall outputs when available")
    run_group.add_argument("-h", "--help", action="store_true", help="Show this help message")

    args = parser.parse_args()
    if args.help:
        show_help()
        sys.exit(0)

    if args.bam and not args.out_prefix:
        print("Error: --out-prefix is required when using --bam", file=sys.stderr)
        sys.exit(1)

    extra_thresholds = []
    if args.threshold_list:
        extra_thresholds = read_file_list(args.threshold_list, "thresholds")
    elif args.threshold:
        extra_thresholds = args.threshold

    thresholds = [DEFAULT_THRESHOLD]
    for threshold in extra_thresholds:
        if threshold not in thresholds:
            thresholds.append(threshold)
    args.threshold = thresholds

    if args.bam_list:
        bam_list = read_file_list(args.bam_list, "bams")
        print(f"Found {len(bam_list)} BAM files to process")
        success_count = 0
        for bam_path, prefix in bam_list:
            if run_single_bam(bam_path, prefix, args):
                success_count += 1
        print(f"\n{'=' * 60}")
        print(f"Completed {success_count}/{len(bam_list)} BAM files successfully")
        if success_count < len(bam_list):
            sys.exit(1)
    else:
        if not run_single_bam(args.bam, args.out_prefix, args):
            sys.exit(1)


if __name__ == "__main__":
    main()
