import argparse
from contextlib import redirect_stdout, redirect_stderr
import gzip
import io
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import HEARTY as hearty


class DownstreamTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.compression = hearty.compression_tools(1)
        self.thresholds = [0.05, 0.25, 0.3]
        self.source = self.root / "sample.basecall.txt.gz"
        rng = random.Random(19)
        counts = [[12, 8, 0, 0], [18, 0, 2, 0], [0, 7, 0, 13], [20, 0, 0, 0],
                  [6, 7, 7, 0], [0, 0, 9, 11], [4, 4, 0, 0], [50, 30, 0, 0],
                  [5, 15, 0, 0], [6, 14, 0, 0], [0, 1, 19, 0], [7, 0, 13, 0]]
        header = ["chr", "pos", "totDepth", *[f"count_{b}" for b in "ACGT"], *"ACGT",
                  "Base_0.05", "Status_0.05", "Base_0.25", "Status_0.25", "Ref"]
        self.rows = []
        for chrom, positions in [("chr1", [*range(1, 75), *range(99995, 100040), *range(700000, 700040)]),
                                 ("chr2", range(1, 25)), ("chr3", range(900000, 900080))]:
            for pos in positions:
                c = counts[rng.randrange(len(counts))]
                freq = hearty.counts_to_freqs(c)
                calls = [v for t in [0.05, 0.25] for v in hearty.classify_counts_for_threshold(freq, t)]
                self.rows.append([chrom, str(pos), str(sum(c)), *map(str, c),
                                  *[f"{v:.2f}" if v > 0 else "0" for v in freq], *calls, "A"])
        self.text = "\t".join(header) + "\n" + "".join("\t".join(row) + "\n" for row in self.rows)
        with gzip.open(self.source, "wt") as handle:
            handle.write(self.text)
        self.original = self.source.read_bytes()

    def outputs(self, name, depth):
        result = hearty.build_outputs(name, self.root / name, 10, 60, self.thresholds[1:], depth, 10)
        result["outdir"].mkdir()
        result["basecall"] = self.source
        return result

    def args(self, depth, workers=1, chunk=13, progress=0):
        return argparse.Namespace(min_depth=10, max_depth=60, downsample_depth=depth,
                                  downsample_reps=10, downsample_seed=17, dry_run=False,
                                  downstream_workers=workers, downstream_chunk_size=chunk,
                                  progress_interval=progress)

    def legacy(self, outputs, depth):
        common = dict(min_depth=10, max_depth=60, compression=self.compression,
                      downsample_depth=depth, downsample_reps=10, downsample_seed=17)
        hearty.minorfreq_from_basecall(self.source, 0.05, outputs["minorfreq"], **common)
        for kind in ("windows", "windows_trv"):
            for threshold, path in outputs[kind].items():
                hearty.build_window_summary(self.source, threshold, path,
                                            transversions_only=kind == "windows_trv", **common)

    def snapshot(self, outputs):
        content = {"minorfreq": outputs["minorfreq"]["table"].read_text()}
        for kind in ("windows", "windows_trv"):
            for threshold, path in outputs[kind].items():
                with gzip.open(path, "rt") as handle:
                    content[(kind, threshold)] = handle.read()
        return content

    def test_exact_legacy_equivalence_across_workers_chunks_and_downsampling(self):
        with patch.object(hearty, "WINDOW_MIN_SITES", 3), redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            for depth in (None, 15):
                oracle = self.outputs(f"old{depth}", depth)
                self.legacy(oracle, depth)
                expected = self.snapshot(oracle)
                for workers, chunk in ((1, 1), (1, 37), (2, 7), (3, 55)):
                    with self.subTest(depth=depth, workers=workers, chunk=chunk):
                        outputs = self.outputs(f"new{depth}_{workers}_{chunk}", depth)
                        hearty.run_downstream("test", outputs, self.args(depth, workers, chunk), self.compression)
                        self.assertEqual(expected, self.snapshot(outputs))
        self.assertEqual(self.source.read_bytes(), self.original)

    def test_one_downsampling_per_passing_site_and_one_input_scan(self):
        outputs = self.outputs("once", 15)
        expected = sum(10 <= int(r[2]) <= 60 and sum(map(int, r[3:7])) >= 15 for r in self.rows)
        with patch.object(hearty, "downsample_mean_freqs", wraps=hearty.downsample_mean_freqs) as sampler, \
                patch.object(hearty, "downstream_input", wraps=hearty.downstream_input) as reader, \
                patch.object(hearty, "iter_compressed_lines", side_effect=AssertionError("extra file pass")), \
                redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            hearty.run_downstream("once", outputs, self.args(15), self.compression)
        self.assertEqual(sampler.call_count, expected)
        self.assertEqual(reader.call_count, 1)

    def test_frequency_only_legacy_table(self):
        legacy = "\n".join("\t".join(field for i, field in enumerate(line.split("\t"))
                                    if i not in range(3, 7)) for line in self.text.splitlines()) + "\n"
        with gzip.open(self.source, "wt") as handle:
            handle.write(legacy)
        self.thresholds = [0.05, 0.25]
        expected = self.outputs("legacy_freqs", None)
        actual = self.outputs("new_freqs", None)
        with patch.object(hearty, "WINDOW_MIN_SITES", 3), redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            self.legacy(expected, None)
            hearty.run_downstream("new_freqs", actual, self.args(None, 2), self.compression)
        self.assertEqual(self.snapshot(expected), self.snapshot(actual))

    def test_cli_basecall_plain_gzip_progress_and_validation(self):
        plain = self.root / "source.txt"
        plain.write_text(self.text)
        for source in (plain, self.source):
            outdir = self.root / ("out_" + source.name)
            command = [sys.executable, hearty.__file__, "--basecall-file", str(source),
                       "-o", "sample", "-d", str(outdir), "-m", "10", "-M", "60",
                       "-t", "0.25", "-t", "0.30", "--downsample-depth", "15",
                       "--downstream-workers", "2", "--downstream-chunk-size", "19",
                       "--progress-interval", "0.01", "--angsd", "/does/not/exist"]
            result = subprocess.run(command, capture_output=True, text=True, timeout=40)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            for message in ("rows processed", "elapsed", "ETA", "100.0%", "(complete)"):
                self.assertIn(message, result.stderr)
            self.assertEqual(len(list(outdir.iterdir())), 5)
            before = {p.name: p.read_bytes() for p in outdir.iterdir()}
            result = subprocess.run(command + ["--force", "--dry-run"], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(before, {p.name: p.read_bytes() for p in outdir.iterdir()})
            for option, value in [("--downstream-workers", "0"), ("--downstream-chunk-size", "0"),
                                  ("--progress-interval", "nan"), ("--progress-interval", "-1")]:
                result = subprocess.run(command + [option, value, "--dry-run"], capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.source.read_bytes(), self.original)

    def test_failed_worker_does_not_replace_outputs(self):
        outputs = self.outputs("failure", None)
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            hearty.run_downstream("failure", outputs, self.args(None), self.compression)
            before = {p: p.read_bytes() for p in outputs["outdir"].iterdir()}
            damaged = [line.split("\t") for line in self.text.splitlines()]
            for row in damaged[1:]:
                row[7] = "not_a_frequency"
            with gzip.open(self.source, "wt") as handle:
                handle.write("\n".join("\t".join(row) for row in damaged) + "\n")
            with self.assertRaises(ValueError):
                hearty.run_downstream("failure", outputs, self.args(None, 2), self.compression)
        self.assertEqual(before, {p: p.read_bytes() for p in outputs["outdir"].iterdir()})

    def test_truncated_gzip_does_not_replace_outputs(self):
        outputs = self.outputs("truncated", None)
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            hearty.run_downstream("truncated", outputs, self.args(None), self.compression)
            before = {p: p.read_bytes() for p in outputs["outdir"].iterdir()}
            self.source.write_bytes(self.source.read_bytes()[:-8])
            with self.assertRaises(subprocess.CalledProcessError):
                hearty.run_downstream("truncated", outputs, self.args(None, 2), self.compression)
        self.assertEqual(before, {p: p.read_bytes() for p in outputs["outdir"].iterdir()})

    def test_empty_or_all_filtered_input_and_minorfreq_only(self):
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            for empty in (False, True):
                if empty:
                    with gzip.open(self.source, "wt") as handle:
                        handle.write(self.text.splitlines()[0] + "\n")
                args = self.args(15, 2)
                args.min_depth = 1000
                outputs = self.outputs(f"empty{empty}", 15)
                hearty.run_downstream("empty", outputs, args, self.compression)
                snapshot = self.snapshot(outputs)
                self.assertEqual(snapshot["minorfreq"].count("\n"), 1)
                self.assertTrue(all(value == "" for key, value in snapshot.items() if key != "minorfreq"))
            self.source = self.root / "plain.txt"
            self.source.write_text(self.text)
            outputs = self.outputs("minoronly", 15)
            outputs["windows"] = outputs["windows_trv"] = {}
            hearty.run_downstream("minoronly", outputs, self.args(15, 2), self.compression)
            self.assertEqual(len(list(outputs["outdir"].iterdir())), 1)


if __name__ == "__main__":
    unittest.main()
