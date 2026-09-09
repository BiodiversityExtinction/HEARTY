import csv
import gzip
import importlib.util
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "HEARTY.py"
spec = importlib.util.spec_from_file_location("hearty", SCRIPT)
hearty = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hearty)


@unittest.skipUnless(shutil.which("gzip"), "gzip is required")
class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.ref = self.root / "reference.fa"
        # Wrapped CRLF sequence, lowercase bases, ambiguous bases, and a second contig.
        self.ref.write_bytes(b">chr1\r\nacgt\r\nNRta\r\n>chr2\r\nTG\r\n")
        Path(str(self.ref) + ".fai").write_text("chr1\t8\t7\t4\t6\nchr2\t2\t26\t2\t4\n")
        self.pos = self.root / "sample.pos.gz"
        self.counts = self.root / "sample.counts.gz"
        self.out = self.root / "sample.basecall.txt.gz"
        with gzip.open(self.pos, "wt") as handle:
            handle.write("chr\tpos\ttotDepth\nchr1\t1\t10\nchr1\t4\t10\nchr1\t5\t10\nchr2\t2\t10\n")
        with gzip.open(self.counts, "wt") as handle:
            handle.write("ind0_A\tind0_C\tind0_G\tind0_T\n" + "6\t4\t0\t0\n" * 4)
        self.compression = {"compress_cmd": ["gzip", "-c"], "decompress_cmd": ["gzip", "-dc"]}

    def build(self, **kwargs):
        hearty.basecall_table(self.pos, self.counts, self.out, [0.05, 0.25], self.compression, **kwargs)

    def rows(self):
        with gzip.open(self.out, "rt") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def cli(self, *args, ok=True):
        result = subprocess.run([sys.executable, str(SCRIPT), *map(str, args)],
                                capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode == 0, ok, result.stdout + result.stderr)
        return result

    def test_indexed_lookup(self):
        with hearty.IndexedReference(self.ref) as ref:
            self.assertEqual("".join(ref.base("chr1", p) for p in range(1, 9)), "ACGTNRTA")
            self.assertEqual(ref.base("chr2", 2), "G")
            for chrom, pos in [("absent", 1), ("chr1", 0), ("chr1", 9)]:
                with self.assertRaises(ValueError):
                    ref.base(chrom, pos)

    def test_build_annotation_and_rebuild(self):
        self.build()
        original = self.rows()
        self.build(reference_path=self.ref)
        annotated = self.rows()
        self.assertEqual([row["Ref"] for row in annotated], list("ATNG"))
        self.assertEqual([{k: v for k, v in row.items() if k != "Ref"} for row in annotated], original)
        self.build(preserve_reference_from=self.out)
        self.assertEqual(self.rows(), annotated)

    def test_annotation_atomic_and_dry_run(self):
        self.build()
        hearty.annotate_basecall_reference(self.out, self.ref, self.compression)
        annotated = self.rows()
        hearty.annotate_basecall_reference(self.out, self.ref, self.compression)
        self.assertEqual(self.rows(), annotated)
        before = self.out.read_bytes()
        hearty.annotate_basecall_reference(self.out, self.ref, self.compression, dry_run=True)
        self.assertEqual(self.out.read_bytes(), before)
        Path(str(self.ref) + ".fai").write_text("wrong_contig\t8\t7\t4\t6\n")
        with self.assertRaises(ValueError):
            hearty.annotate_basecall_reference(self.out, self.ref, self.compression)
        self.assertEqual(self.out.read_bytes(), before)
        with self.assertRaises(ValueError):
            self.build(reference_path=self.ref)
        self.assertEqual(self.out.read_bytes(), before)

    def test_cli_reuse_and_read_only_input(self):
        self.build()
        bam = self.root / "sample.bam"
        bam.touch()
        # ANGSD must never execute in reuse mode.
        args = ["-b", bam, "-o", "sample", "-d", self.root, "--angsd", shutil.which("false"), "--reuse-existing"]
        self.cli(*args, "--ref", self.ref)
        self.assertEqual([row["Ref"] for row in self.rows()], list("ATNG"))
        self.cli(*args, "-t", "0.2")
        self.assertTrue(all("Status_0.2" in row and "Ref" in row for row in self.rows()))
        before = self.out.read_bytes()
        dry_run = self.cli(*args, "--force", "--dry-run", "--ref", self.ref)
        self.assertIn(f"-ref {self.ref}", dry_run.stdout)
        self.assertEqual(self.out.read_bytes(), before)
        self.cli("--basecall-file", self.out, "-o", "summary", "-d", self.root, "-t", "0.3")
        self.assertEqual(self.out.read_bytes(), before)
        self.cli("--basecall-file", self.out, "-o", "summary", "--ref", self.ref, ok=False)
        self.cli(*args, "--ref", self.ref, "--angsd-params", "-ref elsewhere.fa", ok=False)
        Path(str(self.ref) + ".fai").unlink()
        result = self.cli(*args, "--ref", self.ref, ok=False)
        self.assertIn("samtools faidx", result.stderr)

    def test_mismatched_rows_leave_previous_output_intact(self):
        self.build(reference_path=self.ref)
        before = self.out.read_bytes()
        with gzip.open(self.counts, "wt") as handle:
            handle.write("A\tC\tG\tT\n6\t4\t0\t0\n")
        with self.assertRaises(ValueError):
            self.build(reference_path=self.ref)
        self.assertEqual(self.out.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
