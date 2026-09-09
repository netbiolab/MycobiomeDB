#!/usr/bin/env python3
"""MycoProfiler test suite (standard library only).

    python -m unittest discover -s tests -v
    python tests/test_mycoprofiler.py

Everything here runs against `tests/mock_kraken2.py`, NOT a real Kraken2 or a real
database, so the suite finishes in seconds and needs no downloads. It therefore
verifies MycoProfiler's own logic -- argument handling, database resolution, the
stage-1 -> stage-2 hand-off, pair synchronisation, counting units, overwrite
protection and failure propagation -- and says nothing about classification
accuracy. See README.md ("Validation status") for what was checked against real
Kraken2 and real databases.
"""

import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, "data")
MOCK = os.path.join(HERE, "mock_kraken2.py")
sys.path.insert(0, ROOT)

from mycoprofiler import CONFIDENCE                                  # noqa: E402
from mycoprofiler.cli import main as cli_main                        # noqa: E402
from mycoprofiler.dbconfig import ConfigError, resolve_databases     # noqa: E402
from mycoprofiler.kraken import build_classify_command, parse_kraken_stats  # noqa: E402
from mycoprofiler.pipeline import (                                  # noqa: E402
    PipelineError, default_sample_id, prepare_outdir, validate_reads,
)

# The synthetic library in tests/data: 5 BACTERIA + 3 FUNGUS + 2 NOVEL pairs.
TOTAL_PAIRS = 10
BACTERIAL_PAIRS = 5
FUNGAL_PAIRS = 3
NOVEL_PAIRS = 2


def make_mock_db(directory, marker):
    """A directory that passes MycoProfiler's Kraken2-database check, for the mock tool."""
    os.makedirs(directory, exist_ok=True)
    for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        with open(os.path.join(directory, name), "wb") as handle:
            handle.write(b"\0" * 64)
    with open(os.path.join(directory, "mock_marker.txt"), "w") as handle:
        handle.write(marker)
    return directory


def make_db_root(base):
    root = os.path.join(base, "dbs")
    for site in ("gut", "oral", "skin", "vagina"):
        make_mock_db(os.path.join(root, "bacterial", site), "BACTERIA")
        make_mock_db(os.path.join(root, "fungal", site), "FUNGUS")
    return root


def run_cli(argv, env=None):
    """Invoke the CLI in-process, with the mock kraken2 wired in via --kraken2."""
    old = dict(os.environ)
    if env:
        os.environ.update(env)
    try:
        return cli_main(argv)
    finally:
        os.environ.clear()
        os.environ.update(old)


class TempCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="mycoprofiler_test_")
        self.db_root = make_db_root(self.tmp)
        self.out = os.path.join(self.tmp, "out")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def base_argv(self, site="gut", paired=True, outdir=None):
        argv = ["run", "--site", site, "--db-root", self.db_root,
                "-o", outdir or self.out, "--kraken2", MOCK]
        if paired:
            argv += ["-1", os.path.join(DATA, "sample_R1.fastq"),
                     "-2", os.path.join(DATA, "sample_R2.fastq")]
        else:
            argv += ["-U", os.path.join(DATA, "sample_single.fastq")]
        return argv

    def manifest(self, outdir=None):
        path = os.path.join(outdir or self.out, "mycoprofiler_run_manifest.json")
        with open(path) as handle:
            return json.load(handle)

    def summary(self, outdir=None):
        path = os.path.join(outdir or self.out, "mycoprofiler_summary.tsv")
        with open(path) as handle:
            header = handle.readline().rstrip("\n").split("\t")
            values = handle.readline().rstrip("\n").split("\t")
        return dict(zip(header, values))


# --------------------------------------------------------------------------- #
class TestInputValidation(unittest.TestCase):
    def test_rejects_no_reads(self):
        with self.assertRaisesRegex(PipelineError, "no input reads"):
            validate_reads()

    def test_rejects_half_a_pair(self):
        with self.assertRaisesRegex(PipelineError, "both mates"):
            validate_reads(r1=os.path.join(DATA, "sample_R1.fastq"))

    def test_rejects_paired_and_single_together(self):
        with self.assertRaisesRegex(PipelineError, "not both"):
            validate_reads(r1=os.path.join(DATA, "sample_R1.fastq"),
                           r2=os.path.join(DATA, "sample_R2.fastq"),
                           single=os.path.join(DATA, "sample_single.fastq"))

    def test_rejects_missing_file(self):
        with self.assertRaisesRegex(PipelineError, "not found"):
            validate_reads(single="/nonexistent/reads.fastq")

    def test_rejects_same_file_twice(self):
        same = os.path.join(DATA, "sample_R1.fastq")
        with self.assertRaisesRegex(PipelineError, "same file"):
            validate_reads(r1=same, r2=same)

    def test_rejects_empty_file(self):
        tmp = tempfile.mkdtemp()
        try:
            empty = os.path.join(tmp, "empty.fastq")
            open(empty, "w").close()
            with self.assertRaisesRegex(PipelineError, "empty"):
                validate_reads(single=empty)
        finally:
            shutil.rmtree(tmp)

    def test_rejects_mixed_compression(self):
        tmp = tempfile.mkdtemp()
        try:
            gz = os.path.join(tmp, "r1.fastq.gz")
            with gzip.open(gz, "wt") as handle:
                handle.write("@a\nACGT\n+\nIIII\n")
            with self.assertRaisesRegex(PipelineError, "mixed compression"):
                validate_reads(r1=gz, r2=os.path.join(DATA, "sample_R2.fastq"))
        finally:
            shutil.rmtree(tmp)

    def test_sample_id_derivation(self):
        self.assertEqual(default_sample_id("/x/SRR123_R1_001.fastq.gz"), "SRR123")
        self.assertEqual(default_sample_id("/x/SRR123_1.fq.gz"), "SRR123")
        self.assertEqual(default_sample_id("/x/SRR123.fastq"), "SRR123")
        self.assertEqual(default_sample_id("/x/SRR123_1_QC.fastq.gz"), "SRR123_1")


class TestDatabaseResolution(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_db_root_layout(self):
        resolved = resolve_databases("skin", db_root="/db", environ={})
        self.assertTrue(resolved["bacterial"].path.endswith("/db/bacterial/skin"))
        self.assertTrue(resolved["fungal"].path.endswith("/db/fungal/skin"))

    def test_explicit_beats_root(self):
        resolved = resolve_databases("gut", bacterial_db="/explicit", db_root="/db", environ={})
        self.assertEqual(resolved["bacterial"].path, "/explicit")
        self.assertTrue(resolved["fungal"].path.endswith("/db/fungal/gut"))

    def test_config_beats_root(self):
        config = os.path.join(self.tmp, "databases.json")
        with open(config, "w") as handle:
            json.dump({"_note": "ignored", "oral": {"bacterial": "/cfg/bac"}}, handle)
        resolved = resolve_databases("oral", db_config=config, db_root="/db", environ={})
        self.assertEqual(resolved["bacterial"].path, "/cfg/bac")
        self.assertTrue(resolved["fungal"].path.endswith("/db/fungal/oral"))

    def test_environment_fallback(self):
        resolved = resolve_databases("vagina", environ={"MYCOPROFILER_DB_ROOT": "/envdb"})
        self.assertTrue(resolved["bacterial"].path.endswith("/envdb/bacterial/vagina"))

    def test_no_source_is_an_actionable_error(self):
        with self.assertRaises(ConfigError) as ctx:
            resolve_databases("gut", environ={})
        self.assertIn("--db-root", str(ctx.exception))

    def test_bad_json_is_reported(self):
        config = os.path.join(self.tmp, "bad.json")
        with open(config, "w") as handle:
            handle.write("{not json")
        with self.assertRaisesRegex(ConfigError, "not valid JSON"):
            resolve_databases("gut", db_config=config, environ={})

    def test_unknown_site_rejected(self):
        with self.assertRaises(KeyError):
            resolve_databases("lung", db_root="/db", environ={})


class TestCommandConstruction(unittest.TestCase):
    def test_confidence_is_fixed_on_both_stages(self):
        for paired in (True, False):
            reads = ["a.fq", "b.fq"] if paired else ["a.fq"]
            command = build_classify_command("kraken2", "/db", reads, "o", "r", paired=paired)
            self.assertIn("--confidence", command)
            self.assertEqual(command[command.index("--confidence") + 1], str(CONFIDENCE))
            self.assertEqual(str(CONFIDENCE), "0.2")

    def test_paired_unclassified_uses_hash_placeholder(self):
        command = build_classify_command(
            "kraken2", "/db", ["a.fq", "b.fq"], "o", "r", paired=True,
            unclassified_prefix="/tmp/pfx",
        )
        self.assertIn("--paired", command)
        self.assertEqual(command[command.index("--unclassified-out") + 1], "/tmp/pfx#.fastq")

    def test_single_end_has_no_paired_flag(self):
        command = build_classify_command(
            "kraken2", "/db", ["a.fq"], "o", "r", paired=False, unclassified_prefix="/tmp/pfx"
        )
        self.assertNotIn("--paired", command)
        self.assertEqual(command[command.index("--unclassified-out") + 1], "/tmp/pfx.fastq")

    def test_gzip_flag_only_for_gzipped_input(self):
        gz = build_classify_command("kraken2", "/db", ["a.fq.gz", "b.fq.gz"], "o", "r", paired=True)
        plain = build_classify_command("kraken2", "/db", ["a.fq", "b.fq"], "o", "r", paired=True)
        self.assertIn("--gzip-compressed", gz)
        self.assertNotIn("--gzip-compressed", plain)

    def test_stats_unit_depends_on_layout(self):
        lines = ["10 sequences (0.50 Mbp) processed in 1s.",
                 "  5 sequences classified (50.00%)",
                 "  5 sequences unclassified (50.00%)"]
        self.assertEqual(parse_kraken_stats(lines, True)["unit"], "read_pairs")
        self.assertEqual(parse_kraken_stats(lines, False)["unit"], "reads")


class TestOverwriteProtection(TempCase):
    def test_refuses_non_empty_directory(self):
        os.makedirs(self.out)
        with open(os.path.join(self.out, "something.txt"), "w") as handle:
            handle.write("x")
        with self.assertRaisesRegex(PipelineError, "refusing to write"):
            prepare_outdir(self.out)

    def test_force_allows_overwrite(self):
        os.makedirs(self.out)
        with open(os.path.join(self.out, "something.txt"), "w") as handle:
            handle.write("x")
        with self.assertRaisesRegex(PipelineError, "unrelated files"):
            prepare_outdir(self.out, force=True)

    def test_second_run_without_force_exits_nonzero(self):
        self.assertEqual(run_cli(self.base_argv() + ["-q"]), 0)
        self.assertEqual(run_cli(self.base_argv() + ["-q"]), 1)
        self.assertEqual(run_cli(self.base_argv() + ["-q", "--force"]), 0)


class TestPairedEndRun(TempCase):
    def setUp(self):
        super(TestPairedEndRun, self).setUp()
        self.assertEqual(run_cli(self.base_argv() + ["-q", "--keep-intermediate"]), 0)

    def test_stage_counts_and_units(self):
        manifest = self.manifest()
        stage1 = manifest["stage1_bacterial_decontamination"]
        stage2 = manifest["stage2_fungal_classification"]
        self.assertEqual(stage1["counting_unit"], "read_pairs")
        self.assertEqual(stage1["processed"], TOTAL_PAIRS)
        self.assertEqual(stage1["classified_to_prokaryotic_catalogue"], BACTERIAL_PAIRS)
        self.assertEqual(stage1["unclassified_carried_forward"],
                         FUNGAL_PAIRS + NOVEL_PAIRS)
        # Stage 2 sees exactly what stage 1 left unclassified -- nothing else.
        self.assertEqual(stage2["processed"], FUNGAL_PAIRS + NOVEL_PAIRS)
        self.assertEqual(stage2["classified_to_fungal_catalogue"], FUNGAL_PAIRS)
        self.assertEqual(stage2["unclassified"], NOVEL_PAIRS)

    def test_handoff_carries_the_right_reads(self):
        manifest = self.manifest()
        handoff = manifest["handoff"]
        self.assertEqual(handoff["counting_unit"], "read_pairs")
        self.assertEqual(handoff["carried_forward"], FUNGAL_PAIRS + NOVEL_PAIRS)
        self.assertTrue(handoff["mate_record_counts_equal"])
        self.assertEqual(handoff["intermediate_files"], "kept")

        def read_names(path):
            with open(path) as handle:
                return [line[1:].split("/")[0] for line in handle if line.startswith("@read")]

        r1, r2 = sorted(handoff["intermediate_paths"])
        names1, names2 = read_names(r1), read_names(r2)
        # R1/R2 correspondence preserved: same names, same order, no BACTERIA reads.
        self.assertEqual(names1, names2)
        self.assertEqual(len(names1), FUNGAL_PAIRS + NOVEL_PAIRS)
        self.assertFalse(any("BACTERIA" in name for name in names1))
        self.assertEqual(sum("FUNGUS" in name for name in names1), FUNGAL_PAIRS)

    def test_summary_reports_both_units(self):
        summary = self.summary()
        self.assertEqual(summary["counting_unit"], "read_pairs")
        self.assertEqual(summary["input_units"], str(TOTAL_PAIRS))
        self.assertEqual(summary["input_reads_total"], str(TOTAL_PAIRS * 2))
        self.assertEqual(summary["stage2_classified_fungal"], str(FUNGAL_PAIRS))
        self.assertEqual(summary["stage2_classified_percent_of_input"], "30.00")

    def test_manifest_records_reproducibility_fields(self):
        manifest = self.manifest()
        self.assertEqual(manifest["mycoprofiler"]["confidence_stage1"], 0.2)
        self.assertEqual(manifest["mycoprofiler"]["confidence_stage2"], 0.2)
        self.assertIn("Kraken version", manifest["environment"]["kraken2_version"])
        self.assertIn("--confidence 0.2", manifest["stage1_bacterial_decontamination"]["command"])
        self.assertIn("--confidence 0.2", manifest["stage2_fungal_classification"]["command"])
        self.assertTrue(manifest["databases"]["bacterial"]["fingerprint"]["path"])
        self.assertTrue(manifest["databases"]["fungal"]["fingerprint"]["files"]["hash.k2d"])

    def test_expected_output_files_exist(self):
        for relative in (
            "01_bacterial_decontamination/sample.bacterial.kraken2.report",
            "01_bacterial_decontamination/sample.bacterial.kraken2.output",
            "01_bacterial_decontamination/sample.bacterial.kraken2.log",
            "02_fungal_classification/sample.fungal.kraken2.report",
            "02_fungal_classification/sample.fungal.kraken2.output",
            "02_fungal_classification/sample.fungal.kraken2.log",
            "mycoprofiler_summary.tsv",
            "mycoprofiler_run_manifest.json",
            "mycoprofiler.log",
        ):
            self.assertTrue(os.path.isfile(os.path.join(self.out, relative)), relative)


class TestSingleEndRun(TempCase):
    def test_single_end_counts_in_reads(self):
        self.assertEqual(run_cli(self.base_argv(paired=False) + ["-q"]), 0)
        manifest = self.manifest()
        stage1 = manifest["stage1_bacterial_decontamination"]
        self.assertEqual(stage1["counting_unit"], "reads")
        self.assertEqual(stage1["processed"], TOTAL_PAIRS)
        self.assertEqual(manifest["sample"]["library_layout"], "single-end")
        self.assertIsNone(manifest["handoff"]["mate_record_counts_equal"])
        self.assertEqual(self.summary()["input_reads_total"], str(TOTAL_PAIRS))


class TestGzippedInput(TempCase):
    def test_gzipped_paired_input(self):
        for mate in ("R1", "R2"):
            source = os.path.join(DATA, "sample_%s.fastq" % mate)
            target = os.path.join(self.tmp, "sample_%s.fastq.gz" % mate)
            with open(source, "rb") as src, gzip.open(target, "wb") as dst:
                shutil.copyfileobj(src, dst)
        argv = ["run", "--site", "oral", "--db-root", self.db_root, "-o", self.out,
                "--kraken2", MOCK, "-q",
                "-1", os.path.join(self.tmp, "sample_R1.fastq.gz"),
                "-2", os.path.join(self.tmp, "sample_R2.fastq.gz")]
        self.assertEqual(run_cli(argv), 0)
        manifest = self.manifest()
        self.assertIn("--gzip-compressed", manifest["stage1_bacterial_decontamination"]["command"])
        self.assertEqual(manifest["stage1_bacterial_decontamination"]["processed"], TOTAL_PAIRS)
        # Stage 2 reads the UNCOMPRESSED intermediates Kraken2 wrote, so the flag
        # must not be carried over.
        self.assertNotIn("--gzip-compressed", manifest["stage2_fungal_classification"]["command"])


class TestIntermediateHandling(TempCase):
    def test_intermediates_deleted_by_default(self):
        self.assertEqual(run_cli(self.base_argv() + ["-q"]), 0)
        manifest = self.manifest()
        self.assertEqual(manifest["handoff"]["intermediate_files"], "deleted")
        leftovers = os.listdir(os.path.join(self.out, "intermediate_unclassified_reads"))
        self.assertEqual(leftovers, [])

    def test_gzip_intermediate_implies_keep(self):
        self.assertEqual(run_cli(self.base_argv() + ["-q", "--gzip-intermediate"]), 0)
        manifest = self.manifest()
        self.assertEqual(manifest["handoff"]["intermediate_files"], "kept_gzipped")
        for path in manifest["handoff"]["intermediate_paths"]:
            self.assertTrue(path.endswith(".gz"))
            self.assertTrue(os.path.isfile(path))
            with gzip.open(path, "rt") as handle:
                self.assertTrue(handle.readline().startswith("@read"))


class TestFailurePropagation(TempCase):
    def test_stage1_failure_stops_the_run(self):
        code = run_cli(self.base_argv() + ["-q"], env={"MOCK_KRAKEN2_FAIL": "BACTERIA"})
        self.assertEqual(code, 1)
        # Stage 2 must not have run, and no summary may claim success.
        self.assertFalse(os.path.exists(os.path.join(self.out, "mycoprofiler_summary.tsv")))
        self.assertFalse(os.path.exists(os.path.join(self.out, "mycoprofiler_run_manifest.json")))
        stage2 = os.path.join(self.out, "02_fungal_classification")
        self.assertEqual([n for n in os.listdir(stage2)], [])

    def test_stage2_failure_stops_the_run(self):
        code = run_cli(self.base_argv() + ["-q"], env={"MOCK_KRAKEN2_FAIL": "FUNGUS"})
        self.assertEqual(code, 1)
        self.assertFalse(os.path.exists(os.path.join(self.out, "mycoprofiler_run_manifest.json")))
        # Stage 1 did run, so its log is there to diagnose from.
        self.assertTrue(os.path.isfile(os.path.join(
            self.out, "01_bacterial_decontamination", "sample.bacterial.kraken2.log")))

    def test_missing_database_is_reported_before_running(self):
        argv = ["run", "--site", "gut", "--db-root", os.path.join(self.tmp, "nope"),
                "-o", self.out, "--kraken2", MOCK, "-q",
                "-1", os.path.join(DATA, "sample_R1.fastq"),
                "-2", os.path.join(DATA, "sample_R2.fastq")]
        self.assertEqual(run_cli(argv), 1)
        self.assertFalse(os.path.exists(self.out))

    def test_incomplete_database_is_reported(self):
        broken = os.path.join(self.tmp, "broken")
        os.makedirs(broken)
        open(os.path.join(broken, "hash.k2d"), "w").close()
        argv = ["run", "--site", "gut", "--db-root", self.db_root,
                "--bacterial-db", broken, "-o", self.out, "--kraken2", MOCK, "-q",
                "-1", os.path.join(DATA, "sample_R1.fastq"),
                "-2", os.path.join(DATA, "sample_R2.fastq")]
        self.assertEqual(run_cli(argv), 1)


class TestDryRunAndCheck(TempCase):
    def test_dry_run_writes_nothing(self):
        self.assertEqual(run_cli(self.base_argv() + ["--dry-run"]), 0)
        self.assertFalse(os.path.exists(self.out))

    def test_check_passes_with_all_databases(self):
        self.assertEqual(run_cli(["check", "--db-root", self.db_root, "--kraken2", MOCK]), 0)

    def test_check_fails_with_missing_databases(self):
        self.assertEqual(
            run_cli(["check", "--db-root", os.path.join(self.tmp, "absent"), "--kraken2", MOCK]), 1
        )


class TestConsoleEntryPoint(TempCase):
    def test_module_invocation(self):
        argv = [sys.executable, "-m", "mycoprofiler"] + self.base_argv() + ["-q"]
        result = subprocess.run(argv, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        self.assertTrue(os.path.isfile(os.path.join(self.out, "mycoprofiler_summary.tsv")))


if __name__ == "__main__":
    unittest.main(verbosity=2)
