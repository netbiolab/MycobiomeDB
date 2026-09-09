"""The two-stage MycoProfiler workflow.

    stage 1  reads              --Kraken2-->  site prokaryotic catalogue
    stage 2  UNCLASSIFIED reads --Kraken2-->  site MycobiomeDB fungal catalogue

"Unclassified" here means exactly what Kraken2 means by it: a read (or read pair)
whose k-mers did not produce a call against the stage-1 database at the confidence
threshold. It is not an alignment result, and it is not a set of fungal reads --
see `docs/DATABASES.md` and the README's "Scope and interpretation" section.
"""

import datetime
import json
import os
import platform
import shutil
import socket
import sys

from . import CONFIDENCE, __version__
from .dbconfig import db_fingerprint, resolve_databases, validate_db
from .external import ExternalToolError, find_executable, run_command, shell_quote, tool_version
from .kraken import (
    all_gzipped,
    build_classify_command,
    count_fastq_records,
    parse_kraken_stats,
    unclassified_paths,
)
from .sites import get_site

STAGE1_DIR = "01_bacterial_decontamination"
STAGE2_DIR = "02_fungal_classification"
INTERMEDIATE_DIR = "intermediate_unclassified_reads"

SUMMARY_NAME = "mycoprofiler_summary.tsv"
MANIFEST_NAME = "mycoprofiler_run_manifest.json"
LOG_NAME = "mycoprofiler.log"


class PipelineError(Exception):
    """A MycoProfiler-level failure (bad inputs, unusable outputs, refusing to overwrite)."""


class Logger(object):
    """Tiny tee logger: everything the user sees also lands in mycoprofiler.log."""

    def __init__(self, path=None, quiet=False):
        self.quiet = quiet
        self.handle = open(path, "a") if path else None

    def __call__(self, message):
        stamped = "[%s] %s" % (datetime.datetime.now().strftime("%H:%M:%S"), message)
        if not self.quiet:
            sys.stderr.write(stamped + "\n")
            sys.stderr.flush()
        if self.handle:
            self.handle.write(stamped + "\n")
            self.handle.flush()

    def close(self):
        if self.handle:
            self.handle.close()
            self.handle = None


# --------------------------------------------------------------------------- #
# input validation
# --------------------------------------------------------------------------- #

def validate_reads(r1=None, r2=None, single=None):
    """Check the read arguments and return (paths, paired)."""
    if single and (r1 or r2):
        raise PipelineError(
            "give either paired-end reads (-1/-2) or single-end reads (-U), not both."
        )
    if not single and not (r1 or r2):
        raise PipelineError(
            "no input reads. Use -1 R1.fastq.gz -2 R2.fastq.gz for paired-end, "
            "or -U reads.fastq.gz for single-end."
        )
    if (r1 and not r2) or (r2 and not r1):
        raise PipelineError(
            "paired-end input needs both mates: -1 and -2 must be given together. "
            "For unpaired data use -U instead."
        )

    paths = [single] if single else [r1, r2]
    paired = single is None

    for path in paths:
        if not os.path.exists(path):
            raise PipelineError("input FASTQ not found: %s" % path)
        if os.path.isdir(path):
            raise PipelineError("input FASTQ is a directory, not a file: %s" % path)
        if not os.access(path, os.R_OK):
            raise PipelineError("input FASTQ is not readable: %s" % path)
        if os.path.getsize(path) == 0:
            raise PipelineError("input FASTQ is empty: %s" % path)

    if paired:
        if os.path.abspath(paths[0]) == os.path.abspath(paths[1]):
            raise PipelineError(
                "-1 and -2 point at the same file (%s); that is not a read pair." % paths[0]
            )
        if all_gzipped(paths) is None:
            raise PipelineError(
                "mixed compression: one mate is gzip-compressed and the other is not.\n"
                "  -1 %s\n  -2 %s\n"
                "Kraken2 cannot read a mixed pair -- compress or decompress both."
                % (paths[0], paths[1])
            )
    return [os.path.abspath(p) for p in paths], paired


def default_sample_id(path):
    """Derive a sample id from a read filename, stripping the usual suffixes."""
    name = os.path.basename(path)
    for suffix in (".gz", ".gzip"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    for suffix in (".fastq", ".fq", ".fasta", ".fa"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    for suffix in ("_R1_001", "_R2_001", "_R1", "_R2", "_1", "_2", ".1", ".2", "_QC"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name or "sample"


def prepare_outdir(outdir, force=False):
    """Create the output directory, refusing to clobber a previous run."""
    outdir = os.path.abspath(outdir)
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise PipelineError("output path exists and is not a directory: %s" % outdir)
        marker = os.path.join(outdir, MANIFEST_NAME)
        existing = os.listdir(outdir)
        if existing and not force:
            hint = (
                "It already holds a completed MycoProfiler run (%s)." % MANIFEST_NAME
                if os.path.exists(marker)
                else "It is not empty (%d entr%s)."
                % (len(existing), "y" if len(existing) == 1 else "ies")
            )
            raise PipelineError(
                "refusing to write into %s.\n%s\n"
                "Choose a new -o directory, or pass --force to overwrite it."
                % (outdir, hint)
            )
        if existing and force:
            if os.path.islink(outdir):
                raise PipelineError("--force refuses a symlink output directory")
            owned = {STAGE1_DIR, STAGE2_DIR, INTERMEDIATE_DIR, MANIFEST_NAME,
                     SUMMARY_NAME, LOG_NAME}
            unknown = set(existing) - owned
            if unknown:
                raise PipelineError("--force refuses unrelated files: %s" % ", ".join(sorted(unknown)))
            backup = outdir + ".backup-" + datetime.datetime.now().strftime("%Y%m%dT%H%M%S%f")
            os.rename(outdir, backup)
            sys.stderr.write("Previous output preserved at %s\n" % backup)
    os.makedirs(outdir, exist_ok=True)
    for sub in (STAGE1_DIR, STAGE2_DIR, INTERMEDIATE_DIR):
        os.makedirs(os.path.join(outdir, sub), exist_ok=True)
    return outdir


# --------------------------------------------------------------------------- #
# the run itself
# --------------------------------------------------------------------------- #

def run(site, r1=None, r2=None, single=None, outdir=None, sample_id=None, threads=1,
        bacterial_db=None, fungal_db=None, db_config=None, db_root=None,
        kraken2=None, keep_intermediate=False, gzip_intermediate=False,
        memory_mapping=False, force=False, quiet=False, dry_run=False, environ=None):
    """Execute the MycoProfiler workflow for one sample. Returns the manifest dict."""
    record = get_site(site)
    reads, paired = validate_reads(r1=r1, r2=r2, single=single)
    sample_id = sample_id or default_sample_id(reads[0])
    if sample_id in {".", ".."} or any(c in sample_id for c in "/\\\t\r\n"):
        raise PipelineError("sample id must be a filename component without slashes or control characters")
    if outdir is None:
        raise PipelineError("an output directory is required (-o/--outdir).")

    kraken2_exe = find_executable("kraken2", kraken2)
    kraken2_version = tool_version(kraken2_exe)

    databases = resolve_databases(
        site, bacterial_db=bacterial_db, fungal_db=fungal_db,
        db_config=db_config, db_root=db_root, environ=environ,
    )
    for kind in ("bacterial", "fungal"):
        validate_db(databases[kind], kind, site)

    if threads < 1:
        raise PipelineError("--threads must be >= 1 (got %s)" % threads)

    stage1_prefix = os.path.join(INTERMEDIATE_DIR, sample_id + ".unclassified")

    if dry_run:
        outdir = os.path.abspath(outdir)
        plan = _plan(outdir, sample_id, reads, paired, threads, databases,
                     kraken2_exe, memory_mapping, stage1_prefix)
        for line in plan:
            print(line)
        return {"dry_run": True, "commands": plan}

    outdir = prepare_outdir(outdir, force=force)
    log = Logger(os.path.join(outdir, LOG_NAME), quiet=quiet)
    started = datetime.datetime.now()

    try:
        log("MycoProfiler %s  |  site=%s (%s)" % (__version__, site, record["label"]))
        log("sample=%s  layout=%s  threads=%d  confidence=%s (both stages)"
            % (sample_id, "paired-end" if paired else "single-end", threads, CONFIDENCE))
        log("stage 1 database (%s): %s" % (record["bacterial_catalogue"], databases["bacterial"].path))
        log("stage 2 database (%s): %s" % (record["fungal_catalogue"], databases["fungal"].path))

        # ---- stage 1: prokaryotic decontamination ----------------------------
        s1_report = os.path.join(outdir, STAGE1_DIR, sample_id + ".bacterial.kraken2.report")
        s1_output = os.path.join(outdir, STAGE1_DIR, sample_id + ".bacterial.kraken2.output")
        s1_log = os.path.join(outdir, STAGE1_DIR, sample_id + ".bacterial.kraken2.log")
        s1_prefix = os.path.join(outdir, stage1_prefix)

        s1_command = build_classify_command(
            kraken2_exe, databases["bacterial"].path, reads, s1_output, s1_report,
            threads=threads, paired=paired, unclassified_prefix=s1_prefix,
            memory_mapping=memory_mapping,
        )
        log("stage 1/2: Kraken2 vs the %s prokaryotic catalogue" % record["label"].lower())
        s1_lines, s1_elapsed = run_command(s1_command, "bacterial_decontamination", s1_log, echo=log)
        s1_stats = parse_kraken_stats(s1_lines, paired)
        log("stage 1 done in %.1f s: %s %s in, %s classified, %s unclassified"
            % (s1_elapsed, _fmt(s1_stats["processed"]), s1_stats["unit"],
               _fmt(s1_stats["classified"]), _fmt(s1_stats["unclassified"])))

        # ---- hand-off: the unclassified reads become stage 2's input ---------
        carried = [os.path.abspath(p) for p in unclassified_paths(s1_prefix, paired)]
        missing = [p for p in carried if not os.path.isfile(p)]
        if missing:
            raise PipelineError(
                "stage 1 finished but did not produce its unclassified read file(s):\n  %s\n"
                "Kraken2 writes these via --unclassified-out; check %s."
                % ("\n  ".join(missing), s1_log)
            )
        carried_counts = [count_fastq_records(p) for p in carried]
        if paired and carried_counts[0] != carried_counts[1]:
            raise PipelineError(
                "stage 1 unclassified mates are not the same length (R1=%s, R2=%s records).\n"
                "MycoProfiler will not run stage 2 on a de-synchronised pair."
                % (carried_counts[0], carried_counts[1])
            )
        carried_unit = "read_pairs" if paired else "reads"
        carried_n = carried_counts[0]
        log("carried forward: %s %s (%s FASTQ file%s)"
            % (_fmt(carried_n), carried_unit, len(carried), "" if len(carried) == 1 else "s"))
        if carried_n == 0:
            log("WARNING: nothing was left unclassified by stage 1; stage 2 will report "
                "zero fungal reads. This usually means the input was already heavily "
                "depleted, or the wrong body site was selected.")

        # ---- stage 2: fungal classification ----------------------------------
        s2_report = os.path.join(outdir, STAGE2_DIR, sample_id + ".fungal.kraken2.report")
        s2_output = os.path.join(outdir, STAGE2_DIR, sample_id + ".fungal.kraken2.output")
        s2_log = os.path.join(outdir, STAGE2_DIR, sample_id + ".fungal.kraken2.log")

        s2_command = build_classify_command(
            kraken2_exe, databases["fungal"].path, carried, s2_output, s2_report,
            threads=threads, paired=paired, unclassified_prefix=None,
            memory_mapping=memory_mapping,
        )
        log("stage 2/2: Kraken2 vs %s" % record["fungal_catalogue"])
        s2_lines, s2_elapsed = run_command(s2_command, "fungal_classification", s2_log, echo=log)
        s2_stats = parse_kraken_stats(s2_lines, paired)
        log("stage 2 done in %.1f s: %s %s in, %s classified as fungi, %s unclassified"
            % (s2_elapsed, _fmt(s2_stats["processed"]), s2_stats["unit"],
               _fmt(s2_stats["classified"]), _fmt(s2_stats["unclassified"])))

        # ---- intermediates ---------------------------------------------------
        intermediate_state, kept_paths = _handle_intermediates(
            carried, keep_intermediate, gzip_intermediate, log
        )

        finished = datetime.datetime.now()
        manifest = _manifest(
            site=site, record=record, sample_id=sample_id, reads=reads, paired=paired,
            outdir=outdir, threads=threads, databases=databases, kraken2_exe=kraken2_exe,
            kraken2_version=kraken2_version, memory_mapping=memory_mapping,
            s1_command=s1_command, s1_stats=s1_stats, s1_elapsed=s1_elapsed,
            s2_command=s2_command, s2_stats=s2_stats, s2_elapsed=s2_elapsed,
            carried_n=carried_n, carried_unit=carried_unit, carried_counts=carried_counts,
            intermediate_state=intermediate_state, kept_paths=kept_paths,
            started=started, finished=finished,
            outputs={
                "stage1_report": s1_report, "stage1_output": s1_output, "stage1_log": s1_log,
                "stage2_report": s2_report, "stage2_output": s2_output, "stage2_log": s2_log,
            },
        )
        with open(os.path.join(outdir, MANIFEST_NAME), "w") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")
        _write_summary(os.path.join(outdir, SUMMARY_NAME), manifest)
        log("wrote %s and %s" % (SUMMARY_NAME, MANIFEST_NAME))
        log("done: %s" % outdir)
        return manifest
    finally:
        log.close()


def _fmt(value):
    return "?" if value is None else format(value, ",d")


def _plan(outdir, sample_id, reads, paired, threads, databases, kraken2_exe,
          memory_mapping, stage1_prefix):
    prefix = os.path.join(outdir, stage1_prefix)
    s1 = build_classify_command(
        kraken2_exe, databases["bacterial"].path, reads,
        os.path.join(outdir, STAGE1_DIR, sample_id + ".bacterial.kraken2.output"),
        os.path.join(outdir, STAGE1_DIR, sample_id + ".bacterial.kraken2.report"),
        threads=threads, paired=paired, unclassified_prefix=prefix,
        memory_mapping=memory_mapping,
    )
    carried = unclassified_paths(prefix, paired)
    s2 = build_classify_command(
        kraken2_exe, databases["fungal"].path, carried,
        os.path.join(outdir, STAGE2_DIR, sample_id + ".fungal.kraken2.output"),
        os.path.join(outdir, STAGE2_DIR, sample_id + ".fungal.kraken2.report"),
        threads=threads, paired=paired, unclassified_prefix=None,
        memory_mapping=memory_mapping,
    )
    return [
        "# stage 1: prokaryotic decontamination",
        shell_quote(s1),
        "# stage 2: fungal classification of the stage-1 unclassified reads",
        shell_quote(s2),
    ]


def _handle_intermediates(paths, keep, gzip_them, log):
    """Delete, keep, or keep-and-compress the stage-1 unclassified FASTQs."""
    if not keep:
        for path in paths:
            try:
                os.remove(path)
            except OSError as exc:
                log("WARNING: could not remove intermediate %s (%s)" % (path, exc))
        log("removed stage-1 unclassified FASTQs (pass --keep-intermediate to retain them)")
        return "deleted", []
    if not gzip_them:
        log("kept stage-1 unclassified FASTQs (uncompressed) in %s"
            % os.path.dirname(paths[0]))
        return "kept", [os.path.abspath(p) for p in paths]

    gzipped = []
    import gzip as _gzip
    for path in paths:
        target = path + ".gz"
        with open(path, "rb") as src, _gzip.open(target, "wb") as dst:
            shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)
        os.remove(path)
        gzipped.append(os.path.abspath(target))
    log("kept stage-1 unclassified FASTQs, gzip-compressed")
    return "kept_gzipped", gzipped


def _manifest(**kw):
    s1, s2 = kw["s1_stats"], kw["s2_stats"]
    return {
        "mycoprofiler": {
            "version": __version__,
            "confidence_stage1": CONFIDENCE,
            "confidence_stage2": CONFIDENCE,
            "invocation": shell_quote(sys.argv),
            "started": kw["started"].isoformat(timespec="seconds"),
            "finished": kw["finished"].isoformat(timespec="seconds"),
        },
        "environment": {
            "hostname": socket.gethostname(),
            "platform": platform.platform(),
            "python": sys.version.split()[0],
            "kraken2_executable": kw["kraken2_exe"],
            "kraken2_version": kw["kraken2_version"],
            "threads": kw["threads"],
            "memory_mapping": kw["memory_mapping"],
        },
        "sample": {
            "sample_id": kw["sample_id"],
            "body_site": kw["site"],
            "library_layout": "paired-end" if kw["paired"] else "single-end",
            "input_reads": kw["reads"],
            "output_directory": kw["outdir"],
        },
        "databases": {
            "bacterial": {
                "catalogue": kw["record"]["bacterial_catalogue"],
                "resolved_from": kw["databases"]["bacterial"].source,
                "fingerprint": db_fingerprint(kw["databases"]["bacterial"].path),
            },
            "fungal": {
                "catalogue": kw["record"]["fungal_catalogue"],
                "resolved_from": kw["databases"]["fungal"].source,
                "fingerprint": db_fingerprint(kw["databases"]["fungal"].path),
            },
        },
        "stage1_bacterial_decontamination": {
            "command": shell_quote(kw["s1_command"]),
            "elapsed_seconds": round(kw["s1_elapsed"], 2),
            "counting_unit": s1["unit"],
            "processed": s1["processed"],
            "classified_to_prokaryotic_catalogue": s1["classified"],
            "classified_percent": s1["classified_percent"],
            "unclassified_carried_forward": s1["unclassified"],
            "unclassified_percent": s1["unclassified_percent"],
            "megabases_processed": s1["megabases_processed"],
        },
        "handoff": {
            "counting_unit": kw["carried_unit"],
            "carried_forward": kw["carried_n"],
            "fastq_record_counts_per_file": kw["carried_counts"],
            "mate_record_counts_equal": (len(set(kw["carried_counts"])) == 1) if kw["paired"] else None,
            "intermediate_files": kw["intermediate_state"],
            "intermediate_paths": kw["kept_paths"],
        },
        "stage2_fungal_classification": {
            "command": shell_quote(kw["s2_command"]),
            "elapsed_seconds": round(kw["s2_elapsed"], 2),
            "counting_unit": s2["unit"],
            "processed": s2["processed"],
            "classified_to_fungal_catalogue": s2["classified"],
            "classified_percent": s2["classified_percent"],
            "unclassified": s2["unclassified"],
            "unclassified_percent": s2["unclassified_percent"],
            "megabases_processed": s2["megabases_processed"],
        },
        "outputs": kw["outputs"],
    }


SUMMARY_COLUMNS = [
    "sample_id",
    "body_site",
    "library_layout",
    "counting_unit",
    "input_units",
    "input_reads_total",
    "stage1_classified_prokaryotic",
    "stage1_classified_percent",
    "stage1_unclassified",
    "stage1_unclassified_percent",
    "stage2_input_units",
    "stage2_classified_fungal",
    "stage2_classified_percent_of_stage2_input",
    "stage2_classified_percent_of_input",
    "stage1_seconds",
    "stage2_seconds",
]


def _write_summary(path, manifest):
    """One-row TSV whose unit column states what every count is counted in.

    In paired-end mode Kraken2 counts one read PAIR as one sequence, so
    `counting_unit` is `read_pairs` and `input_reads_total` (= 2 x pairs) is given
    beside it so neither unit has to be inferred.
    """
    s1 = manifest["stage1_bacterial_decontamination"]
    s2 = manifest["stage2_fungal_classification"]
    sample = manifest["sample"]
    unit = s1["counting_unit"]
    processed = s1["processed"]
    per_unit = 2 if unit == "read_pairs" else 1

    def percent(numerator, denominator):
        if not numerator or not denominator:
            return "0.00" if numerator == 0 and denominator else "NA"
        return "%.2f" % (100.0 * numerator / denominator)

    row = {
        "sample_id": sample["sample_id"],
        "body_site": sample["body_site"],
        "library_layout": sample["library_layout"],
        "counting_unit": unit,
        "input_units": processed,
        "input_reads_total": None if processed is None else processed * per_unit,
        "stage1_classified_prokaryotic": s1["classified_to_prokaryotic_catalogue"],
        "stage1_classified_percent": s1["classified_percent"],
        "stage1_unclassified": s1["unclassified_carried_forward"],
        "stage1_unclassified_percent": s1["unclassified_percent"],
        "stage2_input_units": s2["processed"],
        "stage2_classified_fungal": s2["classified_to_fungal_catalogue"],
        "stage2_classified_percent_of_stage2_input": s2["classified_percent"],
        "stage2_classified_percent_of_input": percent(
            s2["classified_to_fungal_catalogue"], processed
        ),
        "stage1_seconds": s1["elapsed_seconds"],
        "stage2_seconds": s2["elapsed_seconds"],
    }
    with open(path, "w") as handle:
        handle.write("\t".join(SUMMARY_COLUMNS) + "\n")
        handle.write(
            "\t".join("NA" if row[c] is None else str(row[c]) for c in SUMMARY_COLUMNS) + "\n"
        )
