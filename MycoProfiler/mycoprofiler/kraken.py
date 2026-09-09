"""Building and interpreting the two Kraken2 invocations MycoProfiler performs.

Kraken2 is an *exact k-mer matching* classifier: it breaks each read into k-mers,
looks each one up in a compacted hash table built from the reference genomes, and
assigns the read to the lowest common ancestor of the taxa its k-mers hit. It does
not align reads. Nothing in MycoProfiler uses an aligner.
"""

import os
import re

from . import CONFIDENCE

# "  1234 sequences classified (12.34%)"  /  "... unclassified (87.66%)"
_PROCESSED_RE = re.compile(r"^\s*([\d,]+)\s+sequences\s+\(([\d.]+)\s*Mbp\)\s+processed", re.I)
_CLASSIFIED_RE = re.compile(r"^\s*([\d,]+)\s+sequences\s+classified\s+\(([\d.]+)%\)", re.I)
_UNCLASSIFIED_RE = re.compile(r"^\s*([\d,]+)\s+sequences\s+unclassified\s+\(([\d.]+)%\)", re.I)


def _int(text):
    return int(text.replace(",", ""))


def all_gzipped(paths):
    """True if every input path is gzip-named, False if none is.

    Kraken2 cannot mix compressed and uncompressed inputs in a single paired run,
    so a mixture is rejected by the caller rather than silently mishandled.
    """
    flags = [str(p).endswith((".gz", ".gzip")) for p in paths]
    if all(flags):
        return True
    if not any(flags):
        return False
    return None


def build_classify_command(kraken2, db, reads, output, report, threads=1, paired=False,
                           unclassified_prefix=None, memory_mapping=False):
    """Assemble the argv for one Kraken2 classification.

    `unclassified_prefix` is the *stem* for unclassified reads. In paired mode a
    literal ``#`` is appended, which Kraken2 expands to ``_1`` and ``_2``; the two
    files it writes stay in the same order and hold both mates of every pair it
    emits, so R1/R2 correspondence is preserved for the next step.
    """
    command = [
        kraken2,
        "--db", str(db),
        "--threads", str(threads),
        "--confidence", str(CONFIDENCE),
        "--use-names",
        "--report-zero-counts",
        "--output", str(output),
        "--report", str(report),
    ]
    if memory_mapping:
        command.append("--memory-mapping")

    gz = all_gzipped(reads)
    if gz is True:
        command.append("--gzip-compressed")

    if unclassified_prefix is not None:
        if paired:
            # Kraken2 expands '#' to _1 / _2 for paired input.
            command += ["--unclassified-out", "%s#.fastq" % unclassified_prefix]
        else:
            command += ["--unclassified-out", "%s.fastq" % unclassified_prefix]

    if paired:
        command.append("--paired")
    command += [str(r) for r in reads]
    return command


def unclassified_paths(prefix, paired):
    """The file(s) Kraken2 will write for --unclassified-out with `prefix`."""
    if paired:
        return ["%s_1.fastq" % prefix, "%s_2.fastq" % prefix]
    return ["%s.fastq" % prefix]


def parse_kraken_stats(lines, paired):
    """Extract the processed/classified/unclassified counts from Kraken2 output.

    Kraken2 reports *sequences*. In ``--paired`` mode one "sequence" is one read
    PAIR, so the numbers are returned under explicitly named keys and the unit is
    carried alongside them; nothing downstream has to guess.
    """
    stats = {
        "unit": "read_pairs" if paired else "reads",
        "processed": None,
        "classified": None,
        "unclassified": None,
        "classified_percent": None,
        "unclassified_percent": None,
        "megabases_processed": None,
    }
    for line in lines:
        match = _PROCESSED_RE.match(line)
        if match:
            stats["processed"] = _int(match.group(1))
            stats["megabases_processed"] = float(match.group(2))
            continue
        match = _CLASSIFIED_RE.match(line)
        if match:
            stats["classified"] = _int(match.group(1))
            stats["classified_percent"] = float(match.group(2))
            continue
        match = _UNCLASSIFIED_RE.match(line)
        if match:
            stats["unclassified"] = _int(match.group(1))
            stats["unclassified_percent"] = float(match.group(2))
    return stats


def count_fastq_records(path):
    """Count records in a plain (uncompressed) FASTQ file.

    Kraken2 writes --unclassified-out uncompressed, so this never needs to handle
    gzip. Returns None if the file is absent.
    """
    if not os.path.isfile(path):
        return None
    lines = 0
    last_byte = b""
    with open(path, "rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            lines += block.count(b"\n")
            last_byte = block[-1:]
    # A final line with no trailing newline still counts.
    if last_byte and last_byte != b"\n":
        lines += 1
    return lines // 4
