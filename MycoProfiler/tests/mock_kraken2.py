#!/usr/bin/env python3
"""A stand-in for `kraken2`, used by the MycoProfiler test suite.

It is NOT a classifier. It accepts the same arguments MycoProfiler passes to Kraken2 and
reproduces the behaviour MycoProfiler depends on:

  * `--paired` treats one read PAIR as one "sequence" in the stderr summary;
  * `--unclassified-out PREFIX#.fastq` expands `#` to `_1`/`_2` and writes both
    mates of every unclassified pair, in input order;
  * the stderr summary uses Kraken2's exact wording, which MycoProfiler parses;
  * `--report` / `--output` files are written.

Classification is deterministic: a read is "classified" if its header contains the
token stored in `<db>/mock_marker.txt`. That lets a test assert exact counts at each
stage and prove that the reads carried between stages are the right ones.

`MOCK_KRAKEN2_FAIL=<stage-token>` makes it exit 3 when the database marker matches
that token, so failure propagation can be tested.
"""

import argparse
import gzip
import os
import sys


def open_fastq(path):
    if str(path).endswith((".gz", ".gzip")):
        return gzip.open(path, "rt")
    return open(path)


def read_records(path):
    with open_fastq(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                return
            yield header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def main(argv=None):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--version", action="store_true")
    parser.add_argument("--db")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--confidence", type=float, default=0.0)
    parser.add_argument("--output")
    parser.add_argument("--report")
    parser.add_argument("--unclassified-out", dest="unclassified_out")
    parser.add_argument("--classified-out", dest="classified_out")
    parser.add_argument("--paired", action="store_true")
    parser.add_argument("--use-names", action="store_true")
    parser.add_argument("--report-zero-counts", action="store_true")
    parser.add_argument("--gzip-compressed", action="store_true")
    parser.add_argument("--memory-mapping", action="store_true")
    parser.add_argument("reads", nargs="*")
    args = parser.parse_args(argv)

    if args.version:
        print("Kraken version 2.1.6-mock")
        print("Copyright 2013-2023, Derrick Wood (dwood@cs.jhu.edu)")
        return 0

    if not args.db or not os.path.isdir(args.db):
        sys.stderr.write("mock kraken2: bad --db %r\n" % args.db)
        return 1
    marker_path = os.path.join(args.db, "mock_marker.txt")
    marker = open(marker_path).read().strip() if os.path.isfile(marker_path) else "NEVER"

    if os.environ.get("MOCK_KRAKEN2_FAIL") == marker:
        sys.stderr.write("mock kraken2: simulated failure for marker %s\n" % marker)
        return 3

    expected = 2 if args.paired else 1
    if len(args.reads) != expected:
        sys.stderr.write(
            "mock kraken2: expected %d read file(s), got %d\n" % (expected, len(args.reads))
        )
        return 1

    unclassified_handles = []
    if args.unclassified_out:
        if args.paired:
            if "#" not in args.unclassified_out:
                sys.stderr.write("mock kraken2: paired --unclassified-out needs a '#'\n")
                return 1
            targets = [args.unclassified_out.replace("#", "_1"),
                       args.unclassified_out.replace("#", "_2")]
        else:
            targets = [args.unclassified_out]
        for target in targets:
            os.makedirs(os.path.dirname(os.path.abspath(target)), exist_ok=True)
            unclassified_handles.append(open(target, "w"))

    streams = [read_records(path) for path in args.reads]
    output = open(args.output, "w") if args.output else None

    processed = classified = unclassified = bases = 0
    while True:
        try:
            records = [next(stream) for stream in streams]
        except StopIteration:
            break
        processed += 1
        bases += sum(len(record[1]) for record in records)
        headers = " ".join(record[0] for record in records)
        name = records[0][0][1:].split()[0]
        if marker in headers:
            classified += 1
            if output:
                output.write("C\t%s\t%s (taxid 2)\t%d\t\n"
                             % (name, marker, sum(len(r[1]) for r in records)))
        else:
            unclassified += 1
            if output:
                output.write("U\t%s\tunclassified (taxid 0)\t%d\t\n"
                             % (name, sum(len(r[1]) for r in records)))
            for handle, record in zip(unclassified_handles, records):
                handle.write("\n".join(record) + "\n")

    if output:
        output.close()
    for handle in unclassified_handles:
        handle.close()

    if args.report:
        total = processed or 1
        with open(args.report, "w") as report:
            report.write("%6.2f\t%d\t%d\tU\t0\tunclassified\n"
                         % (100.0 * unclassified / total, unclassified, unclassified))
            report.write("%6.2f\t%d\t0\tR\t1\troot\n"
                         % (100.0 * classified / total, classified))
            report.write("%6.2f\t%d\t%d\tS\t2\t  %s\n"
                         % (100.0 * classified / total, classified, classified, marker))

    percent = (lambda n: 100.0 * n / processed if processed else 0.0)
    sys.stderr.write("Loading database information... done.\n")
    sys.stderr.write("%d sequences (%.2f Mbp) processed in 0.001s (1.0 Kseq/m, 1.0 Mbp/m).\n"
                     % (processed, bases / 1e6))
    sys.stderr.write("  %d sequences classified (%.2f%%)\n" % (classified, percent(classified)))
    sys.stderr.write("  %d sequences unclassified (%.2f%%)\n"
                     % (unclassified, percent(unclassified)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
