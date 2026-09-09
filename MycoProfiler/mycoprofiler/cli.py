"""Command-line interface: `mycoprofiler run` and `mycoprofiler check`."""

import argparse
import os
import sys

from . import CONFIDENCE, __version__
from .dbconfig import ConfigError, ENV_DB_CONFIG, ENV_DB_ROOT, resolve_databases, validate_db
from .external import ExternalToolError, find_executable, tool_version
from .pipeline import PipelineError, run
from .sites import SITE_NAMES, SITES

DESCRIPTION = """\
MycoProfiler -- body-site-specific fungal profiling of human shotgun metagenomes.

Stage 1 classifies the input reads with Kraken2 against the body site's prokaryotic
genome catalogue and sets aside the reads Kraken2 leaves UNCLASSIFIED. Stage 2
classifies exactly those reads with Kraken2 against the same body site's MycobiomeDB
fungal database. Both stages run at --confidence %s. Kraken2 is a k-mer classifier;
no alignment is performed at any point.
""" % CONFIDENCE

EPILOG = """\
database locations (first match wins):
  --bacterial-db / --fungal-db   explicit paths
  --db-config FILE               JSON: {"gut": {"bacterial": "...", "fungal": "..."}}
  $%s                 same, from the environment
  --db-root DIR                  expects DIR/{bacterial,fungal}/<site>/
  $%s                   same, from the environment

examples:
  mycoprofiler run --site gut  -1 S1_R1.fastq.gz -2 S1_R2.fastq.gz --db-root /data/mycoprofiler_db -o out/S1 -t 16
  mycoprofiler run --site skin -U S2.fastq.gz --db-root /data/mycoprofiler_db -o out/S2 -t 8
  mycoprofiler check --site vagina --db-root /data/mycoprofiler_db
""" % (ENV_DB_CONFIG, ENV_DB_ROOT)


def _add_db_arguments(parser):
    group = parser.add_argument_group("databases")
    group.add_argument("--db-root", metavar="DIR",
                       help="root holding DIR/{bacterial,fungal}/<site>/ (env: $%s)" % ENV_DB_ROOT)
    group.add_argument("--db-config", metavar="FILE",
                       help="JSON file mapping body site to database paths (env: $%s)"
                            % ENV_DB_CONFIG)
    group.add_argument("--bacterial-db", metavar="DIR",
                       help="explicit path to the stage-1 prokaryotic Kraken2 database")
    group.add_argument("--fungal-db", metavar="DIR",
                       help="explicit path to the stage-2 MycobiomeDB Kraken2 database")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="mycoprofiler", description=DESCRIPTION, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version="mycoprofiler %s" % __version__)
    sub = parser.add_subparsers(dest="command", metavar="{run,check}")

    # ---- run -------------------------------------------------------------- #
    run_parser = sub.add_parser(
        "run", help="run both Kraken2 stages on one sample",
        description=DESCRIPTION, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    run_parser.add_argument("--site", required=True, choices=SITE_NAMES,
                            help="body site; selects both databases")

    reads = run_parser.add_argument_group("input reads (paired-end OR single-end)")
    reads.add_argument("-1", "--r1", metavar="FASTQ", help="forward mate (.fastq or .fastq.gz)")
    reads.add_argument("-2", "--r2", metavar="FASTQ", help="reverse mate (.fastq or .fastq.gz)")
    reads.add_argument("-U", "--single", metavar="FASTQ", help="single-end reads")

    run_parser.add_argument("-o", "--outdir", required=True, metavar="DIR",
                            help="output directory (must not already hold a run unless --force)")
    run_parser.add_argument("--sample-id", metavar="ID",
                            help="sample name used in output filenames "
                                 "(default: derived from the first read file)")
    run_parser.add_argument("-t", "--threads", type=int, default=1, metavar="N",
                            help="threads passed to both Kraken2 stages (default: 1)")
    _add_db_arguments(run_parser)

    behaviour = run_parser.add_argument_group("behaviour")
    behaviour.add_argument("--kraken2", metavar="PATH",
                           help="kraken2 executable to use (default: first on $PATH)")
    behaviour.add_argument("--memory-mapping", action="store_true",
                           help="pass --memory-mapping to Kraken2: read the hash table from "
                                "disk instead of loading it into RAM (slower, far less memory)")
    behaviour.add_argument("--keep-intermediate", action="store_true",
                           help="keep the stage-1 unclassified FASTQs "
                                "(default: delete them after stage 2 succeeds)")
    behaviour.add_argument("--gzip-intermediate", action="store_true",
                           help="gzip the kept intermediates; implies --keep-intermediate")
    behaviour.add_argument("--force", action="store_true",
                           help="back up recognized previous MycoProfiler outputs, then start fresh")
    behaviour.add_argument("--dry-run", action="store_true",
                           help="print the two Kraken2 commands and exit without running them")
    behaviour.add_argument("-q", "--quiet", action="store_true",
                           help="do not echo progress to stderr (mycoprofiler.log is still written)")

    # ---- check ------------------------------------------------------------ #
    check_parser = sub.add_parser(
        "check", help="verify Kraken2 and the databases without classifying anything",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    check_parser.add_argument("--site", choices=SITE_NAMES, action="append", dest="sites",
                              help="body site to check; repeatable (default: all four)")
    check_parser.add_argument("--kraken2", metavar="PATH", help="kraken2 executable to check")
    _add_db_arguments(check_parser)
    return parser


def _cmd_check(args):
    ok = True
    try:
        kraken2_exe = find_executable("kraken2", args.kraken2)
        print("kraken2: %s" % kraken2_exe)
        print("         %s" % tool_version(kraken2_exe))
    except ExternalToolError as exc:
        print("kraken2: NOT FOUND\n%s" % exc)
        ok = False

    sites = args.sites or SITE_NAMES
    for site in sites:
        record = SITES[site]
        print("\n[%s]" % site)
        for kind, catalogue in (("bacterial", record["bacterial_catalogue"]),
                                ("fungal", record["fungal_catalogue"])):
            try:
                resolved = resolve_databases(
                    site,
                    bacterial_db=args.bacterial_db, fungal_db=args.fungal_db,
                    db_config=args.db_config, db_root=args.db_root,
                )[kind]
                validate_db(resolved, kind, site)
                print("  %-10s OK    %s" % (kind, resolved.path))
                print("  %-10s       %s" % ("", catalogue))
            except ConfigError as exc:
                ok = False
                first = str(exc).splitlines()
                print("  %-10s FAIL  %s" % (kind, first[0]))
                for line in first[1:]:
                    print("  %-10s       %s" % ("", line))
    print("\n%s" % ("all checks passed" if ok else "one or more checks FAILED"))
    return 0 if ok else 1


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.command:
        parser.print_help()
        return 2

    try:
        if args.command == "check":
            return _cmd_check(args)

        if args.gzip_intermediate:
            args.keep_intermediate = True
        run(
            site=args.site, r1=args.r1, r2=args.r2, single=args.single,
            outdir=args.outdir, sample_id=args.sample_id, threads=args.threads,
            bacterial_db=args.bacterial_db, fungal_db=args.fungal_db,
            db_config=args.db_config, db_root=args.db_root, kraken2=args.kraken2,
            keep_intermediate=args.keep_intermediate,
            gzip_intermediate=args.gzip_intermediate,
            memory_mapping=args.memory_mapping, force=args.force,
            quiet=args.quiet, dry_run=args.dry_run,
        )
        return 0
    except (PipelineError, ConfigError) as exc:
        sys.stderr.write("\nmycoprofiler: error: %s\n" % exc)
        return 1
    except ExternalToolError as exc:
        sys.stderr.write("\nmycoprofiler: %s\n" % exc)
        return 1
    except KeyboardInterrupt:
        sys.stderr.write("\nmycoprofiler: interrupted\n")
        return 130


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
