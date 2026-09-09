#!/usr/bin/env python3
"""Build the MycoProfiler *skin* stage-1 Kraken2 database from the SMGC catalogue.

Unlike HRGM2, HROM and VMGC, the Skin Microbial Genome Collection is distributed as
genome FASTA only, so the Kraken2 database has to be built locally. This script does
that end to end:

  1. read the SMGC prokaryotic supplementary table (622 genomes with GTDB lineages);
  2. turn those lineages into a Kraken2-compatible taxonomy (names.dmp / nodes.dmp)
     rooted at taxid 1, with one leaf per species;
  3. rewrite every genome's FASTA headers to Kraken2's `|kraken:taxid|<taxid>` form
     and concatenate them into a single library file;
  4. run `kraken2-build --add-to-library` then `kraken2-build --build`;
  5. verify the result with `kraken2-inspect` and refuse to declare success if the
     hash table came out empty.

Only the 622 PROKARYOTIC genomes are used. The catalogue's 12 eukaryotic and 6,935
viral genomes are excluded on purpose: stage 1 exists to remove prokaryotic reads,
and putting fungi into it would remove exactly the reads stage 2 is looking for.

Inputs (see scripts/download_bacterial_db.sh, which fetches both):
  --catalogue  the extracted 03_finalcatalogue/ directory of SMGC_*.fa files
  --metadata   SMGC.xlsx (sheet "Supplementary Table 4") or a CSV export of it

Requires: kraken2-build on $PATH. pandas + openpyxl only if --metadata is .xlsx.
"""

import argparse
import csv
import datetime
import json
import os
import shutil
import subprocess
import sys

RANK_COLUMNS = [
    ("gtdb_domain", "superkingdom"),
    ("gtdb_phylum", "phylum"),
    ("gtdb_class", "class"),
    ("gtdb_order", "order"),
    ("gtdb_family", "family"),
    ("gtdb_genus", "genus"),
    ("gtdb_species", "species"),
]
XLSX_SHEET = "Supplementary Table 4"


def die(message):
    sys.stderr.write("build_smgc_kraken2_db: error: %s\n" % message)
    raise SystemExit(1)


def read_metadata(path):
    """Return a list of row dicts for the 622 prokaryotic SMGC genomes."""
    if not os.path.isfile(path):
        die("metadata file not found: %s" % path)

    if path.lower().endswith((".xlsx", ".xls")):
        try:
            import pandas
        except ImportError:
            die(
                "reading %s needs pandas + openpyxl (pip install -r requirements.txt).\n"
                "Alternatively export sheet %r to CSV and pass that with --metadata."
                % (path, XLSX_SHEET)
            )
        # Row 1 of the sheet is the table caption; the real header is row 2.
        frame = pandas.read_excel(path, sheet_name=XLSX_SHEET, header=1)
        rows = frame.to_dict("records")
    else:
        with open(path, newline="", encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))

    if not rows:
        die("no rows read from %s" % path)
    missing = [col for col, _ in RANK_COLUMNS if col not in rows[0]]
    if "identifier" not in rows[0]:
        missing.append("identifier")
    if missing:
        die(
            "%s is missing expected column(s): %s\nColumns present: %s"
            % (path, ", ".join(missing), ", ".join(map(str, rows[0])))
        )
    return rows


def lineage_of(row):
    """GTDB lineage as a list of (rank, name), filling blank species from the id."""
    identifier = str(row["identifier"]).strip()
    lineage = []
    for column, rank in RANK_COLUMNS:
        value = row.get(column)
        value = "" if value is None else str(value).strip()
        if not value or value in {"nan", "NaN"}:
            value = "%s__" % rank[0]
        if rank == "species" and value.strip() in {"s__", ""}:
            # SMGC leaves gtdb_species empty for novel species; make the leaf unique
            # by falling back to the genome identifier, exactly as the reference
            # MycoProfiler database build did.
            value = "s__%s" % identifier
        lineage.append((rank, value))
    return identifier, lineage


class Taxonomy(object):
    """Incrementally built NCBI-style taxonomy dump."""

    def __init__(self):
        self.next_taxid = 2
        self.nodes = {1: (1, "no rank")}      # taxid -> (parent, rank)
        self.names = {1: "root"}
        self._by_path = {(): 1}

    def add_lineage(self, lineage):
        path = ()
        parent = 1
        for rank, name in lineage:
            path = path + ((rank, name),)
            if path in self._by_path:
                parent = self._by_path[path]
                continue
            taxid = self.next_taxid
            self.next_taxid += 1
            self.nodes[taxid] = (parent, rank)
            self.names[taxid] = name
            self._by_path[path] = taxid
            parent = taxid
        return parent  # the species-level leaf

    def write(self, directory):
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, "nodes.dmp"), "w") as handle:
            for taxid in sorted(self.nodes):
                parent, rank = self.nodes[taxid]
                handle.write("%d\t|\t%d\t|\t%s\t|\t-\t|\n" % (taxid, parent, rank))
        with open(os.path.join(directory, "names.dmp"), "w") as handle:
            for taxid in sorted(self.names):
                handle.write(
                    "%d\t|\t%s\t|\t\t|\tscientific name\t|\n" % (taxid, self.names[taxid])
                )


def write_library(catalogue, assignments, library_path):
    """Concatenate the genomes, tagging every header with its species taxid."""
    written_genomes = 0
    written_sequences = 0
    missing = []
    with open(library_path, "w") as out:
        for identifier, taxid in sorted(assignments.items()):
            source = None
            for extension in (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz"):
                candidate = os.path.join(catalogue, identifier + extension)
                if os.path.isfile(candidate):
                    source = candidate
                    break
            if source is None:
                missing.append(identifier)
                continue
            opener = open
            mode = "r"
            if source.endswith(".gz"):
                import gzip
                opener, mode = gzip.open, "rt"
            with opener(source, mode) as handle:
                for line in handle:
                    if line.startswith(">"):
                        seqid = line[1:].strip().split()[0]
                        out.write(">%s|kraken:taxid|%d\n" % (seqid, taxid))
                        written_sequences += 1
                    else:
                        out.write(line)
            written_genomes += 1
    return written_genomes, written_sequences, missing


def run(command, description):
    sys.stderr.write("\n== %s\n$ %s\n" % (description, " ".join(command)))
    result = subprocess.run(command)
    if result.returncode != 0:
        die("%s failed (exit status %d): %s" % (description, result.returncode, " ".join(command)))


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--catalogue", required=True,
                        help="extracted SMGC 03_finalcatalogue/ directory")
    parser.add_argument("--metadata", required=True,
                        help="SMGC.xlsx, or a CSV export of its 'Supplementary Table 4'")
    parser.add_argument("--output", required=True, help="Kraken2 database directory to create")
    parser.add_argument("--threads", type=int, default=8, help="threads for kraken2-build")
    parser.add_argument("--kraken2-build", default="kraken2-build",
                        help="kraken2-build executable (default: from $PATH)")
    parser.add_argument("--mask-low-complexity", action="store_true",
                        help="run kraken2-build's dustmasker step (needs BLAST+ dustmasker). "
                             "Off by default, matching how the reference MycoProfiler "
                             "databases were built.")
    parser.add_argument("--keep-library", action="store_true",
                        help="keep library/ and taxonomy/ after a successful build "
                             "(default: removed, saving tens of GB)")
    parser.add_argument("--force", action="store_true",
                        help="back up recognized existing DB files and build in a fresh directory")
    args = parser.parse_args(argv)

    builder = shutil.which(args.kraken2_build)
    if not builder:
        die("'%s' not found on $PATH. Install Kraken2 (see environment.yml)."
            % args.kraken2_build)
    if not os.path.isdir(args.catalogue):
        die("catalogue directory not found: %s" % args.catalogue)

    if os.path.exists(args.output) and os.listdir(args.output) and not args.force:
        die("output directory %s is not empty; pass --force to overwrite" % args.output)
    if args.threads < 1:
        die("--threads must be >= 1")
    if args.force and os.path.isdir(args.output) and os.listdir(args.output):
        if os.path.islink(args.output):
            die("--force refuses a symlink output directory")
        allowed = {"taxonomy", "library", "hash.k2d", "opts.k2d", "taxo.k2d",
                   "seqid2taxid.map", "smgc_prokaryotes.tagged.fna", "build_provenance.json"}
        if set(os.listdir(args.output)) - allowed:
            die("output contains unrecognised files; choose a new --output directory")
        backup = args.output.rstrip(os.sep) + ".backup-" + datetime.datetime.now().strftime("%Y%m%dT%H%M%S%f")
        os.rename(args.output, backup)
        sys.stderr.write("Previous database preserved at %s\n" % backup)
    os.makedirs(args.output, exist_ok=True)

    rows = read_metadata(args.metadata)
    sys.stderr.write("read %d genome records from %s\n" % (len(rows), args.metadata))

    taxonomy = Taxonomy()
    assignments = {}
    for row in rows:
        identifier, lineage = lineage_of(row)
        assignments[identifier] = taxonomy.add_lineage(lineage)
    sys.stderr.write(
        "built taxonomy: %d nodes, %d species leaves\n"
        % (len(taxonomy.nodes), len(set(assignments.values())))
    )

    taxonomy_dir = os.path.join(args.output, "taxonomy")
    taxonomy.write(taxonomy_dir)

    library_path = os.path.join(args.output, "smgc_prokaryotes.tagged.fna")
    genomes, sequences, missing = write_library(args.catalogue, assignments, library_path)
    sys.stderr.write(
        "wrote library: %d genomes, %d sequences -> %s\n" % (genomes, sequences, library_path)
    )
    if missing:
        die("%d genome(s) missing from catalogue: %s" % (len(missing), ", ".join(missing[:5])))
    if genomes == 0:
        die("no genome FASTA files matched the metadata identifiers under %s" % args.catalogue)

    add_command = [builder, "--add-to-library", library_path, "--db", args.output,
                   "--threads", str(args.threads)]
    if not args.mask_low_complexity:
        add_command.append("--no-masking")
    run(add_command, "kraken2-build --add-to-library")

    build_command = [builder, "--build", "--db", args.output, "--threads", str(args.threads)]
    run(build_command, "kraken2-build --build")

    # A Kraken2 build can exit 0 and still leave an EMPTY hash table if the taxonomy
    # and the library headers disagree. Catch that here rather than at classify time.
    inspect = shutil.which("kraken2-inspect")
    if inspect:
        result = subprocess.run(
            [inspect, "--db", args.output], stdout=subprocess.PIPE, universal_newlines=True
        )
        if result.returncode != 0:
            die("kraken2-inspect failed; database is not validated")
        table_size = None
        for line in (result.stdout or "").splitlines():
            if line.startswith("# Table size:"):
                table_size = int(line.split(":", 1)[1].strip())
                break
        if table_size is None:
            for line in (result.stdout or "").splitlines():
                fields = line.split("\t")
                if len(fields) >= 5 and fields[3].strip() == "R" and fields[4].strip() == "1":
                    table_size = int(fields[1].strip())
                    break
        if table_size is None or table_size <= 0:
            die(
                "kraken2-inspect reports 'Table size: 0' -- the database is empty even "
                "though the build exited cleanly. This means the taxonomy did not match "
                "the library headers. Do not use this database."
            )
        sys.stderr.write("kraken2-inspect: Table size = %s\n" % table_size)
    else:
        die("kraken2-inspect is required to validate the database")

    for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        if not os.path.isfile(os.path.join(args.output, name)):
            die("build finished but %s is missing from %s" % (name, args.output))

    with open(os.path.join(args.output, "build_provenance.json"), "w") as handle:
        json.dump({"genomes": assignments, "genome_count": genomes,
                   "sequence_count": sequences, "mask_low_complexity": args.mask_low_complexity,
                   "add_command": add_command, "build_command": build_command,
                   "taxonomy_names": taxonomy.names}, handle, indent=2)

    if not args.keep_library:
        os.remove(library_path)
        shutil.rmtree(os.path.join(args.output, "library"), ignore_errors=True)
        shutil.rmtree(taxonomy_dir, ignore_errors=True)
        sys.stderr.write("removed library/ and taxonomy/ (pass --keep-library to retain)\n")

    sys.stderr.write("\nSkin Kraken2 database ready: %s\n" % args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
