# MycoProfiler

**Body-site-specific fungal profiling of human shotgun metagenomes.**

MycoProfiler runs two Kraken2 passes over one sample. The first classifies the reads
against the *prokaryotic genome catalogue* of the body site the sample came from and
sets aside everything Kraken2 leaves **unclassified**. The second classifies exactly
those set-aside reads against the *MycobiomeDB fungal database* for the same body
site. Both passes run at `--confidence 0.2`.

```
                        ┌──────────────────────────────┐
  input reads  ────────►│  stage 1   Kraken2           │
  (FASTQ, PE or SE)     │  vs the site's prokaryotic   │──► classified  → report only
                        │  catalogue (HRGM2 / HROM /   │
                        │  SMGC / VMGC), confidence 0.2│──► UNCLASSIFIED
                        └──────────────────────────────┘         │
                                                                 │  (FASTQ, pairing preserved)
                                                                 ▼
                        ┌──────────────────────────────┐
                        │  stage 2   Kraken2           │
                        │  vs MycobiomeDB for the same │──► fungal report + per-read output
                        │  site, confidence 0.2        │
                        └──────────────────────────────┘
```

**Kraken2 is a k-mer classifier, not an aligner.** It breaks each read into k-mers,
looks them up in a hash table built from the reference genomes, and assigns the read
to the lowest common ancestor of the taxa those k-mers hit. Nothing in MycoProfiler
performs sequence alignment, and "unmapped" is not a word that applies here — the
reads carried from stage 1 to stage 2 are the ones Kraken2 reported as
**unclassified**.

Scope of this tool: **decontamination and fungal classification only.** It stops at
the stage-2 Kraken2 report. Abundance re-estimation (Bracken), normalisation,
diversity and differential-abundance analysis are deliberately out of scope.

---

## Contents

1. [Supported body sites](#1-supported-body-sites)
2. [Requirements and installation](#2-requirements-and-installation)
3. [Getting the databases from nothing](#3-getting-the-databases-from-nothing)
4. [Database directory layout and configuration](#4-database-directory-layout-and-configuration)
5. [Running MycoProfiler](#5-running-mycoprofiler)
6. [Command-line reference](#6-command-line-reference)
7. [Output files](#7-output-files)
8. [Confidence threshold](#8-confidence-threshold)
9. [Paired-end handling and counting units](#9-paired-end-handling-and-counting-units)
10. [Intermediate files](#10-intermediate-files)
11. [Troubleshooting](#11-troubleshooting)
12. [Scope and interpretation](#12-scope-and-interpretation)
13. [Citation, licences and terms](#14-citation-licences-and-terms)

---

## 1. Supported body sites

Choosing a site with `--site` selects **both** databases. They are never mixed
across sites.

| Body site | Stage 1: prokaryotic decontamination | Stage 2: fungal classification | Fungal species |
|---|---|---|---:|
| `gut` | HRGM2 representatives (4,824 species) | **Myco-HG** (MycobiomeDB Gut) | 766 |
| `oral` | HROM representatives (3,426 species) | **Myco-HO** (MycobiomeDB Oral) | 259 |
| `skin` | SMGC prokaryotic representatives (622 species) | **Myco-HS** (MycobiomeDB Skin) | 642 |
| `vagina` | VMGC prokaryotic SGB representatives (786 species) | **Myco-HV** (MycobiomeDB Vaginal) | 112 |

Every stage-1 catalogue is **prokaryote-only**: bacteria plus a small number of
archaea (26 in HRGM2, 2 in HROM, 1 in SMGC, none in VMGC). None contains human
sequence, fungi, or viruses. Full provenance, citations and composition:
[`docs/DATABASES.md`](docs/DATABASES.md).

Pick the database that matches where the sample came from — feces or intestinal
biopsy → `gut`; saliva, plaque, tongue or oral rinse → `oral`; skin swab → `skin`;
vaginal swab → `vagina`. For a sample type that matches none of the four, consider
carefully whether any of these databases is biologically appropriate.

---

## 2. Requirements and installation

### System requirements

| Item | Requirement |
|---|---|
| OS | Linux x86-64 (developed and tested there) |
| Python | 3.8 or newer |
| External tool | **Kraken2 ≥ 2.1.3** (`kraken2`; `kraken2-build` only if you build the skin database) |
| Python packages | none for the pipeline itself — standard library only |
| RAM | roughly the size of the database's `hash.k2d`, one stage at a time — gut is the largest at ~35 GB (stage 1) and ~40 GB (stage 2). Use `--memory-mapping` to trade speed for a few GB instead. |
| Disk | see the sizing table in [`docs/DATABASES.md` §3](docs/DATABASES.md) |

### Install

```bash
git clone https://github.com/netbiolab/MycobiomeDB.git
cd MycobiomeDB/MycoProfiler

conda env create -f environment.yml     # installs Python + Kraken2
conda activate mycoprofiler

pip install -e .                        # optional: provides the `mycoprofiler` command
```

Without `pip install`, run it as a module from the repository root:

```bash
python -m mycoprofiler --help
```

Both forms are equivalent; `mycoprofiler ...` is used throughout this README.

If you already have Kraken2 from somewhere else, MycoProfiler needs nothing but a Python
3.8+ interpreter and that `kraken2` binary. Point at a specific one with
`--kraken2 /path/to/kraken2`.

### Confirm the setup

```bash
mycoprofiler check --db-root /data/mycoprofiler_db
```

`check` reports the Kraken2 it found, its version, and — for each body site — whether
both databases resolve to a directory containing `hash.k2d`, `opts.k2d` and
`taxo.k2d`. It classifies nothing and exits non-zero if anything is missing.

---

## 3. Getting the databases from nothing

This section assumes you have **no databases and no reference genomes**. Two helper
scripts do most of the work. Both are resumable (`curl -C -`), skip sites that are
already installed, and verify what they download.

> Before you start, check you have the disk space. Gut alone is ~35 GB (stage 1) plus
> ~40 GB (stage 2). All four sites for both stages is roughly **170 GB** of final
> database, with a transient extra ~40 GB while the largest archive is unpacked.

### 3.1 Fungal databases (stage 2) — fully automatic

MycobiomeDB publishes **prebuilt Kraken2 databases**, so nothing has to be built.

```bash
# all four sites
scripts/download_fungal_db.sh --db-root /data/mycoprofiler_db

# or just what you need
scripts/download_fungal_db.sh --db-root /data/mycoprofiler_db --site gut --site vagina
```

The script downloads from Zenodo record `21387880` (MycobiomeDB v1.0.0, DOI
`10.5281/zenodo.21387880`), checks the exact byte count, verifies the **MD5** sum
published in that record, extracts the archive, locates the directory containing
`hash.k2d`, and installs it as `<db-root>/fungal/<site>/`. Genome FASTA and metadata
shipped alongside the database are kept in `<db-root>/fungal/<site>_source/`.

Compressed download sizes: gut 37.9 GiB, oral 10.6 GiB, skin 25.2 GiB, vagina 5.1 GiB.

> The MycobiomeDB README mentions `checksums.sha256`; the Zenodo record actually
> publishes `checksums.md5`. MycoProfiler verifies the MD5 sums that are really there.

### 3.2 Prokaryotic databases (stage 1) — automatic for three sites, one build for skin

```bash
scripts/download_bacterial_db.sh --db-root /data/mycoprofiler_db
```

| Site | Distributed as | What the script does |
|---|---|---|
| `gut` | prebuilt Kraken2 DB | downloads `hash.k2d` / `opts.k2d` / `taxo.k2d` — ready to use |
| `oral` | prebuilt Kraken2 DB | downloads the same three files — ready to use |
| `vagina` | prebuilt Kraken2 DB (`.tar.gz`) | downloads and extracts — ready to use |
| `skin` | **genome FASTA only** | downloads the SMGC catalogue + metadata, then you run one build command |

**Skin needs one extra step**, because SMGC does not publish a Kraken2 database:

```bash
python scripts/build_smgc_kraken2_db.py \
    --catalogue /data/mycoprofiler_db/sources/smgc/03_finalcatalogue \
    --metadata  /data/mycoprofiler_db/sources/smgc/SMGC.xlsx \
    --output    /data/mycoprofiler_db/bacterial/skin \
    --threads 16
```

That script reads the 622-row prokaryotic supplementary table, derives a GTDB-style
Kraken2 taxonomy (one leaf per species), rewrites every contig header into Kraken2's
`<seqid>|kraken:taxid|<taxid>` form, and runs `kraken2-build --add-to-library`
followed by `--build`. Afterwards it runs `kraken2-inspect` and **refuses to report
success if the hash table came out empty** — a Kraken2 build can exit 0 and still
produce an empty database when the taxonomy and the library headers disagree.

It needs `pandas` + `openpyxl` to read the `.xlsx` (`pip install -r requirements.txt`);
if you would rather not install those, export sheet *Supplementary Table 4* to CSV
and pass the CSV to `--metadata`.

Only the **622 prokaryotic** SMGC genomes go into the database. The catalogue's 12
eukaryotic and 6,935 viral genomes are excluded on purpose — putting fungi into the
decontamination database would remove exactly the reads stage 2 is looking for.

#### HRGM2 has two published builds

Both come from the same 4,824 species representatives. MycoProfiler defaults to
`concat` (contigs concatenated per genome, `hash.k2d` 34.7 GiB), which is the build
the reference analyses used. The alternative is smaller:

```bash
scripts/download_bacterial_db.sh --db-root /data/mycoprofiler_db --site gut --hrgm2-variant rep
```

### 3.3 Verify

```bash
mycoprofiler check --db-root /data/mycoprofiler_db
```

---

## 4. Database directory layout and configuration

### Recommended layout

The download scripts produce this, and `--db-root` expects it:

```
/data/mycoprofiler_db/
├── bacterial/
│   ├── gut/           hash.k2d  opts.k2d  taxo.k2d      (HRGM2)
│   ├── oral/          hash.k2d  opts.k2d  taxo.k2d      (HROM)
│   ├── skin/          hash.k2d  opts.k2d  taxo.k2d      (SMGC, built locally)
│   └── vagina/        hash.k2d  opts.k2d  taxo.k2d      (VMGC)
├── fungal/
│   ├── gut/           hash.k2d  opts.k2d  taxo.k2d      (Myco-HG)
│   ├── oral/          …                                 (Myco-HO)
│   ├── skin/          …                                 (Myco-HS)
│   └── vagina/        …                                 (Myco-HV)
├── fungal/<site>_source/    genome FASTA + metadata from the MycobiomeDB archive
└── sources/smgc/            SMGC genomes + metadata (inputs to the skin build)
```

A Kraken2 database is the **directory** holding `hash.k2d`, `opts.k2d` and
`taxo.k2d`. Anything else a build leaves behind (`library/`, `taxonomy/`,
`database.kraken`, `database*mers.kmer_distrib`, `seqid2taxid.map`) is unnecessary
for classification and can be deleted to reclaim a lot of disk.

### Four ways to tell MycoProfiler where the databases are

Checked in this order; the first match wins.

```bash
# 1. explicit paths — for databases that live anywhere at all
mycoprofiler run --site gut --bacterial-db /any/HRGM2 --fungal-db /any/Myco-HG ...

# 2. a JSON config
cp config/databases.example.json config/databases.json   # then edit the paths
mycoprofiler run --site gut --db-config config/databases.json ...

# 3. the conventional root
mycoprofiler run --site gut --db-root /data/mycoprofiler_db ...

# 4. the environment — set once, omit from every command
export MYCOPROFILER_DB_ROOT=/data/mycoprofiler_db
# or: export MYCOPROFILER_DB_CONFIG=/abs/path/to/databases.json
mycoprofiler run --site gut -1 A_R1.fastq.gz -2 A_R2.fastq.gz -o out/A
```

The config is JSON, not YAML, so MycoProfiler has no third-party Python dependency.
Top-level keys starting with `_` are ignored, which gives JSON the comment field it
otherwise lacks. Sites and keys you do not use may be omitted entirely.

```json
{
  "_comment": "paths are directories containing hash.k2d/opts.k2d/taxo.k2d",
  "gut":    {"bacterial": "/data/mycoprofiler_db/bacterial/gut",
             "fungal":    "/data/mycoprofiler_db/fungal/gut"},
  "vagina": {"bacterial": "/data/mycoprofiler_db/bacterial/vagina",
             "fungal":    "/data/mycoprofiler_db/fungal/vagina"}
}
```

**No path is hardcoded anywhere in MycoProfiler.** Every database and input location comes
from an argument, a config file, or an environment variable.

---

## 5. Running MycoProfiler

One command runs both stages. One invocation processes one sample.

### Paired-end

```bash
mycoprofiler run \
    --site gut \
    -1 SRR1234567_1.fastq.gz \
    -2 SRR1234567_2.fastq.gz \
    --db-root /data/mycoprofiler_db \
    -o results/SRR1234567 \
    --threads 16
```

### Single-end

```bash
mycoprofiler run \
    --site gut \
    -U SRR1234567.fastq.gz \
    --db-root /data/mycoprofiler_db \
    -o results/SRR1234567 \
    --threads 16
```

Plain `.fastq` and gzipped `.fastq.gz` are both accepted, and detected from the
filename. In paired-end mode both mates must be compressed the same way; a mixed
pair is rejected up front rather than mishandled.

### One example per body site

```bash
# gut — fecal metagenome
mycoprofiler run --site gut    -1 gut_R1.fastq.gz    -2 gut_R2.fastq.gz \
            --db-root /data/mycoprofiler_db -o results/gut_sample    -t 16

# oral — saliva or dental plaque
mycoprofiler run --site oral   -1 saliva_R1.fastq.gz -2 saliva_R2.fastq.gz \
            --db-root /data/mycoprofiler_db -o results/saliva_sample -t 16

# skin — skin swab
mycoprofiler run --site skin   -1 skin_R1.fastq.gz   -2 skin_R2.fastq.gz \
            --db-root /data/mycoprofiler_db -o results/skin_sample   -t 16

# vagina — vaginal swab
mycoprofiler run --site vagina -1 vag_R1.fastq.gz    -2 vag_R2.fastq.gz \
            --db-root /data/mycoprofiler_db -o results/vaginal_sample -t 16
```

### Several samples

There is no built-in batch mode; loop in the shell, which keeps scheduling and
parallelism under your control.

```bash
export MYCOPROFILER_DB_ROOT=/data/mycoprofiler_db
for r1 in reads/*_1.fastq.gz; do
    sample=$(basename "$r1" _1.fastq.gz)
    mycoprofiler run --site gut \
        -1 "$r1" -2 "reads/${sample}_2.fastq.gz" \
        -o "results/${sample}" --sample-id "$sample" --threads 16 \
      || echo "FAILED: $sample" >> results/failures.txt
done
```

Run samples one after another rather than side by side unless you have the RAM for
several copies of the database. `--memory-mapping` lets concurrent processes share
one on-disk copy at the cost of speed.

### See the commands without running them

```bash
mycoprofiler run --site gut -1 a_1.fq.gz -2 a_2.fq.gz --db-root /data/mycoprofiler_db \
            -o results/a --dry-run
```

`--dry-run` prints the two exact Kraken2 command lines and writes nothing.

---

## 6. Command-line reference

### `mycoprofiler run`

| Argument | Default | Meaning |
|---|---|---|
| `--site {gut,oral,skin,vagina}` | *required* | body site; selects both databases |
| `-1`, `--r1 FASTQ` | — | forward mate (paired-end) |
| `-2`, `--r2 FASTQ` | — | reverse mate (paired-end) |
| `-U`, `--single FASTQ` | — | single-end reads; mutually exclusive with `-1`/`-2` |
| `-o`, `--outdir DIR` | *required* | output directory |
| `--sample-id ID` | derived from the first read filename | name used in output filenames |
| `-t`, `--threads N` | `1` | threads passed to **both** Kraken2 stages |
| `--db-root DIR` | `$MYCOPROFILER_DB_ROOT` | root of the conventional layout |
| `--db-config FILE` | `$MYCOPROFILER_DB_CONFIG` | JSON site→paths config |
| `--bacterial-db DIR` | — | explicit stage-1 database |
| `--fungal-db DIR` | — | explicit stage-2 database |
| `--kraken2 PATH` | first `kraken2` on `$PATH` | Kraken2 executable to use |
| `--memory-mapping` | off | read the hash table from disk instead of RAM |
| `--keep-intermediate` | off (they are deleted) | keep the stage-1 unclassified FASTQs |
| `--gzip-intermediate` | off | gzip the kept intermediates; implies `--keep-intermediate` |
| `--force` | off | overwrite a non-empty output directory |
| `--dry-run` | off | print the two Kraken2 commands and exit |
| `-q`, `--quiet` | off | suppress progress on stderr (`mycoprofiler.log` is still written) |

The confidence threshold is **not** a command-line option; see
[§8](#8-confidence-threshold).

### `mycoprofiler check`

| Argument | Default | Meaning |
|---|---|---|
| `--site NAME` | all four | site to check; repeatable |
| `--db-root` / `--db-config` / `--bacterial-db` / `--fungal-db` | as above | where to look |
| `--kraken2 PATH` | from `$PATH` | executable to check |

Exit codes: `0` success · `1` error (bad input, missing database, failed stage) ·
`2` bad command line · `130` interrupted.

---

## 7. Output files

```
results/SRR1234567/
├── 01_bacterial_decontamination/
│   ├── SRR1234567.bacterial.kraken2.report    # Kraken2 report, stage 1
│   ├── SRR1234567.bacterial.kraken2.output    # per-read calls, stage 1
│   └── SRR1234567.bacterial.kraken2.log       # command, stderr, exit status, runtime
├── 02_fungal_classification/
│   ├── SRR1234567.fungal.kraken2.report       # ← the fungal profile: the main result
│   ├── SRR1234567.fungal.kraken2.output       # per-read fungal calls
│   └── SRR1234567.fungal.kraken2.log
├── intermediate_unclassified_reads/           # empty unless --keep-intermediate
│   ├── SRR1234567.unclassified_1.fastq
│   └── SRR1234567.unclassified_2.fastq
├── mycoprofiler_summary.tsv                        # one row, machine-readable
├── mycoprofiler_run_manifest.json                  # everything needed to reproduce the run
└── mycoprofiler.log                                # human-readable run log
```

| File | What it is |
|---|---|
| `*.kraken2.report` | Standard Kraken2 report: percentage, clade-rooted reads, direct reads, rank code, taxid, name. Written with `--report-zero-counts`, so every taxon in the database appears even at zero. Species rows are rank code `S`. |
| `*.kraken2.output` | One line per read (or read pair): `C`/`U`, read ID, assigned taxon, length, LCA k-mer map. Written with `--use-names`, so taxon names appear alongside taxids. |
| `*.kraken2.log` | The exact command, Kraken2's own stderr, the exit status and the elapsed time for that stage. First place to look when a stage fails. |
| `mycoprofiler_summary.tsv` | Header plus one row. Every count is in the unit named by the `counting_unit` column; see [§9](#9-paired-end-handling-and-counting-units). |
| `mycoprofiler_run_manifest.json` | Both full command lines, Kraken2's version, the resolved database paths with the size and mtime of each `.k2d` file, the confidence used at each stage, hostname, platform, thread count, timestamps, and the per-stage counts. |
| `mycoprofiler.log` | Timestamped narrative of the run, identical to what is printed to stderr. |

The **stage-2 report is the fungal profile** — the primary result. Extract
species-level rows with:

```bash
awk -F'\t' '$4 == "S" && $2 > 0' \
    results/SRR1234567/02_fungal_classification/SRR1234567.fungal.kraken2.report
```

The stage-1 report is kept for QC: it tells you how much of the sample the site's
prokaryotic catalogue accounted for.

### `mycoprofiler_summary.tsv` columns

| Column | Meaning |
|---|---|
| `sample_id`, `body_site`, `library_layout` | run identity |
| `counting_unit` | `read_pairs` (paired-end) or `reads` (single-end) — **the unit for every count below** |
| `input_units` | what Kraken2 processed at stage 1, in `counting_unit` |
| `input_reads_total` | the same input expressed in individual reads (`2 ×` pairs) |
| `stage1_classified_prokaryotic`, `stage1_classified_percent` | assigned to the prokaryotic catalogue and removed |
| `stage1_unclassified`, `stage1_unclassified_percent` | carried forward to stage 2 |
| `stage2_input_units` | what stage 2 received — equals `stage1_unclassified` |
| `stage2_classified_fungal` | assigned to a MycobiomeDB taxon |
| `stage2_classified_percent_of_stage2_input` | as a fraction of what stage 2 saw |
| `stage2_classified_percent_of_input` | as a fraction of the **original** input — usually the number you want |
| `stage1_seconds`, `stage2_seconds` | wall-clock per stage |

---

### Bracken abundance estimation

Bracken can be run separately after MycoProfiler to estimate fungal abundances at the species level.

Use the complete stage-2 fungal Kraken2 report (`*.fungal.kraken2.report`). Do not use the bacterial report, per-read `.kraken2.output`, or a report filtered to species rows.

#### Choose the read length

Set `READ_LEN` to the individual read length in base pairs. For unmerged paired-end reads of 2 × 150 bp, use `150`, not `300`.

For trimmed reads with variable lengths, inspect the lengths of the reads supplied to fungal classification. A representative length is an approximation; document the selection and evaluate its effect if the distribution is broad. Do not automatically assume that the original sequencing-cycle length still applies.

#### Check the Bracken database file

Use the exact fungal database used for MycoProfiler stage 2. The following example assumes gut reads of 150 bp:

```bash
FUNGAL_DB=/data/mycoprofiler_db/fungal/gut
READ_LEN=150

ls -lh "${FUNGAL_DB}/database${READ_LEN}mers.kmer_distrib"
```

If the matching file already exists for this database version, no additional build is needed.

If it is missing, generate it using Bracken:

```bash
# Example only: verify the actual Kraken2 database k-mer length.
KMER_LEN=35

bracken-build \
  -d "${FUNGAL_DB}" \
  -t 16 \
  -k "${KMER_LEN}" \
  -l "${READ_LEN}"
```

`bracken-build` requires the matching database source sequences and taxonomy/mapping resources, including `library/`, `taxonomy/`, and `seqid2taxid.map`. The three `.k2d` files alone are insufficient for this standard build workflow. If these resources are absent, obtain the matching build resources or a precomputed Bracken distribution file from the database provider.

Build once per fungal database version and read length, then reuse the result for matching samples.

#### Run species-level abundance estimation

With Bracken installed and available on `PATH`, run:

```bash
FUNGAL_DB=/data/mycoprofiler_db/fungal/gut
READ_LEN=150
SAMPLE=SRR1234567
RESULT_DIR="results/${SAMPLE}"

mkdir -p "${RESULT_DIR}/03_bracken"

bracken \
  -d "${FUNGAL_DB}" \
  -i "${RESULT_DIR}/02_fungal_classification/${SAMPLE}.fungal.kraken2.report" \
  -o "${RESULT_DIR}/03_bracken/${SAMPLE}.fungal.bracken.tsv" \
  -r "${READ_LEN}" \
  -l S \
  -t 10
```

Here, `-r` selects the read-length-specific distribution, `-l S` requests species-level estimates, and `-t 10` sets the read-count threshold. Bracken's `-t` is not a thread option.

Adjust the site, sample, paths, read length, and threshold for your analysis. Record the Bracken version and parameters alongside the MycoProfiler outputs.

See the [Bracken manual](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual) for installation and parameter details.

---

## 8. Confidence threshold

**Both Kraken2 stages run at `--confidence 0.2`, and this is fixed.** There is no
flag to change it. It is a module-level constant (`mycoprofiler.CONFIDENCE`), it is
written into every `.kraken2.log`, and it is recorded twice in
`mycoprofiler_run_manifest.json`:

```json
"mycoprofiler": { "confidence_stage1": 0.2, "confidence_stage2": 0.2 }
```

Kraken2's confidence score is the fraction of a read's k-mers that must be consistent
with the assigned taxon. At `0.2` a read is left unclassified unless at least 20 % of
its k-mers support the call. Fixing it means every MycoProfiler result is directly
comparable to every other, and that the amount of material carried from stage 1 to
stage 2 is a property of the sample rather than of the settings someone chose.

If you need a different threshold, call Kraken2 yourself — but the result is then not
a MycoProfiler run and should not be described as one.

---

## 9. Paired-end handling and counting units

### How pairing is preserved

In `--paired` mode Kraken2 classifies a read pair **as one unit**, using the k-mers of
both mates together. MycoProfiler passes `--unclassified-out <prefix>#.fastq`; Kraken2
expands the `#` to `_1` and `_2` and writes **both mates of every unclassified pair**,
in input order, to the two files. R1/R2 correspondence is therefore maintained by
construction, and those two files are what stage 2 receives — again with `--paired`.

MycoProfiler does not take this on trust. After stage 1 it counts the records in both
intermediate files and **aborts before stage 2 if they differ**, rather than running
a de-synchronised pair. The check is recorded in the manifest:

```json
"handoff": {
  "counting_unit": "read_pairs",
  "carried_forward": 300,
  "fastq_record_counts_per_file": [300, 300],
  "mates_synchronised": true
}
```

In single-end mode there is one intermediate file and `mates_synchronised` is `null`.

### What the numbers count

This is the one place where paired-end reporting is easy to misread, so MycoProfiler
labels it explicitly rather than leaving it implicit.

> **Kraken2 counts one read PAIR as one "sequence".** A paired-end run of 10 million
> pairs (20 million reads) is reported by Kraken2 as `10000000 sequences processed`.

MycoProfiler therefore carries the unit alongside every count:

* `mycoprofiler_summary.tsv` has a `counting_unit` column — `read_pairs` or `reads`;
* every count in that row is in that unit;
* `input_reads_total` gives the same input in individual reads, so both are available
  without any conversion on your part;
* `mycoprofiler_run_manifest.json` repeats `counting_unit` inside each stage object;
* the progress log spells it out: `500 read_pairs in, 200 classified, 300 unclassified`.

The counts inside the Kraken2 `.report` and `.output` files are Kraken2's own and
follow the same convention: in paired-end mode a "read" there is a pair.

---

## 10. Intermediate files

The stage-1 unclassified reads are the hand-off between the two stages. Kraken2
writes them **uncompressed**, which for a large sample can mean several GB, so:

| Flag | Behaviour |
|---|---|
| *(default)* | deleted after stage 2 completes successfully |
| `--keep-intermediate` | kept as plain FASTQ in `intermediate_unclassified_reads/` |
| `--gzip-intermediate` | kept and gzipped (implies `--keep-intermediate`) |

Deletion happens **only after stage 2 has exited zero**. If either stage fails the
files are left in place so you can inspect them and resume manually. Whichever
happened is recorded in the manifest as `handoff.intermediate_files`
(`deleted` / `kept` / `kept_gzipped`), with the retained paths listed.

Keep them if you intend to re-run stage 2 against a different fungal database, or to
feed the same reads to another tool.

---

## 11. Troubleshooting

**`'kraken2' was not found on $PATH`**
Activate the environment (`conda activate mycoprofiler`), or pass
`--kraken2 /path/to/kraken2`. `mycoprofiler check` reports which binary is found.

**`... database for site 'gut' does not exist`**
The path was resolved but nothing is there. The message names the path *and where
that path came from* (`--db-root`, a config file, `$MYCOPROFILER_DB_ROOT`, …). Run the
download script for that site, or fix the path.

**`... is missing required Kraken2 file(s): hash.k2d`**
A partial download or an interrupted build, or you are pointing one level too high in
the tree. A Kraken2 database is the directory that directly contains `hash.k2d`,
`opts.k2d` and `taxo.k2d`; check whether they are in a subdirectory such as
`Kraken2_DB/`.

**`refusing to write into <dir>`**
The output directory is not empty. Use a different `-o`, or `--force` to overwrite.
This is deliberate: it stops a re-run from silently mixing new outputs with old.

**`mixed compression: one mate is gzip-compressed and the other is not`**
Kraken2 cannot read a mixed pair. Compress or decompress both mates.

**`step 'bacterial_decontamination' failed: kraken2 exited with status 1`**
Read the named `.log` file — it holds the full command and Kraken2's own message. The
usual causes are a truncated database, a corrupt FASTQ, or being killed by the OOM
killer. Stage 2 is never started after a stage-1 failure.

**Killed with no message / `std::bad_alloc`**
Out of memory: Kraken2 loads `hash.k2d` into RAM. Add `--memory-mapping` to read it
from disk instead, or run on a machine with more memory. Gut needs ~35 GB (stage 1)
and ~40 GB (stage 2).

**`stage 1 unclassified mates are not the same length`**
The two intermediate FASTQs disagree in record count. MycoProfiler stops rather than
running a de-synchronised pair. Check the stage-1 log, and whether the run was
interrupted or the disk filled.

**Everything is unclassified at stage 2 / the fungal report is empty**
Genuinely plausible: fungi are a very small fraction of most human metagenomes, and
many samples legitimately yield few or no fungal reads. Before concluding anything,
check `stage1_unclassified` in the summary — if it is near zero, stage 1 removed
almost everything and the sample or the site choice is the issue.

**Very slow**
`--memory-mapping` is much slower than loading the database into RAM, especially on a
network filesystem; drop it if you have the memory. Otherwise raise `--threads`.

---

## 12. Scope and interpretation

Read this before drawing conclusions from a MycoProfiler report.

**What MycoProfiler does.** Two Kraken2 classifications at confidence 0.2, and nothing
else. It stops at the stage-2 report.

**What "unclassified" means, precisely.** The reads carried from stage 1 to stage 2
are the reads Kraken2 did not assign to *the site's prokaryotic catalogue* at
confidence 0.2. That is all. In particular:

* **They are not "fungal reads".** The set still contains host reads, viral reads,
  reads from prokaryotes not represented in the catalogue, low-complexity and adapter
  sequence, and artefacts. Which of them are fungal is what stage 2 decides — and
  stage 2's answer is itself limited by what is in MycobiomeDB.
* **Stage 1 does not remove all bacterial reads.** Each catalogue holds one
  representative genome per species (HRGM2 4,824, HROM 3,426, SMGC 622, VMGC 786).
  Bacteria absent from the catalogue, and sufficiently divergent strains, are not
  classified and pass into stage 2. Decontamination reduces prokaryotic background;
  it does not eliminate it.
* **Archaea are removed too**, in the three catalogues that contain any (26 genomes
  in HRGM2, 2 in HROM, 1 in SMGC). "Bacterial decontamination" is the conventional
  name for this step, but the catalogues are prokaryotic.

**Host reads are not removed.** No stage-1 catalogue contains human sequence. If your
reads are not already host-depleted, human reads will reach stage 2. Remove host
reads upstream (for example with Bowtie2 against GRCh38 or T2T-CHM13) if that matters
for your study.

**A classification is not a detection, and a detection is not a colonisation.** As the
MycobiomeDB documentation itself notes: fungal DNA is often present at very low
abundance; closely related fungal species share genomic regions; reagent and
environmental contamination is a real concern in low-biomass samples; and finding a
taxon does not show it was viable or actively colonising. Low-count species-level
calls in particular deserve scepticism, orthogonal confirmation, and an explicit
filtering threshold chosen for your study design.

**Counts are not abundances.** The Kraken2 report gives read counts, unnormalised for
genome length or for sequencing depth. Cross-sample comparison requires normalisation
that MycoProfiler deliberately does not perform.

**The database bounds the answer.** A body-site database is a deliberate restriction
of the search space to species with evidence from that niche. That is what suppresses
false positives — and it also means a species outside the panel cannot be reported,
whatever is in the sample. Record the MycobiomeDB version you used.

---

## 13. Citation, licences and terms

MycoProfiler is a runner for other people's tools and databases. **Cite them.**

### Databases

> **MycobiomeDB** — Sehun Ahn and Insuk Lee. *MycobiomeDB: Body-site-stratified fungal
> genome database enhancing human mycobiome analysis.* Version 1.0.0. Zenodo.
> [10.5281/zenodo.21387880](https://doi.org/10.5281/zenodo.21387880) · CC BY 4.0 ·
> <https://github.com/netbiolab/MycobiomeDB>
> Record the version you used; MycobiomeDB uses semantic versioning.

> **HRGM2** (gut) — *A human gut metagenome-assembled genome catalogue spanning 41
> countries supports genome-scale metabolic models.* Nature Microbiology (2026)
> 11(1):317-334. [10.1038/s41564-025-02206-1](https://doi.org/10.1038/s41564-025-02206-1)
> · <https://www.decodebiome.org/HRGM2/> · CC BY 4.0

> **HROM** (oral) — *A high-quality genomic catalog of the human oral microbiome
> broadens its phylogeny and clinical insights.* Cell Host & Microbe (2025)
> 33(11):1977-1994.e8.
> [10.1016/j.chom.2025.10.001](https://doi.org/10.1016/j.chom.2025.10.001) ·
> <https://www.decodebiome.org/HROM/>

> **SMGC** (skin) — Saheb Kashaf S. *et al.* *Integrating cultivation and metagenomics
> for a multi-kingdom view of skin microbiome diversity and functions.* Nature
> Microbiology (2022) 7:169-179.
> [10.1038/s41564-021-01011-w](https://doi.org/10.1038/s41564-021-01011-w) ·
> <https://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/skin_microbiome/>

> **VMGC** (vagina) — *A multi-kingdom collection of 33,804 reference genomes for the
> human vaginal microbiome.* Nature Microbiology (2024) 9(8):2185-2200.
> [10.1038/s41564-024-01751-5](https://doi.org/10.1038/s41564-024-01751-5) ·
> [10.5281/zenodo.10457006](https://doi.org/10.5281/zenodo.10457006) ·
> <https://github.com/RChGO/VMGC>

### Software

> **Kraken2** — Wood D.E., Lu J., Langmead B. *Improved metagenomic analysis with
> Kraken 2.* Genome Biology (2019) 20:257.
> [10.1186/s13059-019-1891-0](https://doi.org/10.1186/s13059-019-1891-0) · MIT licence
> · <https://github.com/DerrickWood/kraken2>

### Licences of the resources MycoProfiler uses

Each is distributed by its own provider under its own terms; MycoProfiler downloads them
and does not relicense them. Per-resource terms and the outstanding questions are in
[`docs/DATABASES.md` §4](docs/DATABASES.md). Source genomes remain subject to NCBI and
JGI policies.

### MycoProfiler's own licence and citation

This project is licensed under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/deed.en).

MycobiomeDB and MycoProfiler have not yet been described in a peer-reviewed publication or preprint.

Until a corresponding publication becomes available, please cite the database using its version-specific DOI:

> `Sehun Ahn and Insuk Lee. MycobiomeDB: Body-site-stratified fungal genome database enhancing human mycobiome analysis. Version 1.0.0. Zenodo. 10.5281/zenodo.21387880`

BibTeX:

```bibtex
@dataset{mycobiomedb_v1,
  author    = {Sehun Ahn and Insuk Lee},
  title     = {MycobiomeDB: Human Body Site-Specific Fungal Genome Database},
  year      = {2026},
  version   = {1.0.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.21387880},
  url       = {https://doi.org/10.5281/zenodo.21387880}
}
```

After publication of the associated manuscript, users should cite both the manuscript and the specific database version used in their analysis.

---
