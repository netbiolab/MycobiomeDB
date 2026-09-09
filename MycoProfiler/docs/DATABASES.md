# Databases: what they are, where they come from, and what they are not

MycoProfiler needs **two** Kraken2 databases per body site. This document records where
each one comes from, what it actually contains, and which facts were verified
against the provider versus which are still open.

Everything marked **verified** was confirmed by an HTTP request to the provider
(Zenodo REST API, an HTTP `HEAD` against the file, or the provider's own directory
listing / README) while preparing this repository. Anything not verified is marked
**needs confirmation** and is not asserted as fact.

---

## 1. Stage 2 -- fungal databases (MycobiomeDB)

All four are published together in one Zenodo record.

| Field | Value |
|---|---|
| Resource | MycobiomeDB: Human Body Site-Specific Fungal Genome Database |
| Version | 1.0.0 |
| Version DOI | `10.5281/zenodo.21387880` |
| Concept DOI (always latest) | `10.5281/zenodo.21387879` |
| Creators | Ahn, Sehun; Lee, Insuk |
| License | Creative Commons Attribution 4.0 International (CC BY 4.0) |
| Code repository | <https://github.com/netbiolab/MycobiomeDB> (tag `v1.0.0`) |

**verified** — record metadata, file list, byte sizes and MD5 sums read from
`https://zenodo.org/api/records/21387880`.

| Site | MycoProfiler name | Archive | Bytes (compressed) | MD5 |
|---|---|---|---:|---|
| Gut | Myco-HG | `Gut.tar.gz` | 40,698,406,476 | `401e5a4322c3bf25f14be0727e5c7f9f` |
| Oral | Myco-HO | `Oral.tar.gz` | 11,398,360,153 | `3284577302089ce0786c88216a10f121` |
| Skin | Myco-HS | `Skin.tar.gz` | 27,052,137,193 | `08f327a8cace8702c2beee969f1ceea6` |
| Vagina | Myco-HV | `Vagina.tar.gz` | 5,438,786,935 | `e7808ba3e3129d16752b2e1f5d35eb3a` |

Download URL pattern (verified):
`https://zenodo.org/api/records/21387880/files/<Archive>/content`

Species counts, from the MycobiomeDB README: **gut 766, oral 259, skin 642,
vagina 112**, curated from a shared pool of 4,094 fungal species (16,365 genomes,
completeness > 90 %, contamination < 5 %, from NCBI and JGI).

### Distribution form

The MycobiomeDB README states each site package contains both `Genomes/` (species
representative FASTA + metadata) and `Kraken2_DB/` (a **prebuilt Kraken2 database**),
and adds that "the exact directory structure may differ depending on the release
package". So no extra build step should be required for stage 2.

> **needs confirmation** — the internal layout of the four `.tar.gz` archives has
> not been inspected here; that would mean downloading up to 40 GB, which was not
> done. `scripts/download_fungal_db.sh` therefore does **not** hardcode
> `<Site>/Kraken2_DB/`: it searches the extracted tree for `hash.k2d` and installs
> the directory that contains it, failing loudly with instructions if no such
> directory exists.

> **discrepancy to be aware of** — the MycobiomeDB README tells users to verify
> downloads with `sha256sum -c checksums.sha256`, but the Zenodo record publishes
> `checksums.md5` (185 bytes) containing MD5 sums. MycoProfiler verifies the **MD5**
> sums above, which are the ones actually published.

---

## 2. Stage 1 -- prokaryotic decontamination databases

### 2.1 What these databases contain

This matters for interpretation, so it is stated precisely. Each stage-1 catalogue
is a **prokaryotic genome collection**: bacteria plus a handful of archaea. None of
them contains human/host sequence, fungi, or viruses.

| Site | Catalogue | Genomes in the Kraken2 DB | Composition |
|---|---|---:|---|
| Gut | HRGM2 | 4,824 species representatives | 4,798 bacteria + 26 archaea |
| Oral | HROM | 3,426 species representatives | 3,424 bacteria + 2 archaea |
| Skin | SMGC | 622 prokaryotic representatives | 621 bacteria + 1 archaeon |
| Vagina | VMGC | 786 prokaryotic SGB representatives | bacteria only |

**verified** — counts tallied from each catalogue's own genome-to-taxonomy mapping
table, and cross-checked against the species counts reported in the corresponding
publications (HRGM2 4,824; HROM 3,426; SMGC 622; VMGC 786).

Two consequences follow directly, and are repeated in the README:

* **Reads left unclassified by stage 1 are not "fungal reads".** They are reads
  that did not match *this catalogue* at confidence 0.2. That set still contains
  host reads, viral reads, reads from prokaryotic species the catalogue does not
  represent, low-complexity and adapter sequence, and sequencing artefacts.
* **Stage 1 does not remove all bacterial reads.** It removes reads assigned to the
  site's representative prokaryotic catalogue. Bacteria absent from the catalogue,
  and divergent strains, survive into stage 2 and can only be excluded there by
  the fungal database's own specificity.

Because none of the four catalogues contains human sequence, **host removal is not
performed by MycoProfiler** and must be done upstream if your reads are not already
host-depleted.

### 2.2 Gut -- HRGM2

* Paper: *A human gut metagenome-assembled genome catalogue spanning 41 countries
  supports genome-scale metabolic models.* Nature Microbiology (2026) 11(1):317-334.
  DOI [`10.1038/s41564-025-02206-1`](https://doi.org/10.1038/s41564-025-02206-1),
  PMID 41345261.
* Portal: <https://www.decodebiome.org/HRGM2/> · Code: <https://github.com/netbiolab/HRGM2>
* License stated on the portal: CC BY 4.0.
* Catalogue: 155,211 non-redundant near-complete genomes representing **4,824 species**.

**A ready-to-use Kraken2 database is published**, in two builds from the same 4,824
species representatives (**verified** by HTTP `HEAD`):

| Build | `hash.k2d` bytes | Base URL |
|---|---:|---|
| `HRGMv2_Concat` (MycoProfiler default) | 37,239,042,956 | `…/HRGMv2_kraken2_customdb/HRGMv2_Concat/` |
| `HRGMv2_Rep` | 17,791,303,712 | `…/HRGMv2_kraken2_customdb/HRGMv2_Rep/` |

Base path:
`https://www.decodebiome.org/HRGM2/data/genome_catalog/Taxonomy_Profiling/HRGMv2_kraken2_customdb/`
Each build directory holds `hash.k2d`, `opts.k2d` (64 B), `taxo.k2d` (912,187 B),
plus `seqid2taxid.map`, `library/`, `taxonomy/` and Bracken files, none of which
Kraken2 needs at classification time.

> **needs confirmation from the data producer** — "Concat" and "Rep" are both built
> from the species representatives; the difference is that *Concat* concatenates each
> genome's contigs into a single sequence per genome. MycoProfiler defaults to **Concat**
> because that is the build the reference analyses this pipeline reproduces were run
> against. If your intent is specifically the non-concatenated representative build,
> pass `--hrgm2-variant rep`. The naming of the catalogue itself also varies between
> sources ("HRGM2" on the portal, "HRGMv2" in file paths); they refer to the same
> version-2 catalogue.

### 2.3 Oral -- HROM

* Paper: *A high-quality genomic catalog of the human oral microbiome broadens its
  phylogeny and clinical insights.* Cell Host & Microbe (2025) 33(11):1977-1994.e8.
  DOI [`10.1016/j.chom.2025.10.001`](https://doi.org/10.1016/j.chom.2025.10.001),
  PMID 41167188.
* Portal: <https://www.decodebiome.org/HROM/>
* Catalogue: 72,641 near-complete genomes from **3,426 species**.

**A ready-to-use Kraken2 database is published** (**verified** by HTTP `HEAD`):
`https://www.decodebiome.org/HROM/data/genome_catalog/HROM_kraken2_customdb/`
with `hash.k2d` = 12,503,236,052 B, `opts.k2d` = 64 B, `taxo.k2d` = 337,240 B.
Species-representative FASTA are separately available under
`…/genome_catalog/HROM_representative_genomes/`.

> **needs confirmation** — the HROM portal does not state an explicit reuse licence
> in the directory listing. Check the paper and the portal before redistributing.

### 2.4 Skin -- SMGC

* Paper: Saheb Kashaf S. *et al.* *Integrating cultivation and metagenomics for a
  multi-kingdom view of skin microbiome diversity and functions.* Nature Microbiology
  (2022) 7:169-179. DOI
  [`10.1038/s41564-021-01011-w`](https://doi.org/10.1038/s41564-021-01011-w),
  PMID 34952941, PMCID PMC8732310.
* Data (from the paper's Data Availability statement, **verified** by fetching the
  directory and its README):
  <https://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/skin_microbiome/>
* Code: <https://github.com/Finn-Lab/skin_microbiome>
* Metagenome reads: SRA study accession SRP002480 (plus SRP314237 for new samples).

Relevant files in that directory:

| File | Size | Contents |
|---|---:|---|
| `03_finalcatalogue.tar.gz` | 558 MB | the SMGC itself: **622 prokaryotic** + 12 eukaryotic + 6,935 viral genomes, as `SMGC_<n>.fa` |
| `SMGC.xlsx` | 494 KB | supplementary tables; sheet *Supplementary Table 4* is the 622-row prokaryotic table with GTDB lineages |
| `01_allprok/` | — | all 7,535 non-redundant prokaryotic genomes (superset; **not** what MycoProfiler uses) |
| `README.txt` | 2.1 KB | the provider's own description |

**SMGC does NOT publish a prebuilt Kraken2 database** — this is the one site where
a local build is required. `scripts/build_smgc_kraken2_db.py` does it: it reads the
622-row table, derives a GTDB-style taxonomy (`names.dmp`/`nodes.dmp`, one leaf per
species), tags every contig header as `<seqid>|kraken:taxid|<taxid>`, and runs
`kraken2-build --add-to-library` then `--build`.

**verified locally** — the script's metadata parsing and taxonomy construction were
run against the real `SMGC.xlsx`: 622 genome records, 1,029 taxonomy nodes, 619
distinct species leaves (three genomes share a species name with another, which is
expected and harmless). The `kraken2-build` steps themselves were **not** run here.

Note the difference from MGnify's *Human Skin v1.0* catalogue (579 species-level
clusters): that is a different dereplication of overlapping data and is **not** the
SMGC used here.

Two build details worth knowing:

* The script defaults to `--no-masking`, avoiding a hard dependency on BLAST+
  `dustmasker`. Pass `--mask-low-complexity` to enable masking.
  > **needs confirmation** — whether the published HRGM2/HROM/VMGC Kraken2 databases
  > were built with or without low-complexity masking is not stated by their
  > providers, so the skin database may differ from the other three in this respect.
* A Kraken2 build can exit 0 and still produce an **empty** hash table when the
  taxonomy and library headers disagree. The script runs `kraken2-inspect` afterwards
  and refuses to report success if `Table size` is 0.

### 2.5 Vagina -- VMGC

* Paper: *A multi-kingdom collection of 33,804 reference genomes for the human
  vaginal microbiome.* Nature Microbiology (2024) 9(8):2185-2200. DOI
  [`10.1038/s41564-024-01751-5`](https://doi.org/10.1038/s41564-024-01751-5),
  PMID 38907008.
* Data: Zenodo record `10.5281/zenodo.10457006` · Code: <https://github.com/RChGO/VMGC>
* Catalogue: 33,804 genomes spanning **786 prokaryotic species**, 11 fungal species
  and 4,263 viral OTUs.

**A ready-to-use Kraken2 database is published** for the prokaryotic SGBs
(**verified** via the Zenodo API and by listing the archive's first entries):

| File | Bytes | Contents |
|---|---:|---|
| `VMGC_prokaryote_SGB_KrakenDB.tar.gz` | 2,199,401,770 | Kraken2 DB: `hash.k2d`, `opts.k2d`, `taxo.k2d`, `taxonomy/`, `library/` |
| `VMGC_prokaryote_SGB.tar.gz` | 460,760,000 approx. | the 786 SGB representative genomes |
| `VMGC_eukaryocyte.tar.gz` | 201,200,000 approx. | 38 fungal genomes -- **not used**; stage 2 uses Myco-HV |

URL pattern: `https://zenodo.org/api/records/10457006/files/<file>/content`

MycoProfiler uses **only** `VMGC_prokaryote_SGB_KrakenDB`. The eukaryotic part of VMGC is
deliberately excluded from stage 1: putting fungi into the decontamination database
would remove exactly the reads stage 2 is meant to find.

---

## 3. Size, memory and time

Compressed download sizes and `hash.k2d` sizes below are **verified** byte counts
from the providers. Everything under "RAM" and "build time" is an **estimate** and
is labelled as such; none of it was measured for this repository.

| Site | Stage 1 download | Stage 1 `hash.k2d` | Stage 2 download | Stage 2 DB (approx.) |
|---|---:|---:|---:|---:|
| Gut | 34.7 GiB (concat) / 16.6 GiB (rep) | 34.7 / 16.6 GiB | 37.9 GiB | ~40 GiB |
| Oral | 11.6 GiB | 11.6 GiB | 10.6 GiB | ~12 GiB |
| Skin | 0.6 GiB of genomes, then build | built locally | 25.2 GiB | ~26 GiB |
| Vagina | 2.0 GiB | ~2.0 GiB | 5.1 GiB | ~5.3 GiB |

Stage-2 "DB (approx.)" is inferred from the compressed archive size and is not a
measured figure.

* **RAM (estimate).** Without `--memory-mapping`, Kraken2 loads `hash.k2d` into RAM,
  so a stage needs roughly the size of that file plus a small overhead. Gut is
  therefore the demanding site (~35 GB for stage 1, plus ~40 GB for stage 2, run one
  after the other, not simultaneously). With `--memory-mapping` the hash table is
  read from disk instead: RAM drops to a few GB, at the cost of speed, and it is
  slow on a busy or network filesystem.
* **Disk (estimate).** Allow roughly 2.2x the compressed archive size per site while
  extracting, since the archive and its contents briefly coexist. The download
  scripts delete each archive as soon as extraction succeeds.
* **Skin build time (estimate).** Not measured. `kraken2-build` on ~622 bacterial
  genomes is a modest job by Kraken2 standards; budget on the order of an hour with
  16 threads and re-check with your own hardware.
* **Intermediate reads (calculable, not estimated).** Kraken2 writes
  `--unclassified-out` **uncompressed**, so the intermediate FASTQs are about
  `(1 - stage-1 classified fraction)` of your input, expanded from gzip. For a 5 GB
  gzipped pair with 80 % of reads removed at stage 1, expect roughly 5 GB of plain
  FASTQ. MycoProfiler deletes these after stage 2 succeeds unless you pass
  `--keep-intermediate`.

---

## 4. Licences and terms

Each database is redistributed by its own provider under its own terms. MycoProfiler
downloads them; it does not relicense them.

| Resource | Licence / terms |
|---|---|
| MycobiomeDB | CC BY 4.0 (<https://creativecommons.org/licenses/by/4.0/>). Source genomes remain subject to NCBI and JGI terms. |
| HRGM2 | CC BY 4.0, as stated on <https://www.decodebiome.org/HRGM2/> |
| HROM | **needs confirmation** — no licence stated in the portal listing |
| SMGC | See the paper and <https://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/skin_microbiome/README.txt> |
| VMGC | See the Zenodo record <https://doi.org/10.5281/zenodo.10457006> and <https://github.com/RChGO/VMGC> |
| Kraken2 | MIT (<https://github.com/DerrickWood/kraken2>) |

If you redistribute any derived database, check the upstream terms first.
