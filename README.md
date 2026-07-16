# MycobiomeDB

**MycobiomeDB: A Human Body Site-Specific Fungal Genome Database**

MycobiomeDB is a curated fungal genome resource designed for species-level profiling of the human mycobiome from shotgun metagenomic sequencing data. It provides body-site-specific fungal genome collections and prebuilt Kraken2 databases for four human body sites:

* Gut
* Oral cavity
* Skin
* Vagina

## Background

Fungi are an important but poorly characterized component of the human microbiome. Accurate profiling of the human mycobiome remains challenging because fungal reads are generally scarce in shotgun metagenomic data and currently available reference databases are often incomplete or insufficiently tailored to human-associated fungal communities.

To address these limitations, MycobiomeDB provides curated fungal reference panels stratified by human body site. These databases are intended to reduce false-positive taxonomic assignments while retaining fungal species relevant to each anatomical niche.

## Database Construction

We collected 16,365 publicly available high-quality fungal genomes representing 4,094 fungal species from the National Center for Biotechnology Information, NCBI, and the Joint Genome Institute, JGI.

Only genomes meeting the following quality criteria were included:

* Genome completeness greater than 90%
* Genome contamination lower than 5%

Shotgun metagenomic reads from 24,282 human samples were mapped against species-representative fungal genomes.

The analyzed metagenomes comprised:

| Body site   | Number of metagenomes |
| ----------- | --------------------: |
| Gut         |                 9,241 |
| Oral cavity |                 2,461 |
| Skin        |                 9,187 |
| Vagina      |                 3,393 |
| **Total**   |            **24,282** |

Species satisfying body-site-specific coverage and prevalence criteria were further curated based on nutritional mode and documented host association.

The final MycobiomeDB reference panels contain:

| Body site   | Number of fungal species |
| ----------- | -----------------------: |
| Gut         |                      766 |
| Oral cavity |                      259 |
| Skin        |                      642 |
| Vagina      |                      112 |

## Benchmarking

When used with Kraken2, MycobiomeDB generates species-level fungal taxonomic profiles with substantially fewer false-positive assignments than generic reference databases.

Benchmarking showed that fungal profiles generated using MycobiomeDB with Kraken2 were largely concordant with profiles generated using Bowtie2, while Kraken2 completed the analyses approximately 7- to 24-fold faster.

MycobiomeDB therefore provides an efficient alternative for fungal species profiling in large-scale human shotgun metagenomic studies.

## Repository Contents

This repository provides two types of resources for each body site:

1. **Fungal genome database**

   * Species-representative fungal genome sequences
   * Genome and taxonomic metadata
   * Files that can be used with sequence-alignment or custom database-building workflows

2. **Prebuilt Kraken2 database**

   * Ready-to-use Kraken2 database files
   * Constructed from the body-site-specific MycobiomeDB genome panel
   * Intended for direct species-level fungal profiling of shotgun metagenomic reads

The repository is organized as follows:

```text
MycobiomeDB/
├── Gut/
│   ├── Genomes/
│   ├── Kraken2_DB/
│   └── metadata/
├── Oral/
│   ├── Genomes/
│   ├── Kraken2_DB/
│   └── metadata/
├── Skin/
│   ├── Genomes/
│   ├── Kraken2_DB/
│   └── metadata/
├── Vagina/
│   ├── Genomes/
│   ├── Kraken2_DB/
│   └── metadata/
├── README.md
├── CITATION.cff
├── LICENSE
└── checksums.sha256
```

The exact directory structure may differ depending on the release package.

## Available Databases

### Gut MycobiomeDB

The gut-specific database contains 766 curated fungal species detected in or associated with the human gut.

Recommended for:

* Human fecal metagenomes
* Intestinal metagenomes
* Gut microbiome disease cohorts

### Oral MycobiomeDB

The oral-specific database contains 259 curated fungal species detected in or associated with the human oral cavity.

Recommended for:

* Saliva metagenomes
* Oral rinse metagenomes
* Dental plaque metagenomes
* Tongue or buccal samples

### Skin MycobiomeDB

The skin-specific database contains 642 curated fungal species detected in or associated with human skin.

Recommended for:

* Skin swab metagenomes
* Site-specific cutaneous microbiome studies
* Skin disease cohorts

### Vaginal MycobiomeDB

The vaginal-specific database contains 112 curated fungal species detected in or associated with the human vagina.

Recommended for:

* Vaginal swab metagenomes
* Female reproductive tract microbiome studies
* Vaginal disease cohorts

## Download

The current release can be downloaded from the ZENODO.

* Version-specific DOI: `10.5281/zenodo.21387880`
* Current version: `v1.0.0`

The database archive contains the body-site-specific genome collections and prebuilt Kraken2 databases for the gut, oral cavity, skin, and vagina.

## File Integrity

SHA-256 checksums are provided for the release files.

To verify downloaded files, run:

```bash
sha256sum -c checksums.sha256
```

All files should return `OK`.

## Kraken2 Usage

### Requirements

* Kraken2
* Sufficient storage space for the selected database
* Paired-end or single-end shotgun metagenomic sequencing reads

### Classification of paired-end reads

```bash
kraken2 \
  --db /path/to/MycobiomeDB/Gut/Kraken2_DB \
  --paired \
  --threads 16 \
  --report sample.kraken2.report \
  --output sample.kraken2.output \
  sample_R1.fastq.gz \
  sample_R2.fastq.gz
```

Replace `Gut` with `Oral`, `Skin`, or `Vagina` depending on the sample source.

### Classification of single-end reads

```bash
kraken2 \
  --db /path/to/MycobiomeDB/Gut/Kraken2_DB \
  --threads 16 \
  --report sample.kraken2.report \
  --output sample.kraken2.output \
  sample.fastq.gz
```

### Processing multiple paired-end samples

```bash
for r1 in *_R1.fastq.gz
do
    sample=${r1%_R1.fastq.gz}
    r2=${sample}_R2.fastq.gz

    kraken2 \
      --db /path/to/MycobiomeDB/Gut/Kraken2_DB \
      --paired \
      --threads 16 \
      --report ${sample}.kraken2.report \
      --output ${sample}.kraken2.output \
      ${r1} \
      ${r2}
done
```

## Choosing the Appropriate Database

Users should select the database corresponding to the anatomical origin of the metagenomic sample.

| Sample type                          | Recommended database |
| ------------------------------------ | -------------------- |
| Feces or intestinal sample           | Gut                  |
| Saliva, plaque, tongue, or oral swab | Oral                 |
| Skin swab or cutaneous sample        | Skin                 |
| Vaginal or reproductive tract sample | Vagina               |

Using a body-site-specific database is recommended because it limits the classification search space to fungal species supported by evidence from the corresponding human anatomical niche.

For samples whose anatomical origin does not match one of the four available sites, users should carefully evaluate whether one of the databases is biologically appropriate.

## Output Interpretation

The Kraken2 report contains taxonomic assignments at multiple ranks. Species-level results can be extracted from rows with the Kraken2 species rank code `S`.

Example:

```bash
awk -F '\t' '$4 == "S"' sample.kraken2.report \
  > sample.kraken2.species.report
```

Kraken2 read counts should be interpreted with caution because:

* Fungal DNA may be present at very low abundance.
* Closely related fungal species may share highly similar genomic regions.
* Contamination from laboratory reagents or the environment may affect low-biomass samples.
* Database detection does not by itself demonstrate fungal viability or active colonization.
* Comparisons across samples should account for differences in sequencing depth.

Appropriate filtering and normalization should be selected according to the study design.

## Companion Resource

Using the MycobiomeDB fungal-profiling pipeline, we also generated a companion resource containing body-site-stratified fungal species profiles for publicly available metagenomic disease cohorts.

The companion resource includes:

* Fungal species profiles from publicly available disease cohorts
* Body-site-specific abundance tables
* Comparisons between healthy and diseased samples
* Fungal taxa showing statistically significant differences between study groups

Information and access instructions for the companion resource will be provided separately.

## Versioning

MycobiomeDB uses semantic versioning:

* Major version changes indicate substantial changes to the database structure or construction procedure.
* Minor version changes indicate the addition of genomes, species, or metadata.
* Patch version changes indicate minor corrections that do not substantially alter database contents.

Users should record and report the exact MycobiomeDB version used in each analysis.

## Citation

MycobiomeDB has not yet been described in a peer-reviewed publication or preprint.

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

## License

MycobiomeDB is distributed under the `https://github.com/netbiolab/MycobiomeDB/blob/main/LICENSE` license.

The source genomes remain subject to the terms and policies of their original providers, including NCBI and JGI. Users are responsible for complying with any applicable terms associated with the original genomic resources.

## Data Sources

The database was constructed using publicly available fungal genomes obtained from:

* National Center for Biotechnology Information, NCBI
* Joint Genome Institute, JGI

Users requiring the original genome records should consult the corresponding source accession information provided in the metadata.

## Reporting Issues

Questions, database errors, missing taxa, and taxonomic inconsistencies can be reported through the GitHub Issues page.

When reporting an issue, please include:

* MycobiomeDB version
* Body-site-specific database used
* Relevant species or genome accession
* Description of the issue
* Example command or output, when applicable

## Contact

**Database maintainer:** `Sehun Ahn`
**Affiliation:** `Department of Biotechnology, College of Life Science and Biotechnology, Yonsei University, Seoul 03722, Republic of Korea`
**Email:** `zkskek321@yonsei.ac.kr`

## Disclaimer

MycobiomeDB is intended for research use only. Taxonomic classifications should be independently evaluated in the context of study design, sequencing depth, sample type, and potential contamination.

The database is not intended for direct clinical diagnosis or treatment decisions.
