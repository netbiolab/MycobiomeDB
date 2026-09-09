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

## MycoProfiler: Body-site-specific fungal profiling

MycoProfiler processes human shotgun metagenomic reads in two steps:
1. Kraken2-based bacterial decontamination using a body-site-specific prokaryotic database.
2. Kraken2 classification of the remaining unclassified reads using the corresponding MycobiomeDB fungal database.

MycoProfiler supports gut, oral, skin, and vaginal samples, with a Kraken2 confidence threshold of 0.2 in both steps.

**[MycoProfiler installation and usage guide](MycoProfiler/README.md)**

### Downstream abundance estimation with Bracken

After running MycoProfiler, users can estimate species-level fungal abundances with Bracken using the stage-2 fungal Kraken2 report. The Bracken read-length parameter must match the read length used to generate the corresponding distribution file for the same body-site-specific fungal database.

Bracken is run separately and is not automatically executed by MycoProfiler.

See [Bracken abundance estimation](MycoProfiler/README.md#bracken-abundance-estimation) for instructions.

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

Bracken estimates species-level fungal abundances by redistributing Kraken2 assignments across taxonomic ranks using database-specific, read-length-dependent probabilities. When Bracken is run with `-l S`, its output provides species-level estimates.

Key columns in the Bracken output include:

| Column                  | Description                                                        |
| ----------------------- | ------------------------------------------------------------------ |
| `name`                  | Taxon name                                                         |
| `taxonomy_id`           | Taxonomy identifier used in the database                           |
| `taxonomy_lvl`          | Taxonomic rank (`S` for species)                                   |
| `kraken_assigned_reads` | Count assigned to the taxon before Bracken re-estimation           |
| `added_reads`           | Additional count allocated by Bracken                              |
| `new_est_reads`         | Estimated count after re-estimation                                |
| `fraction_total_reads`  | Relative abundance within the total abundance estimated by Bracken |

Use `new_est_reads` for estimated count tables and `fraction_total_reads` for relative-abundance profiles. For paired-end Kraken2 classification, the input counts represent read pairs; interpret the resulting estimates in the same counting unit.

Because Bracken is applied to the fungal classification report, these relative abundances describe the estimated fungal community represented in the output. They do not measure the fungal fraction of the original metagenome or absolute fungal abundance.

Bracken estimates should be interpreted with caution because:

* Fungal DNA may be present at very low abundance.
* Shared genomic sequences and incomplete references can affect abundance estimates.
* The Bracken distribution file must match the fungal database version and selected read length.
* Low-abundance estimates may be affected by sequencing depth, the Bracken threshold, and contamination.
* Bracken does not independently confirm fungal detection, viability, or active colonization.

Apply study-appropriate filtering and normalization, and use consistent database versions and documented parameter-selection rules across samples.

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

This project is licensed under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/deed.en).

[![CC BY 4.0][cc-by-shield]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

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
