# Corrections — 2026-09-09

## Replacing files

Replace `scripts/download_bacterial_db.sh`, `scripts/build_smgc_kraken2_db.py`,
and `mycoprofiler/pipeline.py` with the supplied files, preserving these exact names.
The complete ZIP additionally includes updated tests and documentation.

## Behavior

- The default SMGC `--no-masking` flag is now passed to `--add-to-library`.
  `--mask-low-complexity` omits that flag and requires the appropriate masking tool.
- Bacterial downloads use `.part` files and resume interrupted transfers. Legacy
  known-size partial files are moved into `.part` for resumption. Unverified
  existing complete files are retained with an `.unverified.*` suffix and fetched
  again. Local SHA-256 receipts detect subsequent byte changes; these are NOT
  provider-published checksums and do not authenticate the source. Known provider
  byte counts are checked. `kraken2-inspect` must succeed and report a positive root
  minimizer count before DB installation is marked complete. Reuse reads the full
  files to check SHA-256, which can take time for large DBs. An HRGM2 URL variant
  change backs up the old directory. Downloads and builds to the same destination
  must not run concurrently.
- VMGC extraction is inspected before installation; existing targets are backed up.
- Pipeline `--force` moves recognized old outputs to an adjacent
  `.backup-<timestamp>` directory, then starts with a fresh output directory.
  Unrecognized files and symlink output directories are refused. Failed reruns
  therefore cannot expose an old success summary/manifest in the current output.
  Backups are retained; account for their disk usage.
- SMGC `--force` likewise backs up recognized previous build files and starts fresh.
  Missing genomes and failed/empty/unparseable inspections stop the build.
  `build_provenance.json` retains genome-to-taxid mappings, taxonomy names and commands.
- Sample IDs cannot contain path separators or tabs/newlines.
- Manifest `handoff.mates_synchronised` was renamed to
  `handoff.mate_record_counts_equal` to describe the actual count-only check.
  Read identifiers/order are not independently verified.

## Validation

49 tests passed with mock external programs; Python compilation and Bash syntax
checks passed. No multi-GB download, actual Kraken2 build or real classification
was performed for these corrections. Run from the MycoProfiler directory:

```bash
chmod +x tests/mock_kraken2.py
python -m unittest discover -s tests -v
```

The fungal downloader and the shallow `mycoprofiler check` DB validation were not
changed in this patch. Their earlier limitations remain. Original DB URLs,
provider byte sizes, and HRGM2 variant assumptions were not reverified.

## Name change to MycoProfiler

The outer project directory, Python package, command entry point, Conda environment,
environment variables, output filenames, tests and documentation use MycoProfiler /
mycoprofiler / MYCOPROFILER. MycobiomeDB and Myco-HG/HO/HS/HV are unchanged.
The old mycomap command and MYCOMAP_* variables are not aliases.
49 mock-based tests passed after renaming; local package installation and
`mycoprofiler --version` also passed. No real DB classification was performed.

Upload the MycoProfiler directory to the repository root. Update the repository's
root README link to `MycoProfiler/README.md`, and remove the old MycoMap directory
on the same working branch. The root README is not included in this package.
