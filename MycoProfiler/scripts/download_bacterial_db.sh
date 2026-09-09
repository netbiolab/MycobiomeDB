#!/usr/bin/env bash
#
# Download the body-site prokaryotic ("bacterial") decontamination databases used by
# MycoProfiler stage 1, and install them into the MycoProfiler database layout.
#
#   scripts/download_bacterial_db.sh --db-root /data/mycoprofiler_db [--site gut] ...
#
# THREE of the four sites publish a ready-to-use Kraken2 database; SKIN does not.
#
#   site    catalogue  distributed as                              action here
#   ------  ---------  ------------------------------------------  ----------------------
#   gut     HRGM2      prebuilt Kraken2 DB (hash/opts/taxo .k2d)    downloaded, ready
#   oral    HROM       prebuilt Kraken2 DB (hash/opts/taxo .k2d)    downloaded, ready
#   vagina  VMGC       prebuilt Kraken2 DB (.tar.gz)                downloaded, ready
#   skin    SMGC       genome FASTA catalogue only                  sources downloaded,
#                                                                   then run
#                                                                   scripts/build_smgc_kraken2_db.py
#
# Every URL and byte size below was confirmed with an HTTP HEAD request against the
# provider at the time of writing. See docs/DATABASES.md for the citations.
set -euo pipefail

DECODE="https://www.decodebiome.org"
HRGM2_CONCAT="${DECODE}/HRGM2/data/genome_catalog/Taxonomy_Profiling/HRGMv2_kraken2_customdb/HRGMv2_Concat"
HRGM2_REP="${DECODE}/HRGM2/data/genome_catalog/Taxonomy_Profiling/HRGMv2_kraken2_customdb/HRGMv2_Rep"
HROM_DB="${DECODE}/HROM/data/genome_catalog/HROM_kraken2_customdb"
VMGC_TAR="https://zenodo.org/api/records/10457006/files/VMGC_prokaryote_SGB_KrakenDB.tar.gz/content"
SMGC_FTP="https://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/skin_microbiome"

# Confirmed Content-Length values (bytes).
HRGM2_CONCAT_HASH_BYTES=37239042956
HRGM2_REP_HASH_BYTES=17791303712
HROM_HASH_BYTES=12503236052
VMGC_TAR_BYTES=2199401770

DB_ROOT=""
SITES=()
HRGM2_VARIANT="concat"
KEEP_ARCHIVE=0

usage() {
  cat <<'USAGE'
Usage: download_bacterial_db.sh --db-root DIR [options]

  --db-root DIR        install into DIR/bacterial/<site>/                (required)
  --site NAME          gut|oral|skin|vagina; repeatable (default: all four)
  --hrgm2-variant V    concat (default) | rep -- which HRGM2 Kraken2 build to use
  --keep-archive       keep downloaded archives after extraction
  -h, --help           this message

HRGM2 variants: both are built from the same 4,824 species-representative genomes.
  concat  contigs concatenated per genome; hash.k2d 34.7 GiB. This is the build the
          published MycoProfiler analyses used -- keep it unless you have a reason not to.
  rep     contigs kept separate;           hash.k2d 16.6 GiB.

Approximate download sizes: gut 34.7 GiB (concat) / 16.6 GiB (rep), oral 11.6 GiB,
vagina 2.0 GiB, skin 0.6 GiB of source genomes (the DB is then built locally).
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --db-root) DB_ROOT="${2:?--db-root needs a value}"; shift 2 ;;
    --site)    SITES+=("${2:?--site needs a value}"); shift 2 ;;
    --hrgm2-variant) HRGM2_VARIANT="${2:?--hrgm2-variant needs a value}"; shift 2 ;;
    --keep-archive) KEEP_ARCHIVE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$DB_ROOT" ]] || { echo "error: --db-root is required" >&2; usage >&2; exit 2; }
[[ ${#SITES[@]} -gt 0 ]] || SITES=(gut oral skin vagina)
command -v curl >/dev/null || { echo "error: curl is required" >&2; exit 1; }

file_bytes() { stat -c %s "$1" 2>/dev/null || stat -f %z "$1"; }

fetch() {  # completed files have a local checksum receipt; partial bytes stay resumable
  local url="$1" dest="$2" expect="${3:-}" part="${2}.part"
  mkdir -p "$(dirname "$dest")"
  if [[ -s "$dest" && -f "${dest}.sha256" ]] &&
     sha256sum -c "${dest}.sha256" >/dev/null 2>&1; then
    if [[ -z "$expect" || "$(file_bytes "$dest")" == "$expect" ]]; then
      echo "    verified existing $(basename "$dest")"
      return 0
    fi
  fi
  # Recover partial known-size files left by the previous downloader.
  if [[ -f "$dest" ]]; then
    if [[ -n "$expect" && ! -e "$part" && "$(file_bytes "$dest")" -lt "$expect" ]]; then
      mv "$dest" "$part"
    else
      mv "$dest" "${dest}.unverified.$(date +%s%N)"
    fi
  fi
  rm -f "${dest}.sha256"
  echo "    GET $(basename "$dest")"
  curl -fL -C - --retry 5 --retry-delay 10 -o "$part" "$url"
  [[ -s "$part" ]] || { echo "ERROR: empty download: $url" >&2; return 1; }
  if [[ -n "$expect" && "$(file_bytes "$part")" != "$expect" ]]; then
    echo "ERROR: unexpected size for $part; partial download retained" >&2
    return 1
  fi
  mv "$part" "$dest"
  sha256sum "$dest" > "${dest}.sha256"
}

validate_k2d() {
  local target="$1" name report
  for name in hash.k2d opts.k2d taxo.k2d; do
    [[ -s "$target/$name" ]] || { echo "ERROR: missing/empty $target/$name" >&2; return 1; }
  done
  report=$(mktemp)
  if ! kraken2-inspect --db "$target" > "$report"; then
    rm -f "$report"; echo "ERROR: Kraken2 DB inspection failed" >&2; return 1
  fi
  if ! awk -F '\t' '$4 == "R" && $5 == "1" && $2+0 > 0 {ok=1} END {exit !ok}' "$report"; then
    rm -f "$report"; echo "ERROR: no positive root minimizer count" >&2; return 1
  fi
  rm -f "$report"
}

install_flat_k2d() {
  local base="$1" target="$2" hash_bytes="$3"
  mkdir -p "$target"
  # A changed URL (e.g. concat -> rep) must not reuse another variant's bytes.
  if [[ -f "$target/.mycoprofiler-source" && "$(cat "$target/.mycoprofiler-source")" != "$base" ]]; then
    mv "$target" "${target}.backup.$(date +%s%N)"
    mkdir -p "$target"
  fi
  printf '%s\n' "$base" > "$target/.mycoprofiler-source"
  rm -f "$target/.mycoprofiler-complete"
  fetch "${base}/opts.k2d" "${target}/opts.k2d" 64
  fetch "${base}/taxo.k2d" "${target}/taxo.k2d"
  fetch "${base}/hash.k2d" "${target}/hash.k2d" "$hash_bytes"
  validate_k2d "$target"
  touch "$target/.mycoprofiler-complete"
}

for tool in sha256sum kraken2-inspect; do
  command -v "$tool" >/dev/null || { echo "ERROR: $tool is required" >&2; exit 1; }
done
case "$HRGM2_VARIANT" in concat|rep) ;; *) echo "ERROR: invalid HRGM2 variant" >&2; exit 2 ;; esac
for site in "${SITES[@]}"; do
  case "$site" in gut|oral|skin|vagina) ;; *) echo "ERROR: unknown site: $site" >&2; exit 2 ;; esac
done

for site in "${SITES[@]}"; do
  target="${DB_ROOT}/bacterial/${site}"

  case "$site" in
    gut)
      case "$HRGM2_VARIANT" in
        concat) base="$HRGM2_CONCAT"; bytes="$HRGM2_CONCAT_HASH_BYTES" ;;
        rep)    base="$HRGM2_REP";    bytes="$HRGM2_REP_HASH_BYTES" ;;
        *) echo "error: --hrgm2-variant must be 'concat' or 'rep'" >&2; exit 2 ;;
      esac
      echo "[gut] HRGM2 (${HRGM2_VARIANT}) prebuilt Kraken2 database"
      install_flat_k2d "$base" "$target" "$bytes"
      echo "[gut] installed -> ${target}"
      ;;

    oral)
      echo "[oral] HROM prebuilt Kraken2 database"
      install_flat_k2d "$HROM_DB" "$target" "$HROM_HASH_BYTES"
      echo "[oral] installed -> ${target}"
      ;;

    vagina)
      echo "[vagina] VMGC prokaryotic SGB prebuilt Kraken2 database"
      if [[ -f "${target}/.mycoprofiler-complete" && -f "${target}/.mycoprofiler-db.sha256" ]] &&
         sha256sum -c "${target}/.mycoprofiler-db.sha256" >/dev/null 2>&1 && validate_k2d "$target"; then
        echo "[vagina] verified existing database -- skipping"
        continue
      fi
      work="${DB_ROOT}/.download/bacterial_vagina"
      archive="${work}/VMGC_prokaryote_SGB_KrakenDB.tar.gz"
      fetch "$VMGC_TAR" "$archive" "$VMGC_TAR_BYTES"
      echo "    extracting"
      extract="${work}/extracted"; rm -rf "$extract"; mkdir -p "$extract"
      tar -xzf "$archive" -C "$extract"
      hash_path="$(find "$extract" -maxdepth 3 -name hash.k2d -print -quit)"
      [[ -n "$hash_path" ]] || { echo "[vagina] ERROR: no hash.k2d in the archive" >&2; exit 1; }
      validate_k2d "$(dirname "$hash_path")"
      mkdir -p "$(dirname "$target")"
      if [[ -e "$target" ]]; then mv "$target" "${target}.backup.$(date +%s%N)"; fi
      mv "$(dirname "$hash_path")" "$target"
      sha256sum "$target/hash.k2d" "$target/opts.k2d" "$target/taxo.k2d" > "$target/.mycoprofiler-db.sha256"
      touch "$target/.mycoprofiler-complete"
      [[ $KEEP_ARCHIVE -eq 0 ]] && rm -rf "$work"
      echo "[vagina] installed -> ${target}"
      ;;

    skin)
      echo "[skin] SMGC publishes genome FASTA only -- no prebuilt Kraken2 database."
      src="${DB_ROOT}/sources/smgc"
      fetch "${SMGC_FTP}/03_finalcatalogue.tar.gz" "${src}/03_finalcatalogue.tar.gz"
      fetch "${SMGC_FTP}/SMGC.xlsx"                "${src}/SMGC.xlsx"
      fetch "${SMGC_FTP}/README.txt"               "${src}/README.txt"
      echo "    extracting genome catalogue"
      tar -xzf "${src}/03_finalcatalogue.tar.gz" -C "${src}"
      cat <<EOSKIN

[skin] Sources are in ${src}. Build the Kraken2 database with:

    python scripts/build_smgc_kraken2_db.py \\
        --catalogue ${src}/03_finalcatalogue \\
        --metadata  ${src}/SMGC.xlsx \\
        --output    ${target} \\
        --threads 16

That step needs pandas + openpyxl (see requirements.txt) and kraken2-build on \$PATH.
It selects only the 622 PROKARYOTIC SMGC genomes; the catalogue's 12 eukaryotic and
6,935 viral genomes are deliberately excluded from the decontamination database.
EOSKIN
      ;;

    *) echo "error: unknown site '$site'" >&2; exit 2 ;;
  esac
done

echo
echo "Done. Verify with:"
echo "  mycoprofiler check --db-root ${DB_ROOT}"
