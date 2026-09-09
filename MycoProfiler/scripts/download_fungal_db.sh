#!/usr/bin/env bash
#
# Download the MycobiomeDB body-site fungal Kraken2 databases (Myco-HG/HO/HS/HV)
# from Zenodo and install them into the MycoProfiler database layout.
#
#   scripts/download_fungal_db.sh --db-root /data/mycoprofiler_db [--site gut] [--site skin]
#
# Everything below (record id, file names, byte sizes, MD5 sums) was read from the
# Zenodo REST API for record 21387880 and from that record's own checksums.md5.
#
# Source : MycobiomeDB v1.0.0, DOI 10.5281/zenodo.21387880 (CC BY 4.0)
#          Sehun Ahn and Insuk Lee
# Note   : the MycobiomeDB GitHub README refers to a `checksums.sha256`; the Zenodo
#          record actually publishes `checksums.md5`. MD5 is what is verified here.
set -euo pipefail

RECORD=21387880
BASE="https://zenodo.org/api/records/${RECORD}/files"

# site : archive : bytes : md5
SITE_GUT="Gut.tar.gz:40698406476:401e5a4322c3bf25f14be0727e5c7f9f"
SITE_ORAL="Oral.tar.gz:11398360153:3284577302089ce0786c88216a10f121"
SITE_SKIN="Skin.tar.gz:27052137193:08f327a8cace8702c2beee969f1ceea6"
SITE_VAGINA="Vagina.tar.gz:5438786935:e7808ba3e3129d16752b2e1f5d35eb3a"

DB_ROOT=""
SITES=()
KEEP_ARCHIVE=0
SKIP_VERIFY=0

usage() {
  cat <<'USAGE'
Usage: download_fungal_db.sh --db-root DIR [options]

  --db-root DIR     install into DIR/fungal/<site>/           (required)
  --site NAME       gut|oral|skin|vagina; repeatable (default: all four)
  --keep-archive    keep the downloaded .tar.gz after extraction
  --skip-verify     skip the MD5 check (not recommended)
  -h, --help        this message

Download sizes (compressed, exact, from Zenodo):
  Gut 37.9 GiB | Oral 10.6 GiB | Skin 25.2 GiB | Vagina 5.1 GiB

Disk: while extracting, both the archive and its contents exist, so allow roughly
2.2x the compressed size per site unless you use --keep-archive=off behaviour
(the archive is deleted right after a successful extraction by default).
Resuming: downloads use `curl -C -`, so re-running continues a partial file.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --db-root) DB_ROOT="${2:?--db-root needs a value}"; shift 2 ;;
    --site)    SITES+=("${2:?--site needs a value}"); shift 2 ;;
    --keep-archive) KEEP_ARCHIVE=1; shift ;;
    --skip-verify)  SKIP_VERIFY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$DB_ROOT" ]] || { echo "error: --db-root is required" >&2; usage >&2; exit 2; }
[[ ${#SITES[@]} -gt 0 ]] || SITES=(gut oral skin vagina)

command -v curl >/dev/null || { echo "error: curl is required" >&2; exit 1; }
command -v tar  >/dev/null || { echo "error: tar is required" >&2; exit 1; }

md5_of() {
  if command -v md5sum >/dev/null; then md5sum "$1" | awk '{print $1}'
  elif command -v md5 >/dev/null; then md5 -q "$1"
  else echo ""; fi
}

for site in "${SITES[@]}"; do
  case "$site" in
    gut)    spec="$SITE_GUT" ;;
    oral)   spec="$SITE_ORAL" ;;
    skin)   spec="$SITE_SKIN" ;;
    vagina) spec="$SITE_VAGINA" ;;
    *) echo "error: unknown site '$site' (use gut|oral|skin|vagina)" >&2; exit 2 ;;
  esac
  archive="${spec%%:*}"; rest="${spec#*:}"
  bytes="${rest%%:*}"; md5="${rest##*:}"

  target="${DB_ROOT}/fungal/${site}"
  if [[ -f "${target}/hash.k2d" ]]; then
    echo "[${site}] already installed at ${target} -- skipping"
    continue
  fi

  work="${DB_ROOT}/.download/fungal_${site}"
  mkdir -p "$work"
  archive_path="${work}/${archive}"

  echo "[${site}] downloading ${archive} (${bytes} bytes) from Zenodo record ${RECORD}"
  curl -fL -C - --retry 5 --retry-delay 10 -o "$archive_path" "${BASE}/${archive}/content"

  actual_bytes=$(stat -c %s "$archive_path" 2>/dev/null || stat -f %z "$archive_path")
  if [[ "$actual_bytes" != "$bytes" ]]; then
    echo "[${site}] ERROR: size mismatch (expected ${bytes}, got ${actual_bytes})" >&2
    echo "         delete ${archive_path} and retry" >&2
    exit 1
  fi

  if [[ $SKIP_VERIFY -eq 0 ]]; then
    echo "[${site}] verifying MD5 (this reads the whole archive and takes a few minutes)"
    actual_md5=$(md5_of "$archive_path")
    if [[ -z "$actual_md5" ]]; then
      echo "[${site}] WARNING: no md5sum/md5 tool found; checksum NOT verified" >&2
    elif [[ "$actual_md5" != "$md5" ]]; then
      echo "[${site}] ERROR: MD5 mismatch (expected ${md5}, got ${actual_md5})" >&2
      exit 1
    else
      echo "[${site}] MD5 OK"
    fi
  fi

  echo "[${site}] extracting"
  extract_dir="${work}/extracted"
  rm -rf "$extract_dir"; mkdir -p "$extract_dir"
  tar -xzf "$archive_path" -C "$extract_dir"

  # The published layout is <Site>/Kraken2_DB/, but the README warns it "may differ
  # depending on the release package", so locate the Kraken2 database by its files.
  hash_path="$(find "$extract_dir" -maxdepth 4 -name hash.k2d -print -quit)"
  if [[ -z "$hash_path" ]]; then
    echo "[${site}] ERROR: no hash.k2d found under ${extract_dir}." >&2
    echo "         The archive layout differs from the expected one. Inspect it and" >&2
    echo "         move the directory holding hash.k2d/opts.k2d/taxo.k2d to ${target}." >&2
    exit 1
  fi
  db_dir="$(dirname "$hash_path")"
  for required in opts.k2d taxo.k2d; do
    [[ -f "${db_dir}/${required}" ]] || {
      echo "[${site}] ERROR: ${db_dir} has hash.k2d but no ${required}" >&2; exit 1; }
  done

  mkdir -p "$(dirname "$target")"
  rm -rf "$target"
  mv "$db_dir" "$target"
  echo "[${site}] installed -> ${target}"

  # Keep any sibling genome/metadata folders next to the DB for provenance.
  siblings="${DB_ROOT}/fungal/${site}_source"
  if [[ -n "$(ls -A "$extract_dir" 2>/dev/null)" ]]; then
    rm -rf "$siblings"; mkdir -p "$siblings"
    mv "$extract_dir"/* "$siblings"/ 2>/dev/null || true
    echo "[${site}] genome/metadata files from the archive kept in ${siblings}"
  fi

  if [[ $KEEP_ARCHIVE -eq 0 ]]; then rm -f "$archive_path"; fi
  rmdir "$extract_dir" 2>/dev/null || true
done

echo
echo "Done. Verify with:"
echo "  mycoprofiler check --db-root ${DB_ROOT}"
