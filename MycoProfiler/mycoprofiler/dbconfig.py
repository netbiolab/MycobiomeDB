"""Resolution and validation of Kraken2 database locations.

A Kraken2 database is a *directory* that must contain at least these three files:

    hash.k2d   opts.k2d   taxo.k2d

Everything else a build leaves behind (``library/``, ``taxonomy/``,
``database.kraken``, ``database*mers.kmer_distrib``, ``seqid2taxid.map``) is not
needed to classify reads and may be deleted to save disk. MycoProfiler only checks
for, and only reads, the three files above.

Resolution order for a given site (first hit wins, highest priority first):

    1. ``--bacterial-db`` / ``--fungal-db``      explicit paths on the command line
    2. ``--db-config FILE``                      a JSON file mapping site -> paths
    3. ``$MYCOPROFILER_DB_CONFIG``                    same JSON file, from the environment
    4. ``--db-root DIR``                         conventional layout (see below)
    5. ``$MYCOPROFILER_DB_ROOT``                      same conventional layout

Conventional layout under a database root::

    <db-root>/bacterial/{gut,oral,skin,vagina}/
    <db-root>/fungal/{gut,oral,skin,vagina}/

JSON is used rather than YAML so that MycoProfiler needs no third-party runtime
dependency. Top-level keys beginning with an underscore are ignored, which gives
JSON the comment field it otherwise lacks. A config file is a plain object, e.g.::

    {
      "gut":  {"bacterial": "/data/db/HRGM2",  "fungal": "/data/db/Myco-HG"},
      "skin": {"bacterial": "/data/db/SMGC",   "fungal": "/data/db/Myco-HS"}
    }
"""

import json
import os

from .sites import get_site

#: Files Kraken2 actually needs at classification time.
REQUIRED_DB_FILES = ("hash.k2d", "opts.k2d", "taxo.k2d")

ENV_DB_ROOT = "MYCOPROFILER_DB_ROOT"
ENV_DB_CONFIG = "MYCOPROFILER_DB_CONFIG"


class ConfigError(Exception):
    """Raised when database locations cannot be resolved or are unusable."""


class ResolvedDb(object):
    """A resolved database path plus a note on where the path came from."""

    __slots__ = ("path", "source")

    def __init__(self, path, source):
        self.path = os.path.abspath(path)
        self.source = source

    def __repr__(self):
        return "ResolvedDb(%r, from=%r)" % (self.path, self.source)


def load_db_config(path):
    """Read and sanity-check a MycoProfiler JSON database config."""
    if not os.path.isfile(path):
        raise ConfigError("database config file not found: %s" % path)
    try:
        with open(path) as handle:
            data = json.load(handle)
    except ValueError as exc:
        raise ConfigError("database config %s is not valid JSON: %s" % (path, exc))
    if not isinstance(data, dict):
        raise ConfigError("database config %s must contain a JSON object at the top level" % path)
    for site, entry in data.items():
        if site.startswith("_"):
            # JSON has no comment syntax; keys starting with '_' are treated as notes.
            continue
        if not isinstance(entry, dict):
            raise ConfigError(
                "database config %s: entry for site %r must be an object with "
                "'bacterial' and/or 'fungal' keys" % (path, site)
            )
        unknown = set(entry) - {"bacterial", "fungal"}
        if unknown:
            raise ConfigError(
                "database config %s: site %r has unsupported key(s) %s "
                "(expected 'bacterial' and/or 'fungal')"
                % (path, site, ", ".join(sorted(unknown)))
            )
    return data


def resolve_databases(site, bacterial_db=None, fungal_db=None, db_config=None, db_root=None,
                      environ=None):
    """Resolve the (bacterial, fungal) database directories for `site`.

    Returns a dict with keys ``bacterial`` and ``fungal`` mapping to ResolvedDb.
    Raises ConfigError with an actionable message if either cannot be resolved.
    Existence of the directories is *not* checked here; see `validate_db`.
    """
    environ = os.environ if environ is None else environ
    record = get_site(site)

    config_path = db_config or environ.get(ENV_DB_CONFIG)
    config = load_db_config(config_path) if config_path else {}
    config_entry = config.get(site, {}) if isinstance(config, dict) else {}

    root = db_root or environ.get(ENV_DB_ROOT)

    resolved = {}
    for kind, explicit in (("bacterial", bacterial_db), ("fungal", fungal_db)):
        if explicit:
            resolved[kind] = ResolvedDb(explicit, "--%s-db" % kind)
            continue
        if config_entry.get(kind):
            resolved[kind] = ResolvedDb(
                config_entry[kind], "db config %s (site '%s')" % (config_path, site)
            )
            continue
        if root:
            subdir = record["%s_dir" % kind]
            resolved[kind] = ResolvedDb(
                os.path.join(root, kind, subdir),
                "db root %s (%s/%s)" % (root, kind, subdir),
            )
            continue
        raise ConfigError(
            "cannot locate the %s database for site '%s'.\n"
            "Provide one of:\n"
            "  --%s-db /path/to/db            (explicit path)\n"
            "  --db-config databases.json     (or set $%s)\n"
            "  --db-root /path/to/databases   (or set $%s), which expects "
            "<db-root>/%s/%s/"
            % (kind, site, kind, ENV_DB_CONFIG, ENV_DB_ROOT, kind, record["%s_dir" % kind])
        )
    return resolved


def missing_db_files(path):
    """Return the list of REQUIRED_DB_FILES absent from `path`."""
    return [name for name in REQUIRED_DB_FILES if not os.path.isfile(os.path.join(path, name))]


def validate_db(resolved, kind, site):
    """Check that a resolved database directory is a usable Kraken2 database."""
    path = resolved.path
    if not os.path.exists(path):
        raise ConfigError(
            "%s database for site '%s' does not exist: %s\n"
            "(resolved from %s)\n"
            "Download or build it first -- see scripts/download_%s_db.sh and README.md."
            % (kind, site, path, resolved.source, kind)
        )
    if not os.path.isdir(path):
        raise ConfigError(
            "%s database for site '%s' is not a directory: %s\n"
            "A Kraken2 database is the directory holding hash.k2d/opts.k2d/taxo.k2d, "
            "not a single file." % (kind, site, path)
        )
    missing = missing_db_files(path)
    if missing:
        raise ConfigError(
            "%s database for site '%s' at %s is missing required Kraken2 file(s): %s\n"
            "(resolved from %s)\n"
            "If the download or build was interrupted, re-run it; if the files are one "
            "level deeper, point at that subdirectory instead."
            % (kind, site, path, ", ".join(missing), resolved.source)
        )
    return True


def db_fingerprint(path):
    """Small, cheap, reproducibility-oriented description of a database directory.

    Records size and mtime of the three required files rather than hashing the
    (tens of GB) hash table, which would dominate runtime.
    """
    entry = {"path": os.path.abspath(path), "files": {}}
    for name in REQUIRED_DB_FILES:
        full = os.path.join(path, name)
        try:
            stat = os.stat(full)
        except OSError:
            entry["files"][name] = None
            continue
        entry["files"][name] = {"bytes": stat.st_size, "mtime_epoch": int(stat.st_mtime)}
    return entry
