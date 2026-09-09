"""Body-site definitions.

Each supported body site pins exactly one prokaryotic ("bacterial") decontamination
database and one MycobiomeDB fungal database. Selecting a site therefore selects
both databases; they are never mixed across sites.
"""

from collections import OrderedDict

# NOTE ON WORDING
# ---------------
# The decontamination databases below are *prokaryotic genome catalogues*. Each
# one contains bacteria plus a small number of archaea, and nothing else -- no
# human/host sequence, no fungi, no viruses. See `docs/DATABASES.md` for the
# per-catalogue composition and what that implies for interpretation.

SITES = OrderedDict(
    (
        (
            "gut",
            {
                "label": "Gut",
                "bacterial_catalogue": "HRGM2 (Human Reference Gut Microbiome v2), species representatives",
                "bacterial_dir": "gut",
                "fungal_catalogue": "Myco-HG (MycobiomeDB Gut)",
                "fungal_dir": "gut",
                "mycobiomedb_archive": "Gut.tar.gz",
            },
        ),
        (
            "oral",
            {
                "label": "Oral",
                "bacterial_catalogue": "HROM (Human Reference Oral Microbiome), species representatives",
                "bacterial_dir": "oral",
                "fungal_catalogue": "Myco-HO (MycobiomeDB Oral)",
                "fungal_dir": "oral",
                "mycobiomedb_archive": "Oral.tar.gz",
            },
        ),
        (
            "skin",
            {
                "label": "Skin",
                "bacterial_catalogue": "SMGC (Skin Microbial Genome Collection), prokaryotic representatives",
                "bacterial_dir": "skin",
                "fungal_catalogue": "Myco-HS (MycobiomeDB Skin)",
                "fungal_dir": "skin",
                "mycobiomedb_archive": "Skin.tar.gz",
            },
        ),
        (
            "vagina",
            {
                "label": "Vagina",
                "bacterial_catalogue": "VMGC (Vaginal Microbial Genome Collection), prokaryotic SGB representatives",
                "bacterial_dir": "vagina",
                "fungal_catalogue": "Myco-HV (MycobiomeDB Vaginal)",
                "fungal_dir": "vagina",
                "mycobiomedb_archive": "Vagina.tar.gz",
            },
        ),
    )
)

SITE_NAMES = list(SITES)


def get_site(name):
    """Return the site record for `name`, or raise KeyError with a usable message."""
    try:
        return SITES[name]
    except KeyError:
        raise KeyError(
            "unknown body site %r; supported sites are: %s" % (name, ", ".join(SITE_NAMES))
        )
