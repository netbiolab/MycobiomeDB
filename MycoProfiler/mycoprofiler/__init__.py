"""MycoProfiler: body-site-specific bacterial decontamination followed by fungal
classification of human shotgun metagenomes with Kraken2."""

__version__ = "0.1.0"

# Both Kraken2 stages of MycoProfiler run at this confidence threshold. It is a
# constant, not a tunable, so that every MycoProfiler result is comparable.
CONFIDENCE = 0.2

__all__ = ["__version__", "CONFIDENCE"]
