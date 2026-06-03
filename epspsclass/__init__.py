"""
EPSPSClass
==========
Open-source reimplementation of the EPSPS glyphosate-sensitivity classifier.

Based on the algorithm described in:
    Leino et al. (2021) Classification of the glyphosate target enzyme
    (5-enolpyruvylshikimate-3-phosphate synthase) for assessing sensitivity
    of organisms to the herbicide. Environment International 149:106334.
    https://doi.org/10.1016/j.envint.2020.106334

The original web server (http://ppuigbo.me/programs/EPSPSClass) is no longer
reliably available.  This package provides a fully reproducible, locally
runnable, open-source alternative with identical classification logic.

Differences from original tool
-------------------------------
- Fully open source (MIT licence) with all logic exposed
- Runs locally; no internet required (reference sequences are bundled)
- Batch FASTA input with TSV output; AWS-deployment ready
- Marker table is versioned and human-readable in classifier.py
- Validates reference sequence integrity on load
- CLI plus importable Python API

Citation
--------
If you use EPSPSClass in published research, please cite both:
    1. Leino et al. (2021): the original classification framework
    2. [Your paper]: for this reimplementation

Quickstart
----------
    # Classify sequences (reference sequences are bundled (no setup needed))
    epspsclass classify -i my_sequences.fasta -o results.tsv

    # Python API
    from epspsclass import EPSPSClassifier
    clf = EPSPSClassifier()
    result = clf.classify("seq1", seq_object)
"""

from .classifier import EPSPSClassifier, ClassificationResult
from .ref_downloader import download_all_references

__version__ = "1.0.3"
__all__ = ["EPSPSClassifier", "ClassificationResult", "download_all_references"]
