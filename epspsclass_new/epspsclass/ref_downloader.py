"""
epspsclass/ref_downloader.py
============================
Fetch the four canonical EPSPS reference sequences used for classification.

Sources
-------
    vcEPSPS (class I):  UniProt Q9KNE7  (Vibrio cholerae O1 N16961)
    cbEPSPS (class II): UniProt Q83EH4  (Coxiella burnetii RSA 493)
    bvEPSPS (class III): NCBI   CAA73210 (Brevundimonas vesicularis)
    sdEPSPS (class IV): UniProt H6WNZ5  (Streptomyces davawensis JCM 4913)

These accessions are taken directly from Leino et al. (2021) supplementary
table 2.  The sequences are downloaded from public databases and cached in
the package data directory.  No API key is required.
"""

from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

_DATA_DIR = Path(__file__).parent / "data" / "reference_sequences"

_UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
_NCBI_FASTA_URL    = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=protein&id={accession}&rettype=fasta&retmode=text"
)

_REFERENCES = [
    # (class_label, filename, source, accession)
    ("I",   "vcEPSPS_Q9KNE7.fasta",   "uniprot", "Q9KNE7"),
    ("II",  "cbEPSPS_Q83EH4.fasta",   "uniprot", "Q83EH4"),
    ("III", "bvEPSPS_CAA73210.fasta", "ncbi",    "CAA73210"),
    ("IV",  "sdEPSPS_H6WNZ5.fasta",   "uniprot", "H6WNZ5"),
]


def _fetch_uniprot(accession: str) -> str:
    url = _UNIPROT_FASTA_URL.format(accession=accession)
    with urllib.request.urlopen(url, timeout=30) as resp:
        return resp.read().decode("utf-8")


def _fetch_ncbi(accession: str) -> str:
    url = _NCBI_FASTA_URL.format(accession=accession)
    with urllib.request.urlopen(url, timeout=30) as resp:
        return resp.read().decode("utf-8")


def download_all_references(force: bool = False) -> None:
    """
    Download all four reference FASTA files into the package data directory.

    Parameters
    ----------
    force : bool
        If True, re-download even if the file already exists.
    """
    _DATA_DIR.mkdir(parents=True, exist_ok=True)

    for cls, filename, source, accession in _REFERENCES:
        dest = _DATA_DIR / filename
        if dest.exists() and not force:
            print(f"  [skip] {filename} already present.", file=sys.stderr)
            continue

        print(f"  Downloading class {cls} reference ({accession}) …",
              file=sys.stderr, end=" ", flush=True)
        try:
            if source == "uniprot":
                data = _fetch_uniprot(accession)
            else:
                data = _fetch_ncbi(accession)
            dest.write_text(data)
            print("OK", file=sys.stderr)
        except Exception as exc:
            print(f"FAILED: {exc}", file=sys.stderr)
            raise RuntimeError(
                f"Could not download reference sequence for class {cls} "
                f"({accession}). Check your internet connection or download "
                f"manually from UniProt/NCBI and place the FASTA file at "
                f"{dest}."
            ) from exc
