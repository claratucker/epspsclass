"""
epspsclass/groupby.py
=====================
Organism-level aggregation of per-sequence EPSPS classifications.

When a genome carries multiple EPSPS copies, classifying each sequence
independently and then aggregating by organism accession reveals:

  - Organisms with a single class (all copies agree)
  - Organisms with mixed classes (copies belong to different classes)
  - Organisms where some copies are unclassified

This is the organism-level "mixed" signal described in Leino et al. (2021)
Table 1 (intraspecific variation) and Rainio et al. (2021), which is distinct
from sequence-level mixed classification (a single sequence matching markers
from more than one class simultaneously).

Usage
-----
    epspsclass classify -i sequences.fasta -o results.tsv
    epspsclass group-by-organism \\
        --results results.tsv \\
        --metadata genome_metadata.tsv \\
        --output organism_summary.tsv

Metadata format (tab-separated, with header)
--------------------------------------------
    sequence_id    organism_id    [optional extra columns]

The sequence_id column must match the query_id values in the results TSV.
Any additional columns (e.g. taxonomy, habitat, exposure tier) are passed
through to the output.
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, TextIO


def _load_results(results_path: str) -> List[Dict[str, str]]:
    with open(results_path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _load_metadata(metadata_path: str) -> Dict[str, Dict[str, str]]:
    """Returns {sequence_id: {organism_id: ..., ...}}."""
    mapping = {}
    with open(metadata_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "sequence_id" not in (reader.fieldnames or []):
            raise ValueError(
                "Metadata file must have a 'sequence_id' column. "
                f"Found columns: {reader.fieldnames}"
            )
        if "organism_id" not in (reader.fieldnames or []):
            raise ValueError(
                "Metadata file must have an 'organism_id' column. "
                f"Found columns: {reader.fieldnames}"
            )
        for row in reader:
            mapping[row["sequence_id"]] = row
    return mapping


def group_by_organism(
    results: List[Dict[str, str]],
    metadata: Dict[str, Dict[str, str]],
) -> List[Dict[str, str]]:
    """
    Aggregate per-sequence classifications by organism.

    Returns one row per organism with the following columns:
      organism_id         - from metadata
      n_sequences         - total EPSPS sequences classified for this organism
      classes_found       - semicolon-separated set of all primary classes seen
      primary_class       - single class if all sequences agree; "Mixed" otherwise
      is_mixed            - True if sequences disagree on primary class
      n_unclassified      - count of unclassified sequences
      n_too_divergent     - count of sequences below identity threshold
      sensitive_count     - count of class I sequences
      resistant_count     - count of class II/III/IV sequences
      sequence_ids        - semicolon-separated list of contributing sequence IDs
      [extra metadata columns passed through from the metadata file]
    """
    # Group results by organism
    by_organism: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    unmapped: List[str] = []

    for row in results:
        seq_id = row["query_id"]
        if seq_id in metadata:
            org_id = metadata[seq_id]["organism_id"]
            by_organism[org_id].append(row)
        else:
            unmapped.append(seq_id)

    if unmapped:
        print(
            f"WARNING: {len(unmapped)} sequence(s) not found in metadata "
            f"and will be excluded from organism-level output.",
            file=sys.stderr,
        )

    # Get extra metadata columns to pass through
    extra_cols = []
    if metadata:
        first = next(iter(metadata.values()))
        extra_cols = [k for k in first.keys() if k not in ("sequence_id", "organism_id")]

    output_rows = []
    for org_id, rows in sorted(by_organism.items()):
        classes = [r["primary_class"] for r in rows]
        unique_classes = sorted(set(c for c in classes if c != "Unclassified"))

        n_unclassified = sum(1 for c in classes if c == "Unclassified")
        n_too_divergent = sum(
            1 for r in rows if r.get("is_too_divergent", "False") == "True"
        )
        sensitive_count = sum(1 for c in classes if c == "I")
        resistant_count = sum(1 for c in classes if c in ("II", "III", "IV"))

        is_mixed = len(unique_classes) > 1
        if is_mixed:
            primary = "Mixed"
        elif len(unique_classes) == 1:
            primary = unique_classes[0]
        else:
            primary = "Unclassified"

        # Pass through extra metadata from the first row for this organism
        meta_row = metadata.get(rows[0]["query_id"], {})
        extra = {k: meta_row.get(k, "") for k in extra_cols}

        output_rows.append({
            "organism_id":     org_id,
            "n_sequences":     str(len(rows)),
            "classes_found":   ";".join(unique_classes) if unique_classes else "Unclassified",
            "primary_class":   primary,
            "is_mixed":        str(is_mixed),
            "n_unclassified":  str(n_unclassified),
            "n_too_divergent": str(n_too_divergent),
            "sensitive_count": str(sensitive_count),
            "resistant_count": str(resistant_count),
            "sequence_ids":    ";".join(r["query_id"] for r in rows),
            **extra,
        })

    return output_rows


ORGANISM_FIELDNAMES = [
    "organism_id",
    "n_sequences",
    "classes_found",
    "primary_class",
    "is_mixed",
    "n_unclassified",
    "n_too_divergent",
    "sensitive_count",
    "resistant_count",
    "sequence_ids",
]


def write_organism_summary(
    rows: List[Dict[str, str]],
    output: TextIO,
    delimiter: str = "\t",
) -> None:
    if not rows:
        return
    extra_cols = [k for k in rows[0].keys() if k not in ORGANISM_FIELDNAMES]
    fieldnames = ORGANISM_FIELDNAMES + extra_cols
    writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter=delimiter)
    writer.writeheader()
    writer.writerows(rows)
