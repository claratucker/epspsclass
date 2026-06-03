"""
epspsclass/cli.py
=================
Command-line interface for EPSPSClass.

Pipeline integration
--------------------
EPSPSClass is designed to slot into standard bioinformatics pipelines:

    # Pipe from stdin, write to stdout
    cat sequences.fasta | epspsclass classify -i - -o -

    # Snakemake rule — see workflow/rules/epspsclass.smk
    # Nextflow module  — see workflow/modules/epspsclass.nf

    # Exit codes:
    #   0 — success (including sequences with Unclassified result)
    #   1 — input/output error
    #   2 — no sequences in input
    #   3 — all sequences below identity threshold (likely wrong input)

Usage
-----
    epspsclass classify -i sequences.fasta -o results.tsv
    epspsclass classify -i sequences.fasta -o results.tsv --threshold 40
    epspsclass classify -i sequences.fasta -o results.tsv --format tsv
    epspsclass classify -i sequences.fasta -o results.tsv --format csv
    epspsclass classify -i - -o -                    # stdin → stdout
    epspsclass classify -i seqs.fasta -o out.tsv --summary  # print summary to stderr
    epspsclass download-refs                          # fetch reference sequences
    epspsclass validate-markers                       # print marker table and exit
    epspsclass --version
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path

from .classifier import EPSPSClassifier
from .ref_downloader import download_all_references

logger = logging.getLogger(__name__)

# Exit codes — documented so pipeline wrappers can act on them
EXIT_OK             = 0
EXIT_IO_ERROR       = 1
EXIT_NO_SEQUENCES   = 2
EXIT_ALL_UNRELIABLE = 3


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

FIELDNAMES = [
    "query_id",
    "primary_class",
    "all_classes",
    "sensitivity",
    "identity_I",
    "identity_II",
    "identity_III",
    "identity_IV",
    "is_unclassified",
    "notes",
]


def _row(r):
    return {
        "query_id":        r.query_id,
        "primary_class":   r.primary_class,
        "all_classes":     ";".join(r.classes) if r.classes else "None",
        "sensitivity":     r.sensitivity,
        "identity_I":      f"{r.identity_pct.get('I',  0):.2f}",
        "identity_II":     f"{r.identity_pct.get('II', 0):.2f}",
        "identity_III":    f"{r.identity_pct.get('III',0):.2f}",
        "identity_IV":     f"{r.identity_pct.get('IV', 0):.2f}",
        "is_unclassified": str(r.is_unclassified),
        "notes":           " | ".join(r.notes),
    }


def _print_summary(results, file=sys.stderr):
    from collections import Counter
    counts = Counter(r.primary_class for r in results)
    n = len(results)
    mixed = [r for r in results if len(r.classes) > 1]
    print(f"\nEPSPSClass summary ({n} sequences):", file=file)
    for cls in ["I", "II", "III", "IV", "Unclassified"]:
        pct = counts.get(cls, 0) / n * 100 if n else 0
        print(f"  Class {cls:12s}: {counts.get(cls, 0):6d}  ({pct:.1f}%)", file=file)
    if mixed:
        pct = len(mixed) / n * 100
        print(f"  Mixed (>1 class) : {len(mixed):6d}  ({pct:.1f}%)", file=file)
    unreliable = sum(1 for r in results if any("below threshold" in n for n in r.notes))
    if unreliable:
        print(f"  Flagged unreliable: {unreliable} ({unreliable/n*100:.1f}%)", file=file)


# ---------------------------------------------------------------------------
# Subcommand: classify
# ---------------------------------------------------------------------------

def cmd_classify(args: argparse.Namespace) -> int:
    clf = EPSPSClassifier(identity_threshold=args.threshold)

    # Handle stdin
    if args.input == "-":
        import io
        from Bio import SeqIO
        raw = sys.stdin.read()
        records = list(SeqIO.parse(io.StringIO(raw), "fasta"))
        results = [clf.classify(r.id, r.seq) for r in records]
    else:
        try:
            results = clf.classify_fasta(args.input)
        except FileNotFoundError:
            print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
            return EXIT_IO_ERROR

    if not results:
        print("ERROR: no sequences found in input.", file=sys.stderr)
        return EXIT_NO_SEQUENCES

    # Warn if everything is unreliable
    n_unreliable = sum(
        1 for r in results
        if all("below threshold" in n for n in r.notes) and r.notes
    )
    if n_unreliable == len(results):
        print(
            "WARNING: all sequences are below the identity threshold. "
            "Check that the input contains EPSPS protein sequences.",
            file=sys.stderr,
        )
        # Don't exit with error — still write the output file

    # Choose delimiter
    delimiter = "\t" if args.format == "tsv" else ","

    # Handle stdout
    if args.output == "-":
        out_fh = sys.stdout
        close_fh = False
    else:
        try:
            out_fh = open(args.output, "w", newline="")
            close_fh = True
        except OSError as e:
            print(f"ERROR: cannot open output file: {e}", file=sys.stderr)
            return EXIT_IO_ERROR

    try:
        writer = csv.DictWriter(out_fh, fieldnames=FIELDNAMES, delimiter=delimiter)
        writer.writeheader()
        for r in results:
            writer.writerow(_row(r))
    finally:
        if close_fh:
            out_fh.close()

    if args.summary or args.output == "-":
        _print_summary(results)
    elif args.output != "-":
        print(
            f"Classified {len(results)} sequences → {args.output}",
            file=sys.stderr,
        )

    if args.show_mixed:
        mixed = [r for r in results if len(r.classes) > 1]
        if mixed:
            print(
                f"\n--- Mixed-class sequences ({len(mixed)}) ---",
                file=sys.stderr,
            )
            for r in mixed:
                combo = ";".join(sorted(r.classes))
                ids = ", ".join(f"{k}:{v:.1f}%" for k, v in r.identity_pct.items())
                print(
                    f"  {r.query_id}: classes=[{combo}] "
                    f"primary={r.primary_class} identities=[{ids}]",
                    file=sys.stderr,
                )
        else:
            print("No mixed-class sequences found.", file=sys.stderr)

    return EXIT_OK


# ---------------------------------------------------------------------------
# Subcommand: download-refs
# ---------------------------------------------------------------------------

def cmd_download_refs(args: argparse.Namespace) -> int:
    try:
        download_all_references(force=args.force)
        print("Reference sequences ready.", file=sys.stderr)
        return EXIT_OK
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return EXIT_IO_ERROR


# ---------------------------------------------------------------------------
# Subcommand: validate-markers
# ---------------------------------------------------------------------------

def cmd_validate_markers(args: argparse.Namespace) -> int:
    from .classifier import (
        CLASS_I_MARKERS, CLASS_II_MARKERS, CLASS_IV_MARKERS,
        CAROZZI_DOMAINS,
    )
    print("=== EPSPSClass marker table ===\n")
    print(f"Class I   (vcEPSPS, {len(CLASS_I_MARKERS)} positions) — Sensitive")
    for pos, aa in sorted(CLASS_I_MARKERS.items()):
        print(f"  pos {pos:4d}: {aa}")

    print(f"\nClass II  (cbEPSPS, {len(CLASS_II_MARKERS)} positions) — Resistant")
    for pos, aa in sorted(CLASS_II_MARKERS.items()):
        print(f"  pos {pos:4d}: {aa}")

    print(f"\nClass IV  (sdEPSPS, {len(CLASS_IV_MARKERS)} positions, ≥10 required) — Resistant")
    for pos, aa in sorted(CLASS_IV_MARKERS.items()):
        print(f"  pos {pos:4d}: {aa}")

    print(f"\nClass III (bvEPSPS, {len(CAROZZI_DOMAINS)} Carozzi domains, >=1 required) — Resistant")
    for name, domain in CAROZZI_DOMAINS.items():
        print(f"  Domain {name:5s}: {len(domain)} positions")

    print(
        "\nAll positions in reference sequence coordinates (1-based).\n"
        "Derived from MAFFT L-INS-i alignment of Leino et al. (2021)\n"
        "Supplementary Table 2 reference sequences (464 alignment columns).\n"
        "Canonical sensitivity marker: vcEPSPS[101]=P / cbEPSPS[100]=L."
    )
    return EXIT_OK


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main() -> None:
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(
        prog="epspsclass",
        description=(
            "EPSPSClass — open-source reimplementation of the EPSPS "
            "glyphosate-sensitivity classifier (Leino et al. 2021).\n\n"
            "Pipeline use:\n"
            "  cat seqs.fasta | epspsclass classify -i - -o - | downstream_tool\n"
            "  epspsclass classify -i seqs.fasta -o results.tsv --summary"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.2")

    sub = parser.add_subparsers(dest="command", required=True)

    # classify
    p_clf = sub.add_parser(
        "classify",
        help="Classify EPSPS sequences from a FASTA file.",
        description=(
            "Classify protein sequences as EPSPS class I/II/III/IV.\n"
            "Use '-' for stdin/stdout to integrate with pipelines."
        ),
    )
    p_clf.add_argument(
        "-i", "--input", required=True,
        help="Input FASTA file (protein sequences). Use '-' for stdin.",
    )
    p_clf.add_argument(
        "-o", "--output", required=True,
        help="Output file path. Use '-' for stdout.",
    )
    p_clf.add_argument(
        "--threshold", type=float, default=40.0,
        help="Minimum %% identity to reference for marker checking (default: 40.0).",
    )
    p_clf.add_argument(
        "--format", choices=["tsv", "csv"], default="tsv",
        help="Output format: tsv (default) or csv.",
    )
    p_clf.add_argument(
        "--summary", action="store_true",
        help="Print a class-count summary to stderr after classification.",
    )
    p_clf.add_argument(
        "--show-mixed", action="store_true",
        help=(
            "Print sequences matching more than one class to stderr. "
            "These are sequences where markers from multiple classes are "
            "simultaneously present — biologically meaningful and reported "
            "in the original EPSPSClass tool. The primary_class column still "
            "shows the highest-priority class (IV > III > II > I); the "
            "all_classes column shows all matches."
        ),
    )
    p_clf.set_defaults(func=cmd_classify)

    # download-refs
    p_dl = sub.add_parser(
        "download-refs",
        help="Download/verify reference sequences from UniProt/NCBI.",
    )
    p_dl.add_argument(
        "--force", action="store_true",
        help="Re-download even if files already exist.",
    )
    p_dl.set_defaults(func=cmd_download_refs)

    # validate-markers
    p_vm = sub.add_parser(
        "validate-markers",
        help="Print the full marker table and exit.",
    )
    p_vm.set_defaults(func=cmd_validate_markers)

    args = parser.parse_args()
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
