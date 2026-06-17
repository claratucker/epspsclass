#!/usr/bin/env python3
# classify_and_calibrate.py
#
# Reads results/epsps/class1_calibration_benchmark.faa (built by
# calibrate_class1_subset.py), aligns each benchmark sequence against the
# vcEPSPS reference using epspsclass's own pairwise alignment function, and
# reports per-organism counts against the existing 148-position
# CLASS_I_MARKERS set. Then tests whether scoring against a SUBSET of those
# 148 positions, rather than the full set, gives a cleaner separation
# between the benchmark class I and class II organisms.
#
# Subset selection method: for each of the 148 positions, check whether the
# benchmark class I organisms match that position's expected residue more
# often than the benchmark class II organisms do. Positions where this
# difference is large are more discriminating; positions where class I and
# class II organisms match about equally often add noise rather than
# signal. This mirrors how Leino et al. (2021) single out one canonical
# position (101, Pro106 in E. coli numbering) as the primary discriminator
# rather than weighting all 148 positions equally (Funke et al. 2009;
# Sammons and Gaines 2014).
#
# This script only measures and reports. It does not modify classifier.py.
# Review the output before deciding whether any subset gives a safe,
# explainable threshold.

import os
import sys

# Import the local repo's epspsclass package (this script lives in
# scripts/epsps/ inside the repo, so the package is two directories up).
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _REPO_ROOT)

from Bio import SeqIO
from epspsclass.classifier import (
    CLASS_I_MARKERS,
    _pairwise_align,
    _ref_pos_to_aligned,
    _load_reference,
)

BENCHMARK_FAA = os.path.join(_REPO_ROOT, "results", "class1_calibration_benchmark.faa")


def score_against_class_i(seq):
    """Align seq against vcEPSPS and return {position: matched_bool} for
    every position in CLASS_I_MARKERS, using the classifier's own alignment
    and position-mapping functions so results are directly comparable to
    what the live classifier computes."""
    ref = _load_reference("I")
    aln = _pairwise_align(seq, ref)
    aq, ar = str(aln.seqA), str(aln.seqB)
    ref_map = _ref_pos_to_aligned(ar)

    matches = {}
    for pos, expected_aa in CLASS_I_MARKERS.items():
        aligned_col = ref_map.get(pos)
        if aligned_col is None:
            matches[pos] = False
            continue
        matches[pos] = aq[aligned_col] == expected_aa
    return matches


def main():
    if not os.path.exists(BENCHMARK_FAA):
        sys.exit(f"Benchmark FASTA not found at {BENCHMARK_FAA}. Run calibrate_class1_subset.py first.")

    records = list(SeqIO.parse(BENCHMARK_FAA, "fasta"))
    if not records:
        sys.exit("Benchmark FASTA is empty.")

    class_i_records = [r for r in records if "benchmark_classI_" in r.id]
    class_ii_records = [r for r in records if "benchmark_classII_" in r.id]

    print(f"Benchmark: {len(class_i_records)} class I organisms, {len(class_ii_records)} class II organisms\n")

    all_match_results = {}
    for r in records:
        matches = score_against_class_i(str(r.seq))
        all_match_results[r.id] = matches
        n_matched = sum(matches.values())
        print(f"{r.id}: {n_matched}/{len(CLASS_I_MARKERS)} class I markers matched")

    print("\n--- Per-position discrimination ---")
    print("For each position, compare match rate in class I vs class II benchmark organisms.")
    position_scores = []
    for pos in CLASS_I_MARKERS:
        i_matches = [all_match_results[r.id][pos] for r in class_i_records]
        ii_matches = [all_match_results[r.id][pos] for r in class_ii_records]
        i_rate = sum(i_matches) / len(i_matches) if i_matches else 0
        ii_rate = sum(ii_matches) / len(ii_matches) if ii_matches else 0
        discrimination = i_rate - ii_rate
        position_scores.append((pos, i_rate, ii_rate, discrimination))

    position_scores.sort(key=lambda x: -x[3])
    print("\nTop 20 most discriminating positions (class I match rate minus class II match rate):")
    print(f"{'pos':>5} {'classI_rate':>12} {'classII_rate':>13} {'discrimination':>15}")
    for pos, i_rate, ii_rate, disc in position_scores[:20]:
        print(f"{pos:>5} {i_rate:>12.2f} {ii_rate:>13.2f} {disc:>15.2f}")

    print("\nBottom 10 (least discriminating, or inverted):")
    for pos, i_rate, ii_rate, disc in position_scores[-10:]:
        print(f"{pos:>5} {i_rate:>12.2f} {ii_rate:>13.2f} {disc:>15.2f}")

    for subset_size in [10, 20, 30, 50, 148]:
        subset = [p for p, _, _, _ in position_scores[:subset_size]]
        print(f"\n--- Testing subset of top {subset_size} discriminating positions ---")
        i_counts, ii_counts = [], []
        for r in class_i_records:
            n = sum(1 for p in subset if all_match_results[r.id][p])
            i_counts.append(n)
            print(f"  class I  {r.id}: {n}/{subset_size}")
        for r in class_ii_records:
            n = sum(1 for p in subset if all_match_results[r.id][p])
            ii_counts.append(n)
            print(f"  class II {r.id}: {n}/{subset_size}")
        if i_counts and ii_counts:
            gap = min(i_counts) - max(ii_counts)
            print(f"  Gap (min class I count - max class II count): {gap}")
            print(f"  {'SAFE TO THRESHOLD' if gap > 0 else 'NOT SAFE, overlapping ranges'}")

    print("\nReview the gaps above. A subset is only usable for thresholding")
    print("if the gap is positive and comfortably wide, not a 1-2 marker margin.")
    print("If no subset size gives a safe gap, the count-threshold approach")
    print("does not work for Class I and a different method is needed.")


if __name__ == "__main__":
    main()
