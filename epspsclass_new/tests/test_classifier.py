"""
tests/test_classifier.py
========================
Unit tests for EPSPSClass classifier logic.

These tests use synthetic marker sequences so they do not require
internet access or the reference FASTA files.  Integration tests
(requiring downloaded references) are in tests/test_integration.py.
"""

import pytest
from unittest.mock import patch, MagicMock
from Bio.Seq import Seq

from epspsclass.classifier import (
    _percent_identity,
    _ref_pos_to_aligned,
    _check_markers,
    _check_class_iii_domains,
    CLASS_I_MARKERS,
    CLASS_II_MARKERS,
    CAROZZI_DOMAINS,
    ClassificationResult,
)


# ---------------------------------------------------------------------------
# _percent_identity
# ---------------------------------------------------------------------------

class TestPercentIdentity:
    def test_identical(self):
        assert _percent_identity("ACDEF", "ACDEF") == 100.0

    def test_no_match(self):
        assert _percent_identity("AAAAA", "CCCCC") == 0.0

    def test_gaps_ignored(self):
        # Gap columns are excluded from both numerator and denominator
        pct = _percent_identity("A-C", "A-C")
        assert pct == 100.0

    def test_partial(self):
        pct = _percent_identity("ACDEF", "ACXXX")
        assert abs(pct - 40.0) < 0.01

    def test_all_gaps(self):
        assert _percent_identity("---", "---") == 0.0


# ---------------------------------------------------------------------------
# _ref_pos_to_aligned
# ---------------------------------------------------------------------------

class TestRefPosToAligned:
    def test_no_gaps(self):
        mapping = _ref_pos_to_aligned("ACDEF")
        assert mapping == {1: 0, 2: 1, 3: 2, 4: 3, 5: 4}

    def test_with_gaps(self):
        # ref: A-CDE  → ref positions 1,2,3,4 at cols 0,2,3,4
        mapping = _ref_pos_to_aligned("A-CDE")
        assert mapping[1] == 0
        assert mapping[2] == 2
        assert 1 not in mapping.values()  # col 1 is a gap, not in mapping values

    def test_leading_gap(self):
        mapping = _ref_pos_to_aligned("-ACDE")
        assert mapping[1] == 1  # first ref residue at col 1


# ---------------------------------------------------------------------------
# _check_markers
# ---------------------------------------------------------------------------

class TestCheckMarkers:
    def _make_alignment(self, ref_seq: str, query_seq: str):
        """Helper: aligned ref and query (same length, no gaps for simplicity)."""
        return query_seq, ref_seq

    def test_all_markers_present(self):
        # Build a synthetic alignment where CLASS_I_MARKERS positions are set
        # We test with a small custom marker dict to avoid dependency on
        # actual reference coordinates
        markers = {1: "A", 2: "C"}
        # ref:   A C ...
        # query: A C ...
        ref   = "AC" + "X" * 100
        query = "AC" + "X" * 100
        ok, found = _check_markers(query, ref, markers)
        assert ok is True
        assert found == {1: "A", 2: "C"}

    def test_marker_missing(self):
        markers = {1: "A", 2: "C"}
        ref   = "AC" + "X" * 100
        query = "AX" + "X" * 100  # pos 2 wrong
        ok, found = _check_markers(query, ref, markers)
        assert ok is False
        assert 1 in found
        assert 2 not in found

    def test_position_out_of_range(self):
        markers = {999: "A"}  # position far beyond alignment length
        ref   = "ACDEF"
        query = "ACDEF"
        ok, found = _check_markers(query, ref, markers)
        assert ok is False


# ---------------------------------------------------------------------------
# ClassificationResult
# ---------------------------------------------------------------------------

class TestClassificationResult:
    def _make(self, classes):
        return ClassificationResult(
            query_id="test",
            classes=classes,
            identity_pct={"I": 90.0, "II": 30.0, "III": 25.0, "IV": 20.0},
            aligned_markers={},
            is_unclassified=(len(classes) == 0),
        )

    def test_primary_class_I(self):
        r = self._make(["I"])
        assert r.primary_class == "I"
        assert r.sensitivity == "Sensitive"

    def test_primary_class_II(self):
        r = self._make(["II"])
        assert r.primary_class == "II"
        assert r.sensitivity == "Resistant"

    def test_primary_class_IV_priority(self):
        # IV should take priority over I or II when multiple match
        r = self._make(["I", "IV"])
        assert r.primary_class == "IV"

    def test_unclassified(self):
        r = self._make([])
        assert r.primary_class == "Unclassified"
        assert r.sensitivity == "Unknown"
        assert r.is_unclassified is True

    def test_str_repr(self):
        r = self._make(["I"])
        s = str(r)
        assert "test" in s
        assert "Sensitive" in s


# ---------------------------------------------------------------------------
# Integration smoke test (skipped if refs not downloaded)
# ---------------------------------------------------------------------------

def test_classify_requires_refs(tmp_path):
    """Classifier raises FileNotFoundError if reference files are missing."""
    from epspsclass.classifier import EPSPSClassifier
    clf = EPSPSClassifier()
    # Patch _load_reference to raise FileNotFoundError
    with patch("epspsclass.classifier._load_reference",
               side_effect=FileNotFoundError("missing")):
        with pytest.raises(FileNotFoundError):
            clf.classify("q1", Seq("MSPRQ" * 80))
