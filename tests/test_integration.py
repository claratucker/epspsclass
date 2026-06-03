"""
tests/test_integration.py
=========================
Integration tests for EPSPSClass classifier using the bundled reference
sequences. These tests verify that each reference sequence self-classifies
correctly — the most basic validation that the marker tables, alignment,
and classification logic are internally consistent.

These tests require the bundled reference FASTA files to be present in
epspsclass/data/reference_sequences/ (included with the package).
"""

import pytest
from Bio import SeqIO
from pathlib import Path

from epspsclass.classifier import EPSPSClassifier, REF_FASTA

DATA_DIR = Path(__file__).parent.parent / "epspsclass" / "data" / "reference_sequences"


def _load_ref(cls: str):
    path = REF_FASTA[cls]
    records = list(SeqIO.parse(str(path), "fasta"))
    assert records, f"Empty or missing reference FASTA for class {cls}"
    return records[0]


@pytest.fixture(scope="module")
def clf():
    return EPSPSClassifier(identity_threshold=40.0)


class TestSelfClassification:
    """Each reference sequence must classify as its own class."""

    def test_class_i_reference(self, clf):
        rec = _load_ref("I")
        result = clf.classify(rec.id, rec.seq)
        assert result.primary_class == "I", (
            f"vcEPSPS (class I reference) classified as {result.primary_class}. "
            f"identities={result.identity_pct}"
        )
        assert result.sensitivity == "Sensitive"
        assert not result.is_unclassified

    def test_class_ii_reference(self, clf):
        rec = _load_ref("II")
        result = clf.classify(rec.id, rec.seq)
        assert result.primary_class == "II", (
            f"cbEPSPS (class II reference) classified as {result.primary_class}. "
            f"identities={result.identity_pct}"
        )
        assert result.sensitivity == "Resistant"
        assert not result.is_unclassified

    def test_class_iii_reference(self, clf):
        rec = _load_ref("III")
        result = clf.classify(rec.id, rec.seq)
        assert result.primary_class == "III", (
            f"bvEPSPS (class III reference) classified as {result.primary_class}. "
            f"identities={result.identity_pct}"
        )
        assert result.sensitivity == "Resistant"
        assert not result.is_unclassified

    def test_class_iv_reference(self, clf):
        rec = _load_ref("IV")
        result = clf.classify(rec.id, rec.seq)
        assert result.primary_class == "IV", (
            f"sdEPSPS (class IV reference) classified as {result.primary_class}. "
            f"identities={result.identity_pct}"
        )
        assert result.sensitivity == "Resistant"
        assert not result.is_unclassified


class TestCrossClassification:
    """Reference sequences should NOT classify as each other's class."""

    def test_class_i_not_class_ii(self, clf):
        rec = _load_ref("I")
        result = clf.classify(rec.id, rec.seq)
        assert "II" not in result.classes, (
            "vcEPSPS incorrectly matched class II markers"
        )

    def test_class_ii_not_class_i(self, clf):
        rec = _load_ref("II")
        result = clf.classify(rec.id, rec.seq)
        assert "I" not in result.classes, (
            "cbEPSPS incorrectly matched class I markers"
        )

    def test_class_i_not_class_iii(self, clf):
        """Class I reference must not trigger any Carozzi class III domains."""
        rec = _load_ref("I")
        result = clf.classify(rec.id, rec.seq)
        assert "III" not in result.classes, (
            "vcEPSPS incorrectly matched class III Carozzi domains — "
            "this would indicate a false positive in the domain set."
        )

    def test_class_i_not_class_iv(self, clf):
        rec = _load_ref("I")
        result = clf.classify(rec.id, rec.seq)
        assert "IV" not in result.classes


class TestIdentityValues:
    """Sanity checks on reported identity percentages."""

    def test_class_i_self_identity_is_100(self, clf):
        rec = _load_ref("I")
        result = clf.classify(rec.id, rec.seq)
        assert result.identity_pct.get("I", 0) == pytest.approx(100.0, abs=0.1)

    def test_class_ii_self_identity_is_100(self, clf):
        rec = _load_ref("II")
        result = clf.classify(rec.id, rec.seq)
        assert result.identity_pct.get("II", 0) == pytest.approx(100.0, abs=0.1)

    def test_cross_identity_below_self(self, clf):
        """A reference should always have higher identity to itself than to others."""
        for cls in ["I", "II", "III", "IV"]:
            rec = _load_ref(cls)
            result = clf.classify(rec.id, rec.seq)
            self_id = result.identity_pct.get(cls, 0)
            others = [v for k, v in result.identity_pct.items() if k != cls]
            assert all(self_id >= o for o in others), (
                f"Class {cls} reference has higher identity to another class "
                f"than to itself: self={self_id:.1f}%, others={others}"
            )


class TestNewFeatures:
    """Test features added in v1.0.1 and v1.0.2."""

    def test_no_mixed_class_in_references(self, clf):
        """None of the four reference sequences should be mixed-class."""
        for cls in ["I", "II", "III", "IV"]:
            rec = _load_ref(cls)
            result = clf.classify(rec.id, rec.seq)
            assert len(result.classes) == 1, (
                f"Class {cls} reference matched multiple classes: {result.classes}"
            )

    def test_is_too_divergent_false_for_references(self, clf):
        """Reference sequences are never too divergent."""
        for cls in ["I", "II", "III", "IV"]:
            rec = _load_ref(cls)
            result = clf.classify(rec.id, rec.seq)
            assert not result.is_too_divergent

    def test_nonstandard_aa_detected(self, clf):
        """Non-standard amino acid codes should generate a note."""
        from Bio.Seq import Seq
        # X is a common ambiguity code in database sequences
        seq_with_x = Seq("MXSPRQ" + "A" * 400)
        result = clf.classify("test_nonstandard", seq_with_x)
        assert any("Non-standard" in n for n in result.notes), (
            "Non-standard amino acid X not detected in notes"
        )

    def test_fasta_error_isolation(self, clf, tmp_path):
        """A bad sequence in a batch should not abort classification of others."""
        from Bio.Seq import Seq
        # Create a FASTA with one valid-ish sequence and one empty sequence
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">good\nMSPRQITL" + "A" * 400 + "\n>empty\n\n")
        # Should not raise
        results = clf.classify_fasta(str(fasta))
        assert len(results) >= 1
