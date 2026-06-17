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


class TestGroupByOrganism:
    """Test organism-level aggregation."""

    def test_single_class_organism(self):
        from epspsclass.groupby import group_by_organism

        results = [
            {"query_id": "seq1", "primary_class": "I",
             "is_too_divergent": "False"},
            {"query_id": "seq2", "primary_class": "I",
             "is_too_divergent": "False"},
        ]
        metadata = {
            "seq1": {"sequence_id": "seq1", "organism_id": "org1"},
            "seq2": {"sequence_id": "seq2", "organism_id": "org1"},
        }
        rows = group_by_organism(results, metadata)
        assert len(rows) == 1
        assert rows[0]["primary_class"] == "I"
        assert rows[0]["is_mixed"] == "False"
        assert rows[0]["n_sequences"] == "2"

    def test_mixed_organism(self):
        from epspsclass.groupby import group_by_organism

        results = [
            {"query_id": "seq1", "primary_class": "I",
             "is_too_divergent": "False"},
            {"query_id": "seq2", "primary_class": "II",
             "is_too_divergent": "False"},
        ]
        metadata = {
            "seq1": {"sequence_id": "seq1", "organism_id": "org1"},
            "seq2": {"sequence_id": "seq2", "organism_id": "org1"},
        }
        rows = group_by_organism(results, metadata)
        assert rows[0]["primary_class"] == "Mixed"
        assert rows[0]["is_mixed"] == "True"
        assert rows[0]["classes_found"] == "I;II"

    def test_unmapped_sequences_excluded(self, capsys):
        from epspsclass.groupby import group_by_organism

        results = [
            {"query_id": "seq1", "primary_class": "I",
             "is_too_divergent": "False"},
            {"query_id": "seq_orphan", "primary_class": "II",
             "is_too_divergent": "False"},
        ]
        metadata = {
            "seq1": {"sequence_id": "seq1", "organism_id": "org1"},
        }
        rows = group_by_organism(results, metadata)
        assert len(rows) == 1
        captured = capsys.readouterr()
        assert "WARNING" in captured.err

    def test_extra_metadata_passthrough(self):
        from epspsclass.groupby import group_by_organism

        results = [
            {"query_id": "seq1", "primary_class": "I",
             "is_too_divergent": "False"},
        ]
        metadata = {
            "seq1": {
                "sequence_id": "seq1",
                "organism_id": "org1",
                "exposure_tier": "Tier1",
                "habitat": "agricultural_soil",
            },
        }
        rows = group_by_organism(results, metadata)
        assert rows[0]["exposure_tier"] == "Tier1"
        assert rows[0]["habitat"] == "agricultural_soil"


class TestSummaryUnreliableCount:
    """Regression test for a bug where the CLI summary and warning miscounted
    'unreliable' sequences by string-matching on notes text instead of using
    the classifier's own is_too_divergent property.

    A sequence that fails the identity threshold against some references but
    passes against at least one, and gets a real primary classification, was
    being counted as unreliable because it still had a 'below threshold' note
    for the references it failed. is_too_divergent is False for any sequence
    that received a primary classification; only sequences that failed the
    threshold against ALL FOUR references and got no classification should
    count as unreliable.

    Caught while running this classifier on EPSPS proteins fetched from NCBI
    for 164 genera: every sequence was flagged unreliable in the printed
    summary despite the per-sequence is_too_divergent column correctly
    reading False for sequences with a valid primary_class.
    """

    def test_partial_threshold_failure_not_counted_unreliable(self, clf):
        # A real Pseudomonas EPSPS sequence. Passes identity threshold
        # against class I and class III references, fails against class II
        # and class IV, and gets a valid primary classification (III).
        seq = (
            "MSLPTVTVVPPGTPLQGRVAPPGSKSITNRALLLAALAKGTSRLTGALKSDDTRHMSEALRQMGVTIE"
            "EPDATTFVVTSSGVLQPPAGPLFLGNAGTAVRFLTAAVATVQGEVVLDGDAYMQKRPIGPLLEALRAA"
            "GVDASSPTGCPPVTVRGNGKLQVSRLEIDGGLSSQYVSALLMLAALGEAPIEVALTGKDIGARGYVD"
            "LTLACMQAFGAQVETVDDSTWRVAATGYRANDYLIEPDASAATYLWAAEALTHGRIDLGVANEAFTQ"
            "PDAQAQAHIARFPHMDAVIDGSQMQDAVPTLAVLAAFNTTPVRFTQLANLRVKECDRVQALFDGLNA"
            "IRPGLATVEGDDLLVAADPALAGTATDAQIDTHSDHRMAMCFALAGLKVSGIRIQDPACVGKTYPGY"
            "WDALRELGVELVFA"
        )
        result = clf.classify("test_pseudomonas", seq)

        # Sanity check on the scenario this test depends on: identity must
        # fail for at least one reference and pass for at least one, with a
        # real classification, or this test is not exercising the bug.
        assert result.primary_class != "Unclassified"
        assert any(pct < 40.0 for pct in result.identity_pct.values())
        assert any(pct >= 40.0 for pct in result.identity_pct.values())

        # The actual regression: this sequence must not be flagged
        # too_divergent, since it received a real classification.
        assert not result.is_too_divergent

    def test_print_summary_unreliable_count_uses_is_too_divergent(self, clf, capsys):
        from epspsclass.cli import _print_summary

        seq = (
            "MSLPTVTVVPPGTPLQGRVAPPGSKSITNRALLLAALAKGTSRLTGALKSDDTRHMSEALRQMGVTIE"
            "EPDATTFVVTSSGVLQPPAGPLFLGNAGTAVRFLTAAVATVQGEVVLDGDAYMQKRPIGPLLEALRAA"
            "GVDASSPTGCPPVTVRGNGKLQVSRLEIDGGLSSQYVSALLMLAALGEAPIEVALTGKDIGARGYVD"
            "LTLACMQAFGAQVETVDDSTWRVAATGYRANDYLIEPDASAATYLWAAEALTHGRIDLGVANEAFTQ"
            "PDAQAQAHIARFPHMDAVIDGSQMQDAVPTLAVLAAFNTTPVRFTQLANLRVKECDRVQALFDGLNA"
            "IRPGLATVEGDDLLVAADPALAGTATDAQIDTHSDHRMAMCFALAGLKVSGIRIQDPACVGKTYPGY"
            "WDALRELGVELVFA"
        )
        result = clf.classify("test_pseudomonas", seq)
        _print_summary([result])
        captured = capsys.readouterr()

        # Before the fix, this printed "Flagged unreliable: 1 (100.0%)" for
        # a sequence that actually received a valid primary classification.
        assert "Flagged unreliable" not in captured.err

    def test_cmd_classify_no_warning_for_partial_threshold_failure(self, tmp_path, capsys):
        import argparse

        from epspsclass.cli import cmd_classify

        seq = (
            "MSLPTVTVVPPGTPLQGRVAPPGSKSITNRALLLAALAKGTSRLTGALKSDDTRHMSEALRQMGVTIE"
            "EPDATTFVVTSSGVLQPPAGPLFLGNAGTAVRFLTAAVATVQGEVVLDGDAYMQKRPIGPLLEALRAA"
            "GVDASSPTGCPPVTVRGNGKLQVSRLEIDGGLSSQYVSALLMLAALGEAPIEVALTGKDIGARGYVD"
            "LTLACMQAFGAQVETVDDSTWRVAATGYRANDYLIEPDASAATYLWAAEALTHGRIDLGVANEAFTQ"
            "PDAQAQAHIARFPHMDAVIDGSQMQDAVPTLAVLAAFNTTPVRFTQLANLRVKECDRVQALFDGLNA"
            "IRPGLATVEGDDLLVAADPALAGTATDAQIDTHSDHRMAMCFALAGLKVSGIRIQDPACVGKTYPGY"
            "WDALRELGVELVFA"
        )
        in_path = tmp_path / "test.faa"
        in_path.write_text(f">test_pseudomonas\n{seq}\n")
        out_path = tmp_path / "out.tsv"

        args = argparse.Namespace(
            input=str(in_path),
            output=str(out_path),
            threshold=40.0,
            format="tsv",
            summary=True,
            show_mixed=False,
        )

        # capsys captures sys.stdout/sys.stderr at the point pytest reads it,
        # which works correctly with _print_summary's late-bound sys.stderr
        # default, unlike contextlib.redirect_stderr.
        cmd_classify(args)
        captured = capsys.readouterr()

        # Before the fix, this warning fired for a single-sequence input
        # that failed the identity threshold against any one reference,
        # even with a valid classification from the others.
        assert "WARNING: all sequences are below the identity threshold" not in captured.err
        assert "Flagged unreliable" not in captured.err


class TestClassICoreMarkerThreshold:
    """Regression tests for the Class I classification fix.

    Class I previously required all 148 CLASS_I_MARKERS positions to match
    exactly, which meant no real, evolutionarily diverged organism was ever
    classified as class I. An Ia/Ib marker split was suspected as the cause
    and checked against Leino et al.'s real Supplementary Table 1; the Ia
    and Ib columns are identical at all 464 positions, so there was no
    split to recover. A simple count threshold on the full 148 markers was
    also tried and rejected: real class II organisms found 18-22 of 148,
    only a 2-marker gap below the real class I floor of 24-28, too narrow
    to threshold safely.

    The actual fix selects the 20 most discriminating positions from
    CLASS_I_MARKERS (CLASS_I_CORE_MARKERS) based on real benchmark testing,
    and thresholds on that subset instead (CLASS_I_CORE_THRESHOLD = 4).

    These are the exact sequences used to derive and validate that
    threshold, fetched from NCBI (esearch/esummary/efetch, aroA or EPSP
    synthase annotated, RefSeq preferred). Hardcoded here so this test does
    not require network access to run. If this benchmark is expanded,
    rerun scripts/epsps/classify_and_calibrate.py and update both the
    CLASS_I_CORE_MARKERS derivation comment and these sequences together.
    """

    # Confirmed class II organisms (resistant to glyphosate).
    BENCHMARK_CLASS_II = {
        "Staphylococcus_aureus|YFF31137":
            "MVNEQIIDISGPLKGEIEVPGDKSMTHRAIMLASLAEGVSTIYKPLLGEDCRRTMDIFRLLGVEIKEDDEKLVVTSPGYQSFNTPHQVLYTGNSGTTTRLLAGLLSGLGIESVLSGDVSIGKRPMDRVLRPLKSMNANIEGIEDNYTPLIIKPSVIKGINYKMEVASAQVKSAILFASLFSKEPTIIKELDVSRNHTETMFRHFNIPIEAEGLSITTTPEAIRYIKPADFHVPGDISSAAFFIVAALITPGSDVTIHNVGINPTRSGIIDIVEKMGGNIQLFNQTTGAEPTASIRIQYTPMLQPITIEGELVPKAIDELPVIALLCTQAVGTSTIKDAEELKVKETNRIDTTADMLNLLGFELQPTNDGLIIHPSEFKTNATVDSLTDHRIGMMLAVASLLSSEPVKIKQFDAVNVSFPGFLPKLKLLENEG",
        "Ruminococcus_gnavus|WP_268803945":
            "MELTSITGLKGEVSIPGDKSISHRSVMFASLANGMTEIHNFLNGADCLATIDCFRKMGIKIEEHQNRILVHGKGLHGLTAPSETLQVKNSGTTTRLLSGILAGQPFSTSLSGDESLNSRPMKRIIEPLTQMGAHISSLHGNGCAPLIIEPGKLHGIHYTSPVASAQVKSCILLAGLYAEGETSVTEPILSRNHTELMLKEFGADIRTVHQLAGSEATSVIQPCPELHGQKITVPGDISSAAYFIAAGLLVPDSEILVKNVGINPTRAGLLKVCEDMGGNITLLNERTEAGEKMADILVRSSQLHGISIHGDIIPTLIDEIPIIAVMAACAEGTTIIRDAQELRVKETDRIETITDNLIAMGCSVLPTEDGMVIKGGEPLKGATIHTLLDHRIAMAFSIAALVADGRTKILDSHCIDVSYPGFYDAFEHLL",
        "Dorea_formicigenerans|WP_390371098":
            "MSDRIICPCKGLHGEIMIPGDKSISHRSIMLGALALGTTEITNFLEGADCLSTIGCFQSMGIQIDRTPEKIIVHGKGMHGLSAPKDILNVGNSGTTTRLMSGILSAQDFTSVMSGDASLNSRPMGRVITPLTQMGAHITSVNGDLCAPLKIEPGMLHGIDYTSPVASAQVKSAILLAGLYADGETSVTEPALSRNHTELMLKSFGADITSTVNPDGTATAHVKPCQELYGQSICVPGDISSAAYFIAAGLLTLDSELLVKNVGINKTRAGFLEVCQNMGADITLVNESLEGGEPRADILVRTSKLHGTTIEGALIPTLIDEIPMIAVMAACAEGTTIIKDAAELKVKETNRIDTTTEALRSMGTDITPTDDGMIIQGGHTLHGAKINSYLDHRIAMAFAIAALSADGDTIIHDSQCVDVSYPEFFEILDGCR",
    }

    # Confirmed class I organisms (sensitive to glyphosate), from
    # phylogenetically distant taxa (Bacteroidetes, Firmicutes), the same
    # taxa where this issue was originally observed.
    BENCHMARK_CLASS_I = {
        "Bacteroides_fragilis|WP_220391824":
            "MRYLLSAPSHIKATIQLPASKSISNRALIIHALSKGNDVLSNLSDCDDTRVMVKALTEGGEVIDILAAGTAMRFLTAYLSSTPGTHIITGTERMQQRPIQILVNALRELGASIEYTRNEGFPPLRIEGAPLAGNEITLKGNVSSQYISALLMIGPILKNGLQLRLTGEVVSRPYINLTLQLMKDFGASARWTSDQSISVEPEPYRCVPFTVESDWSAASYWYQMAALSSEADIELTGLFRHSYQGDSRGAEVFARLGIETEYTEEGIRLRKNGSYVKRLDEDFVDIPDLAQTFVVTCALLDVPFRFTGLQSLKIKETDRIEALKAEMKKLGYVLHDEDNSILYWNGERIEPQACPVIKTYEDHRMAMAFAPAAIHYPTIQIDEPQVVSKSYPGYWDDLRKAGFMIEAHVTES",
        "Bacteroides_thetaiotaomicron|WP_434276370":
            "MLYKLISPSMVKATIQLPASKSISNRALIINALGKGIYPPENLSDCDDTQVMIKALTEGKETIDIMAAGTAMRFLTAYLSATSGERIITGTARMQQRPIQILVNALRELGAEIEYTHNEGYPPLRIKGAELKGNEITLKGNVSSQYISALLMIGPVLKDGLTLHLTGEIISRPYINLTLQLMQDFGAKAAWTSPSSISVAPQPYQSVPFTVESDWSAASYWYQIAALSPEAEIELLGLFRNSYQGDSRGAEVFSRLGITTEFTPQGVKIKKTGKAPERLEEDFVDIPDLAQTFVVTCALLNIPFRFTGLQSLKIKETDRIAALRTELKKLGYLIEEENDSVLMWNGERCEPEAVPVIATYEDHRMAMAFAPAVITFPKLLIADPQVVSKSYPGYWEDLKLARFQVINEG",
        "Clostridium_perfringens|WP_483510405":
            "MKKVIITPSKLKGSVKIPPSKSMAHRAIICASLSKGESVISNIDFSEDIIATMEGMKSLGANIKVEKDKLIINGENILKDSNYKVIDCNESGSTLRFLVPISLIKDNRVNFIGRGNLGKRPLKTYYEIFEEQEIKYSYEEENLDLNIEGSLKGGEFKVKGNISSQFISGLLFTLPLLKEDSKIIITTELESKGYIDLTLDMIEKFGVTIKNNNYREFLIKGNQSYKPMNYKVEGDYSQAAFYFSAGALGSEINCLDLDLSSYQGDKECIEILEGMGARLIESQKRSLSIIHGDLNGTIIDASQCPDIIPVLTVVAALSKGETRIINGERLRIKECDRLNAIFTELNKLGADIKELKDGLIINGVKDLIGGEVYSHKDHRIAMSLAIASTRCKEEVIIKEPDCVKKSYPGFWEDFKSLGGILREE",
        "Prevotella_copri|EFB34131":
            "MKYTIKAPRQLNASINLPASKSISNRALVINAMAGCKLQPRNLSDCDDTEVIIAALRDMPDVINIKAAGTAMRFMTAYLSATPGEHTITGTERMQNRPIAILVDALRYLGADIQYEKKEGYPPLHIVGKPLEGGHLEVVGNISSQYISALLMIGPILKNGLELKLTGEIASRPYIDLTLWTMQNFGASAEWTDVDTITVKPQPYSCVADYTIENDWSASSYWYEMMALNGNPDSEVRLEGLFDSSKQGDSVVKYIFSLLGVKSEFENRDVLSPVKLKVQRCLLPRFDYDFSGSPDLAQTIVVACCALGVKFKFTGLASLKIKETDRIEALKKELKKVGYVIYDENDNTLIWDGETCEPSFEPIDTYEDHRMALAFAPLAFKFPQIEINNPEVVSKSYPHYWEDLKKVGFEIVES",
    }

    def test_benchmark_class_i_organisms_classify_as_class_i(self, clf):
        for name, seq in self.BENCHMARK_CLASS_I.items():
            result = clf.classify(name, seq)
            assert result.primary_class == "I", (
                f"{name} should classify as Class I, got {result.primary_class}. "
                f"identities={result.identity_pct}"
            )
            assert result.sensitivity == "Sensitive"

    def test_benchmark_class_ii_organisms_do_not_classify_as_class_i(self, clf):
        for name, seq in self.BENCHMARK_CLASS_II.items():
            result = clf.classify(name, seq)
            assert "I" not in result.classes, (
                f"{name} incorrectly matched class I markers. "
                f"classes={result.classes}, identities={result.identity_pct}"
            )

    def test_benchmark_class_ii_organisms_still_classify_as_class_ii(self, clf):
        for name, seq in self.BENCHMARK_CLASS_II.items():
            result = clf.classify(name, seq)
            assert result.primary_class == "II", (
                f"{name} should classify as Class II, got {result.primary_class}. "
                f"identities={result.identity_pct}"
            )
            assert result.sensitivity == "Resistant"
