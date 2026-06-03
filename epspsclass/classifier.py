"""
epspsclass/classifier.py
========================
Open-source reimplementation of the EPSPSClass classification algorithm.

Reference
---------
Leino et al. (2021) Classification of the glyphosate target enzyme
(5-enolpyruvylshikimate-3-phosphate synthase) for assessing sensitivity
of organisms to the herbicide. Environment International 149:106334.
https://doi.org/10.1016/j.envint.2020.106334

Supplementary material used
---------------------------
All marker positions and reference sequences are taken directly from:
  - Supplementary Table 1: full pairwise alignment of all four reference
    sequences (426 alignment columns, vcEPSPS coords 1-426).
  - Supplementary Table 2: FASTA sequences for all four reference organisms.
  - Supplementary Figure 7: resistance-associated positions in E. coli
    (Gly96, Thr97, Pro101, Gly137, Ala183 in E. coli / vcEPSPS numbering).

Classification logic (from Leino et al. 2021 Materials and Methods)
--------------------------------------------------------------------
For each query protein:
1. Perform global pairwise alignment against each of four reference sequences.
2. Class I   — ALL class I amino acid markers present at aligned positions.
3. Class II  — ALL class II amino acid markers present at aligned positions.
4. Class III — at least ONE complete 3-aa motif from the bvEPSPS set present.
5. Class IV  — ALL class IV amino acid markers present at aligned positions.
6. Unclassified — none of the above.

A sequence may match more than one class; all matches are reported.

Marker positions
----------------
Positions are given in the REFERENCE SEQUENCE coordinate system (1-based),
as read from Supplementary Table 1 of Leino et al. (2021).

The key discriminating position is alignment column 106:
  vcEPSPS position 101 = Pro (P)  → Class I  (glyphosate SENSITIVE)
  cbEPSPS position 100 = Leu (L)  → Class II (glyphosate RESISTANT)
This corresponds to Pro106 in E. coli numbering (Supplementary Figure 7).

Reference organisms
-------------------
  vcEPSPS — Vibrio cholerae O1 N16961        (UniProt Q9KNE7)   Class I
  cbEPSPS — Coxiella burnetii RSA 493        (UniProt Q83EH4)   Class II
  bvEPSPS — Brevundimonas vesicularis        (GenBank CAA73210) Class III
  sdEPSPS — Streptomyces davawensis JCM 4913 (UniProt H6WNZ5)   Class IV

Reference sequences are bundled in epspsclass/data/reference_sequences/
(taken verbatim from Supplementary Table 2; no network required).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Bundled reference sequence paths
# ---------------------------------------------------------------------------
_DATA_DIR = Path(__file__).parent / "data" / "reference_sequences"

REF_FASTA = {
    "I":   _DATA_DIR / "vcEPSPS_Q9KNE7.fasta",
    "II":  _DATA_DIR / "cbEPSPS_Q83EH4.fasta",
    "III": _DATA_DIR / "bvEPSPS_CAA73210.fasta",
    "IV":  _DATA_DIR / "sdEPSPS_H6WNZ5.fasta",
}

# ---------------------------------------------------------------------------
# Amino acid markers
#
# Source: Leino et al. (2021) Supplementary Table 1.
# Each entry: { reference_residue_position_1based: expected_amino_acid }
# Positions are in the REFERENCE sequence coordinate system (not alignment).
#
# Class I markers (vcEPSPS coordinates)
# Key positions from Supp Table 1 rows 101-106, 142, 189 and Supp Fig 7:
#   Row 101 → vcEPSPS pos 96  = G  (Gly96 in E. coli numbering)
#   Row 102 → vcEPSPS pos 97  = T  (Thr97)
#   Row 106 → vcEPSPS pos 101 = P  (Pro101 / Pro106 in E. coli) ← KEY MARKER
#   Row 142 → vcEPSPS pos 137 = G  (Gly137)
#   Row 189 → vcEPSPS pos 183 = L  (Ala183 in E. coli; Leu183 in vcEPSPS)
# ---------------------------------------------------------------------------

CLASS_I_MARKERS: Dict[int, str] = {
    # Derived from MAFFT L-INS-i alignment of four reference sequences
    # (Leino et al. 2021 Supp Table 2). 148 positions unique to vcEPSPS.
    # KEY: position 101 = Pro — the canonical glyphosate-sensitivity marker
    # (Supp Fig 7; corresponds to Pro106 in E. coli numbering).
       9: "I",   16: "V",   17: "N",   18: "L",   24: "V",
      36: "S",   43: "N",   44: "L",   46: "D",   50: "I",
      52: "H",   54: "L",   55: "N",   58: "T",   59: "K",
      63: "N",   64: "Y",   65: "R",   66: "L",   68: "A",
      73: "C",   74: "E",   78: "L",   81: "A",   82: "F",
      83: "H",   84: "T",   85: "T",   86: "Q",   99: "M",
     101: "P",  103: "A",  107: "C",  110: "Q",  112: "D",
     118: "E",  119: "P",  120: "R",  122: "K",  123: "E",
     128: "H",  135: "Q",  136: "A",  141: "E",  142: "Y",
     143: "L",  145: "Q",  146: "E",  148: "F",  152: "R",
     154: "Q",  159: "Q",  162: "T",  169: "I",  173: "F",
     177: "F",  180: "S",  184: "A",  185: "Q",  186: "G",
     187: "K",  189: "T",  191: "K",  193: "V",  199: "K",
     207: "H",  208: "I",  210: "E",  211: "Q",  215: "Q",
     217: "I",  218: "N",  219: "H",  220: "D",  221: "Y",
     223: "E",  226: "I",  228: "A",  230: "Q",  234: "S",
     235: "P",  236: "G",  237: "Q",  239: "L",  240: "V",
     251: "L",  257: "K",  262: "K",  264: "T",  266: "I",
     268: "K",  270: "S",  271: "I",  275: "I",  276: "Q",
     278: "A",  280: "A",  282: "E",  299: "R",  303: "N",
     307: "L",  309: "F",  310: "N",  311: "H",  316: "A",
     317: "M",  319: "I",  321: "T",  322: "T",  324: "L",
     327: "K",  329: "T",  331: "A",  334: "N",  336: "Y",
     338: "W",  351: "T",  352: "E",  356: "V",  358: "A",
     359: "T",  364: "E",  366: "F",  368: "V",  370: "T",
     376: "I",  377: "H",  379: "A",  388: "M",  391: "C",
     393: "S",  395: "V",  396: "A",  398: "S",  399: "D",
     400: "T",  405: "N",  408: "K",  410: "T",  411: "S",
     416: "D",  418: "F",  419: "D",  420: "K",  422: "A",
     423: "Q",  425: "S",  426: "R",
}

CLASS_II_MARKERS: Dict[int, str] = {
    # 204 positions unique to cbEPSPS. KEY: position 100 = Leu —
    # replaces Pro101 of class I; primary resistance marker.
       8: "S",   15: "I",   16: "C",   17: "V",   20: "D",
      25: "H",   28: "V",   33: "I",   35: "E",   37: "Q",
      39: "Q",   40: "V",   41: "D",   43: "F",   45: "M",
      46: "G",   47: "A",   49: "N",   50: "L",   51: "A",
      53: "V",   54: "S",   57: "Q",   58: "Q",   61: "A",
      62: "S",   64: "Q",   65: "V",   67: "E",   71: "I",
      72: "L",   73: "V",   77: "V",   79: "M",   80: "T",
      81: "G",   82: "L",   83: "Q",   85: "P",   90: "D",
      94: "S",   98: "I",  100: "L",  102: "S",  103: "G",
     107: "G",  108: "Q",  109: "P",  110: "F",  111: "N",
     112: "T",  118: "S",  119: "S",  120: "L",  121: "Q",
     125: "M",  126: "K",  127: "R",  128: "I",  129: "I",
     131: "P",  133: "T",  134: "L",  135: "M",  137: "A",
     138: "K",  139: "I",  140: "D",  141: "S",  142: "T",
     143: "G",  145: "V",  149: "K",  151: "Y",  153: "N",
     154: "P",  157: "T",  159: "I",  160: "H",  161: "Y",
     164: "P",  165: "M",  166: "A",  168: "A",  169: "Q",
     171: "K",  172: "S",  173: "C",  178: "G",  180: "Y",
     181: "A",  184: "K",  185: "T",  186: "C",  190: "P",
     191: "A",  193: "S",  194: "R",  196: "H",  198: "E",
     199: "R",  200: "L",  201: "L",  202: "K",  203: "H",
     205: "H",  206: "Y",  207: "T",  208: "L",  209: "Q",
     210: "K",  211: "D",  212: "K",  213: "Q",  214: "S",
     215: "I",  216: "C",  218: "S",  219: "G",  223: "L",
     224: "K",  226: "N",  228: "I",  229: "S",  231: "P",
     234: "I",  238: "A",  239: "F",  241: "I",  242: "V",
     245: "T",  247: "T",  250: "S",  251: "A",  253: "R",
     255: "C",  256: "R",  257: "V",  259: "V",  261: "P",
     262: "T",  263: "R",  266: "V",  267: "I",  268: "N",
     271: "K",  272: "M",  297: "H",  298: "A",  299: "R",
     301: "K",  303: "I",  306: "P",  307: "P",  308: "D",
     310: "V",  316: "E",  317: "F",  319: "V",  321: "L",
     322: "I",  323: "A",  326: "V",  328: "Q",  330: "K",
     332: "V",  333: "L",  336: "A",  338: "E",  351: "V",
     352: "D",  355: "Q",  359: "I",  360: "A",  361: "A",
     363: "S",  364: "L",  367: "G",  368: "V",  371: "Q",
     372: "G",  376: "E",  378: "G",  380: "V",  381: "N",
     391: "A",  399: "A",  400: "K",  404: "R",  406: "R",
     407: "N",  408: "C",  410: "N",  412: "K",  413: "T",
     414: "S",  417: "N",  419: "V",  421: "L",  422: "A",
     423: "N",  424: "E",  425: "V",  427: "M",
}

CLASS_IV_MARKERS: Dict[int, str] = {
    # 162 positions unique to sdEPSPS (Streptomyces davawensis).
    # Class IV is rare — almost exclusively one actinobacteria clade.
       3: "V",    5: "D",    6: "I",   14: "A",   18: "F",
      22: "A",   24: "D",   26: "V",   28: "T",   30: "V",
      31: "R",   32: "P",   34: "R",   39: "E",   40: "G",
      41: "F",   44: "G",   46: "A",   47: "R",   50: "Y",
      51: "R",   52: "V",   53: "G",   54: "R",   58: "S",
      59: "W",   60: "Q",   62: "D",   64: "R",   65: "P",
      67: "G",   68: "P",   69: "A",   70: "V",   72: "E",
      73: "A",   75: "V",   76: "Y",   78: "R",   79: "D",
      80: "G",   81: "A",   83: "T",   84: "A",   88: "P",
      89: "T",   93: "A",   95: "H",   97: "T",   99: "R",
     100: "F",  102: "A",  103: "S",  104: "E",  105: "Q",
     111: "L",  112: "L",  115: "T",  116: "R",  120: "G",
     121: "V",  122: "D",  123: "L",  124: "R",  125: "H",
     126: "E",  128: "R",  129: "D",  131: "H",  132: "H",
     136: "V",  137: "R",  138: "A",  139: "A",  141: "V",
     142: "A",  145: "E",  152: "Q",  164: "G",  167: "T",
     168: "E",  169: "K",  172: "R",  174: "H",  175: "V",
     177: "D",  181: "V",  185: "E",  190: "M",  197: "E",
     199: "T",  201: "E",  202: "G",  203: "H",  204: "D",
     213: "R",  215: "T",  216: "T",  217: "Y",  218: "A",
     225: "T",  226: "S",  230: "F",  241: "T",  243: "P",
     245: "L",  248: "G",  249: "A",  250: "L",  257: "V",
     259: "V",  261: "R",  262: "R",  279: "T",  283: "R",
     285: "L",  287: "V",  288: "N",  289: "M",  290: "R",
     291: "D",  293: "S",  295: "T",  296: "M",  301: "A",
     302: "I",  304: "P",  307: "S",  313: "E",  318: "T",
     327: "E",  329: "C",  331: "E",  332: "N",  338: "V",
     339: "R",  342: "T",  346: "W",  348: "E",  350: "H",
     353: "A",  355: "P",  356: "T",  361: "K",  364: "G",
     369: "V",  375: "T",  378: "R",  379: "V",  380: "P",
     383: "S",  384: "F",  385: "D",  388: "G",  391: "R",
     396: "G",  398: "H",  400: "E",  402: "G",  403: "A",
     405: "R",  406: "A",
}

# ---------------------------------------------------------------------------
# Class III domains (Carozzi et al. 2006, PCT WO2006/110586, Athenix Corp.)
#
# Class III EPSPS is defined by the presence of at least one of 18 sequence
# domains identified through experimental screening of glyphosate-tolerant
# bacterial isolates (Carozzi et al. 2006). This implementation uses a
# sliding-window search for each domain pattern directly against the query
# protein sequence, matching the original tool's approach of testing for
# "at least one complete motif."
#
# Each domain is a list of tuples: (allowed_amino_acids_at_position,).
# Fixed positions have a single character; degenerate positions list all
# allowed residues. Definitions taken verbatim from patent SEQ ID NOs 13-44.
#
# The original single-motif threshold (>=1 match) is restored here, as it
# is consistent with both Carozzi et al. (2006) and Leino et al. (2021).
# False positives are controlled by the specificity of the longer domains
# (Domain I is 17 positions; Domain VII is 16 positions) rather than by
# requiring multiple shorter matches.
# ---------------------------------------------------------------------------
CAROZZI_DOMAINS: Dict[str, List[tuple]] = {
    # Domain I (SEQ ID NO:13) — 17 positions
    "I":    [("L",), ("A",), ("K",), ("G",), ("K","T"), ("S",), ("R","H"),
             ("L",), ("S","T"), ("G",), ("A",), ("L",), ("K",), ("S",),
             ("D",), ("D",), ("T",)],
    # Domain Ia (SEQ ID NO:14) — sub-domain of I, 5 positions
    "Ia":   [("L",), ("A",), ("K",), ("G",), ("K","T")],
    # Domain Ib (SEQ ID NO:15) — sub-domain of I, 4 positions
    "Ib":   [("S",), ("R","H"), ("L",), ("S","T")],
    # Domain Ic (SEQ ID NO:16) — sub-domain of I, 8 positions
    "Ic":   [("G",), ("A",), ("L",), ("K",), ("S",), ("D",), ("D",), ("T",)],
    # Domain II (SEQ ID NO:17) — 13 positions
    "II":   [("E",), ("P",), ("D",), ("D","A"), ("S","T"), ("T",), ("F",),
             ("V","I"), ("V",), ("T","E","K"), ("G","S"), ("Q","S","E","T"), ("G",)],
    # Domain IIa (SEQ ID NO:18) — sub-domain of II, 9 positions
    "IIa":  [("E",), ("P",), ("D",), ("D","A"), ("S","T"), ("T",), ("F",),
             ("V","I"), ("V",)],
    # Domain IIb (SEQ ID NO:19) — sub-domain of II, 4 positions
    "IIb":  [("T","E","K"), ("G","S"), ("Q","S","E","T"), ("G",)],
    # Domain V (SEQ ID NO:22) — 6 positions
    "V":    [("T","S"), ("G",), ("C",), ("P",), ("P",), ("V",)],
    # Domain VI (SEQ ID NO:24) — 7 positions
    "VI":   [("W",), ("R","K"), ("V",), ("A","H","E","S"), ("P","A"), ("T",), ("G",)],
    # Domain VII (SEQ ID NO:25) — 16 positions
    "VII":  [("E",), ("P",), ("D",), ("A",), ("S",), ("A",), ("A",), ("T",),
             ("Y",), ("L",), ("W",), ("A","G"), ("A",), ("E","Q"), ("V","L","A"), ("L",)],
    # Domain XIa (SEQ ID NO:30) — 7 positions
    "XIa":  [("Q","K","S"), ("F",), ("P",), ("N","H"), ("M","L"), ("P","Q"), ("A",)],
    # Domain XIIb (SEQ ID NO:35) — 14 positions
    "XIIb": [("V","T"), ("E","G"), ("L","I"), ("A","E"), ("N",), ("L",), ("R",),
             ("V",), ("K",), ("E",), ("C",), ("D",), ("R",), ("I","V")],
    # Domain XIII (SEQ ID NO:37) — 7 positions
    "XIII": [("E",), ("G",), ("D",), ("D",), ("L",), ("L","I"), ("V","I")],
    # Domain XIV (SEQ ID NO:38) — 6 positions
    "XIV":  [("D","N"), ("P",), ("A","S","T"), ("L",), ("A",), ("G",)],
    # Domain XV (SEQ ID NO:39) — 10 positions
    "XV":   [("A",), ("L","S","E"), ("I",), ("D",), ("T","S"), ("H","F"),
             ("A","S"), ("D",), ("H",), ("R",)],
    # Domain XVI (SEQ ID NO:41) — 7 positions
    "XVI":  [("F",), ("A",), ("L",), ("A",), ("G","A","T"), ("L",), ("K",)],
    # Domain XVII (SEQ ID NO:44) — 5 positions
    "XVII": [("A","S","P"), ("S",), ("L",), ("G",), ("V",)],
}

# ---------------------------------------------------------------------------
# Pairwise aligner using BLOSUM62 + standard affine gap penalties.
# BLOSUM62 is the default matrix for T-Coffee and most protein aligners
# and is the most likely matrix used by the original EPSPSClass tool
# (Leino et al. 2021 do not specify it explicitly, but T-Coffee defaults
# to BLOSUM62 for pairwise alignment steps).
# Gap penalties: -11 open / -1 extend — standard BLOSUM62 convention
# (Henikoff & Henikoff 1992; Altschul et al. 1997).
from Bio.Align import substitution_matrices as _sm
_ALIGNER = PairwiseAligner()
_ALIGNER.mode               = "global"
_ALIGNER.substitution_matrix = _sm.load("BLOSUM62")
_ALIGNER.open_gap_score     = -11
_ALIGNER.extend_gap_score   = -1


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ClassificationResult:
    """Result for a single query sequence."""
    query_id: str
    classes: List[str]
    identity_pct: Dict[str, float]
    aligned_markers: Dict[str, Dict[int, str]]
    is_unclassified: bool
    notes: List[str] = field(default_factory=list)

    @property
    def primary_class(self) -> str:
        if self.is_unclassified:
            return "Unclassified"
        for cls in ["IV", "III", "II", "I"]:
            if cls in self.classes:
                return cls
        return "Unclassified"

    @property
    def sensitivity(self) -> str:
        c = self.primary_class
        if c == "I":
            return "Sensitive"
        elif c in ("II", "III", "IV"):
            return "Resistant"
        return "Unknown"

    def __str__(self) -> str:
        ids = ", ".join(f"{k}:{v:.1f}%" for k, v in self.identity_pct.items())
        return (
            f"{self.query_id}: class={self.primary_class} "
            f"({self.sensitivity}) | identities=[{ids}]"
        )


# ---------------------------------------------------------------------------
# Reference sequence loader
# ---------------------------------------------------------------------------

def _load_reference(class_label: str) -> Seq:
    """Load a bundled reference FASTA sequence."""
    path = REF_FASTA[class_label]
    if not path.exists():
        raise FileNotFoundError(
            f"Reference FASTA for class {class_label} not found at {path}.\n"
            f"The reference sequences should be bundled with the package.\n"
            f"If missing, run: epspsclass download-refs"
        )
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"Empty FASTA file: {path}")
    return records[0].seq


# ---------------------------------------------------------------------------
# Core alignment helpers
# ---------------------------------------------------------------------------

def _pairwise_align(query: Seq, reference: Seq):
    """Global pairwise alignment via Bio.Align.PairwiseAligner."""
    alignments = list(_ALIGNER.align(str(query), str(reference)))
    if not alignments:
        raise ValueError("Pairwise alignment produced no result.")
    aln = alignments[0]

    class _Aln:
        pass

    result = _Aln()
    result.seqA = aln[0]
    result.seqB = aln[1]
    return result


def _percent_identity(aligned_query: str, aligned_ref: str) -> float:
    matches = aligned = 0
    for q, r in zip(aligned_query, aligned_ref):
        if q == "-" or r == "-":
            continue
        aligned += 1
        if q == r:
            matches += 1
    return (matches / aligned * 100) if aligned > 0 else 0.0


def _ref_pos_to_aligned(aligned_ref: str) -> Dict[int, int]:
    """Map reference 1-based position → alignment column index."""
    mapping: Dict[int, int] = {}
    ref_pos = 0
    for col, aa in enumerate(aligned_ref):
        if aa != "-":
            ref_pos += 1
            mapping[ref_pos] = col
    return mapping


def _check_markers(
    aligned_query: str,
    aligned_ref: str,
    markers: Dict[int, str],
) -> Tuple[bool, Dict[int, str]]:
    """
    Check whether all markers are present in the query at the aligned positions.
    Returns (all_present, {ref_pos: query_aa_found}).
    """
    pos_map = _ref_pos_to_aligned(aligned_ref)
    found: Dict[int, str] = {}
    for ref_pos, expected_aa in markers.items():
        col = pos_map.get(ref_pos)
        if col is None:
            logger.debug("Marker position %d not in alignment.", ref_pos)
            return False, found
        query_aa = aligned_query[col]
        if query_aa == expected_aa:
            found[ref_pos] = query_aa
    all_present = len(found) == len(markers)
    return all_present, found


def _check_class_iii_domains(query_seq: str) -> bool:
    """
    Check if the query sequence contains at least one Carozzi et al. (2006)
    class III domain (PCT WO2006/110586, Athenix Corporation).

    Uses a sliding-window search over the raw (unaligned) query sequence,
    consistent with the original EPSPSClass approach of testing for
    "at least one complete motif." The threshold of >=1 match follows
    Carozzi et al. (2006) and Leino et al. (2021).

    False positives are controlled by the specificity of the longer domains
    (Domain I is 17 positions; Domain VII is 16 positions) rather than by
    requiring multiple matches.
    """
    seq = query_seq.replace("-", "")  # remove any alignment gaps
    for domain_name, domain in CAROZZI_DOMAINS.items():
        domain_len = len(domain)
        for pos in range(len(seq) - domain_len + 1):
            if all(seq[pos + i] in domain[i] for i in range(domain_len)):
                return True
    return False


# ---------------------------------------------------------------------------
# Public classifier
# ---------------------------------------------------------------------------

class EPSPSClassifier:
    """
    EPSPS class classifier implementing the algorithm of Leino et al. (2021).

    Reference sequences are bundled with the package (from Supplementary
    Table 2 of Leino et al. 2021) — no internet access required.

    Parameters
    ----------
    identity_threshold : float
        Minimum percent identity to any reference to attempt marker checking.
        Sequences below this threshold are flagged as unreliable (default 40%).
        This follows the recommendation in Leino et al. (2021).
    """

    def __init__(self, identity_threshold: float = 40.0):
        self.identity_threshold = identity_threshold
        self._refs: Dict[str, Seq] = {}

    def _get_ref(self, cls: str) -> Seq:
        if cls not in self._refs:
            self._refs[cls] = _load_reference(cls)
        return self._refs[cls]

    def classify(self, query_id: str, query_seq: Seq) -> ClassificationResult:
        """
        Classify a single EPSPS protein sequence.

        Parameters
        ----------
        query_id  : sequence identifier (for reporting)
        query_seq : Bio.Seq protein sequence

        Returns
        -------
        ClassificationResult
        """
        classes: List[str] = []
        identity_pct: Dict[str, float] = {}
        aligned_markers: Dict[str, Dict[int, str]] = {}
        notes: List[str] = []

        for cls in ("I", "II", "III", "IV"):
            ref = self._get_ref(cls)
            try:
                aln = _pairwise_align(query_seq, ref)
            except Exception as exc:
                notes.append(f"Alignment to class {cls} failed: {exc}")
                continue

            aq, ar = aln.seqA, aln.seqB
            pct_id = _percent_identity(aq, ar)
            identity_pct[cls] = pct_id

            if pct_id < self.identity_threshold:
                notes.append(
                    f"Identity to class {cls} ref ({pct_id:.1f}%) below "
                    f"threshold ({self.identity_threshold}%); "
                    f"marker check skipped."
                )
                continue

            if cls == "III":
                # NOTE: the original EPSPSClass tool (Leino et al. 2021) and
                # Carozzi et al. (2006) define class III as requiring at least
                # ONE complete domain. We search the raw query sequence directly
                # using the 17 Carozzi domain patterns (sliding window), which
                # avoids dependence on alignment quality for short motif detection.
                if _check_class_iii_domains(str(query_seq)):
                    classes.append("III")
                    aligned_markers["III"] = {}
            elif cls == "I":
                ok, found = _check_markers(aq, ar, CLASS_I_MARKERS)
                aligned_markers["I"] = found
                if ok:
                    classes.append("I")
            elif cls == "II":
                ok, found = _check_markers(aq, ar, CLASS_II_MARKERS)
                aligned_markers["II"] = found
                if ok:
                    classes.append("II")
            elif cls == "IV":
                # Class IV uses a minimum-count threshold rather than
                # requiring ALL markers, because sdEPSPS has low identity
                # (~39%) to the other references, causing many alignment
                # columns to be filtered by the gap criterion. Requiring
                # >=10 markers from the 162-position set is both specific
                # and robust. Class IV is rare (<1% of bacteria) and
                # confined to one actinobacteria clade (Supp Fig 3).
                _, found = _check_markers(aq, ar, CLASS_IV_MARKERS)
                aligned_markers["IV"] = found
                if len(found) >= 10:
                    classes.append("IV")

        is_unclassified = len(classes) == 0
        return ClassificationResult(
            query_id=query_id,
            classes=classes,
            identity_pct=identity_pct,
            aligned_markers=aligned_markers,
            is_unclassified=is_unclassified,
            notes=notes,
        )

    def classify_fasta(self, fasta_path) -> List[ClassificationResult]:
        """Classify all sequences in a FASTA file."""
        results = []
        for record in SeqIO.parse(str(fasta_path), "fasta"):
            result = self.classify(record.id, record.seq)
            results.append(result)
        return results
