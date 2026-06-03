#!/bin/bash
# derive_markers.sh
# =================
# Run on EC2 (Amazon Linux 2023, t3.micro is sufficient — finishes in <1 min).
#
# Aligns the four EPSPS reference sequences using T-Coffee (matching the
# original EPSPSClass tool — Leino et al. 2021 Supplementary Figure 2),
# derives class-discriminating marker positions from the alignment, and
# outputs Python code ready to paste into epspsclass/classifier.py.
#
# Usage:
#   chmod +x derive_markers.sh
#   ./derive_markers.sh 2>&1 | tee marker_derivation_output.txt
#
# Output files:
#   refs/aligned_tcoffee.fasta   — T-Coffee multiple sequence alignment
#   marker_derivation_output.txt — marker positions + validation

set -euo pipefail
echo "=== EPSPSClass marker derivation (T-Coffee alignment) ==="

# ── 1. Dependencies ───────────────────────────────────────────────
echo "[1/5] Installing dependencies..."

# T-Coffee
if ! command -v t_coffee &>/dev/null; then
    echo "  Installing T-Coffee..."
    # T-Coffee static binary for Linux x86_64
    curl -fsSL \
      "https://s3.eu-west-1.amazonaws.com/tcoffee-packages/Stable/Latest/linux/T-COFFEE_installer_Version_13.46.0.919e8c6b_linux_x64.tar.gz" \
      -o /tmp/tcoffee.tar.gz
    tar -xzf /tmp/tcoffee.tar.gz -C /tmp/
    sudo cp /tmp/T-COFFEE*/bin/t_coffee /usr/local/bin/
    echo "  T-Coffee installed: $(t_coffee --version 2>&1 | head -1)"
else
    echo "  T-Coffee already available: $(t_coffee --version 2>&1 | head -1)"
fi

# Python / biopython
pip3 install biopython --break-system-packages -q
python3 -c "from Bio import SeqIO; print('  biopython ok')"

# ── 2. Write the four reference sequences ─────────────────────────
echo "[2/5] Writing reference FASTA files..."
mkdir -p refs

cat > refs/all_four_refs.fasta << 'FASTA'
>vcEPSPS_ClassI|Q9KNE7|Vibrio_cholerae
MESLTLQPIELISGEVNLPGSKSVSNRALLLAALASGTTRLTNLLDSDDIRHMLNALTKL
GVNYRLSADKTTCEVEGLGQAFHTTQPLELFLGNAGTAMRPLAAALCLGQGDYVLTGEPR
MKERPIGHLVDALRQAGAQIEYLEQENFPPLRIQGTGLQAGTVTIDGSISSQFLTAFLMS
APLAQGKVTIKIVGELVSKPYIDITLHIMEQFGVQVINHDYQEFVIPAGQSYVSPGQFLV
EGDASSASYFLAAAAIKGGEVKVTGIGKNSIQGDIQFADALEKMGAQIEWGDDYVIARRG
ELNAVDLDFNHIPDAAMTIATTALFAKGTTAIRNVYNWRVKETDRLAAMATELRKVGATV
EEGEDFIVITPPTKLIHAAIDTYDDHRMAMCFSLVALSDTPVTINDPKCTSKTFPDYFDK
FAQLSR
>cbEPSPS_ClassII|Q83EH4|Coxiella_burnetii
MDYQTIPSQGLSGEICVPGDKSISHRAVLLAAIAEGQTQVDGFLMGADNLAMVSALQQMG
ASIQVIEDENILVVEGVGMTGLQAPPEALDCGNSGTAIRLLSGLLAGQPFNTVLTGDSSL
QRRPMKRIIDPLTLMGAKIDSTGNVPPLKIYGNPRLTGIHYQLPMASAQVKSCLLLAGLY
ARGKTCITEPAPSRDHTERLLKHFHYTLQKDKQSICVSGGGKLKANDISIPGDISSAAFF
IVAATITPGSAIRLCRVGVNPTRLGVINLLKMMGADIEVTHYTEKNEEPTADITVRHARL
KGIDIPPDQVPLTIDEFPVLLIAAAVAQGKTVLRDAAELRVKETDRIAAMVDGLQKLGIA
AESLPDGVIIQGGTLEGGEVNSYDDHRIAMAFAVAGTLAKGPVRIRNCDNVKTSFPNFVE
LANEVGMNVKGVRGRGGF
>bvEPSPS_ClassIII|CAA73210|Brevundimonas_vesicularis
MMMGRAKLTIIPPGKPLTGRAMPPGSKSITNRALLLAGLAKGTSRLTGALKSDDTRYMAE
ALRAMGVTIDEPDDTTFIVKGSGKLQPPAAPLFLGNAGTATRFLTAAAALVDGKVIVDGD
AHMRKRPIGPLVDALRSLGIDASAETGCPPVTINGTGRFEASRVQIDGGLSSQYVSALLM
MAAGGDRAVDVELLGEHIGALGYIDLTVAAMRAFGAKVERVSPVAWRVEPTGYHAADFVI
EPDASAATYLWAAEVLSGGKIDLGTPAEQFSQPDAKAYDLISKFPHLPAVIDGSQMQDAI
PTLAVLAAFNEMPVRFVGIENLRVKECDRIRALSSGLSRIVPNLGTEEGDDLIIASDPSL
AGKILTAEIDSFADHRIAMSFALAGLKIGGITILDPDCVAKTFPSYWNVLSSLGVAYED
>sdEPSPS_ClassIV|H6WNZ5|Streptomyces_davawensis
MPVADIPGSKSITARALFLAAAADGVTTLVRPLRSDDTEGFAEGLARLGYRVGRTPDSWQ
VDGRPQGPAVAEADVYCRDGATTARFLPTLAAAGHGTYRFDASEQMRRRPLLPLTRALRE
LGVDLRHEERDGHHPLTVRAAGVAGGEVTLDAGQSSQYLTALLLLGPLTEKGLRIHVTDL
VSVPYIEITLAMMRAFGVEVTREGHDFVVPPGGYRATTYAIEPDASTSSYFFAAAALSGG
EVTVPGLGEGALQGDLGFVDVLRRMGAEVEIGADRTTVRGTGELRGLTVNMRDISDTMPT
LAAIAPFASGPVRIEDVANTRVKECDRLEACAENLRRLGVRVETGPDWIEIHPGATPTGA
EIKTYGDHRIVMSFAVTGLRVPGISFDDPGCVRKTFPGFHEEFGALRARL
FASTA

echo "  Reference FASTA written: refs/all_four_refs.fasta"

# ── 3. T-Coffee alignment ─────────────────────────────────────────
echo "[3/5] Running T-Coffee alignment (matching Leino et al. 2021 Supp Fig 2)..."
t_coffee refs/all_four_refs.fasta \
    -output fasta_aln \
    -outfile refs/aligned_tcoffee.fasta \
    -quiet
echo "  Alignment written: refs/aligned_tcoffee.fasta"

# Print alignment summary
python3 -c "
from Bio import SeqIO
records = list(SeqIO.parse('refs/aligned_tcoffee.fasta', 'fasta'))
print(f'  Sequences: {len(records)}')
for r in records:
    cls = [c for c in [\"I\",\"II\",\"III\",\"IV\"] if c in r.id.split(\"|\")[0]]
    print(f'    Class {cls[0] if cls else \"?\"}: {r.id[:30]}, aligned length={len(r.seq)}')
"

# ── 4. Derive markers ─────────────────────────────────────────────
echo "[4/5] Deriving class-discriminating marker positions..."

python3 << 'PYEOF'
"""
Derives EPSPS class marker positions from the T-Coffee multiple sequence
alignment of the four reference sequences.

Methodology (matching Leino et al. 2021):
- Each class's markers are alignment columns where that class has a unique,
  non-gap amino acid not shared by any other class.
- Class III motifs are consecutive triplets of such exclusive positions in
  bvEPSPS, matching Carozzi et al. (2006)'s motif definition.
- All positions are reported in REFERENCE SEQUENCE coordinates (1-based),
  as used by the classification algorithm at runtime.

The key canonical marker (Pro101 in vcEPSPS / Leu100 in cbEPSPS) must
appear in both CLASS_I_MARKERS and CLASS_II_MARKERS — this is the primary
validation check. If it does not appear, the alignment differs from
Supplementary Table 1 and should be inspected manually.
"""

from Bio import SeqIO
from typing import Dict, List, Tuple

ALN_FILE = "refs/aligned_tcoffee.fasta"

records = list(SeqIO.parse(ALN_FILE, "fasta"))
assert len(records) == 4, f"Expected 4 sequences, got {len(records)}"

# Map to class labels
def get_class(rec_id):
    for cls in ["ClassI", "ClassII", "ClassIII", "ClassIV"]:
        if cls in rec_id:
            return cls.replace("Class", "")
    raise ValueError(f"Cannot determine class from ID: {rec_id}")

aln = {get_class(r.id): str(r.seq) for r in records}
aln_len = len(next(iter(aln.values())))
assert all(len(s) == aln_len for s in aln.values()), "Alignment length mismatch"

def col_to_seqpos(aligned_seq: str) -> Dict[int, int]:
    """Map alignment column → reference sequence position (1-based)."""
    m = {}
    pos = 0
    for col, aa in enumerate(aligned_seq):
        if aa != "-":
            pos += 1
            m[col] = pos
    return m

pos_maps = {cls: col_to_seqpos(seq) for cls, seq in aln.items()}

def find_unique_markers(target_cls: str) -> Dict[int, str]:
    """
    Find positions where target_cls has a unique non-gap residue
    not shared by any other class. Returns {seqpos: aa}.
    """
    other = [c for c in ["I", "II", "III", "IV"] if c != target_cls]
    markers = {}
    for col in range(aln_len):
        taa = aln[target_cls][col]
        if taa == "-":
            continue
        others = [aln[c][col] for c in other]
        if any(aa == "-" for aa in others):
            continue  # skip gapped columns
        if taa not in others:
            seqpos = pos_maps[target_cls].get(col)
            if seqpos:
                markers[seqpos] = taa
    return markers

def find_class_iii_motifs() -> Tuple[List, int]:
    """
    Find consecutive exclusive triplets in bvEPSPS (class III).
    Returns (motif_list, n_exclusive_cols).
    """
    other = ["I", "II", "IV"]
    bv = aln["III"]
    bv_map = pos_maps["III"]

    # All alignment columns where bvEPSPS has a unique residue
    excl = []
    for col in range(aln_len):
        baa = bv[col]
        if baa == "-":
            continue
        others = [aln[c][col] for c in other]
        if any(aa == "-" for aa in others):
            continue
        if baa not in others:
            excl.append(col)

    # Build (col, seqpos, aa) list for exclusive columns
    pts = [(col, bv_map[col], bv[col]) for col in excl if col in bv_map]

    # Find non-overlapping consecutive triplets (adjacent in bvEPSPS sequence)
    motifs = []
    i = 0
    while i <= len(pts) - 3:
        c1, p1, a1 = pts[i]
        c2, p2, a2 = pts[i+1]
        c3, p3, a3 = pts[i+2]
        if p2 == p1 + 1 and p3 == p1 + 2:
            motifs.append((p1, a1, p2, a2, p3, a3))
            i += 3
        else:
            i += 1

    return motifs, len(excl)

# ── Run derivation ─────────────────────────────────────────────────
print("\n" + "=" * 65)
print("DERIVED MARKERS — paste into epspsclass/classifier.py")
print("=" * 65)

for cls, label in [("I","vcEPSPS"), ("II","cbEPSPS"), ("IV","sdEPSPS")]:
    markers = find_unique_markers(cls)
    print(f"\n# Class {cls} markers ({label} coordinates)")
    print(f"# {len(markers)} discriminating positions found")
    print(f"CLASS_{cls}_MARKERS: Dict[int, str] = {{")
    for pos, aa in sorted(markers.items()):
        print(f"    {pos:4d}: \"{aa}\",")
    print("}")

motifs, n_excl = find_class_iii_motifs()
print(f"\n# Class III motifs (bvEPSPS coordinates)")
print(f"# {n_excl} exclusive columns; {len(motifs)} consecutive exclusive triplets")
print("CLASS_III_MOTIFS: List[Tuple[int,str,int,str,int,str]] = [")
for p1,a1,p2,a2,p3,a3 in motifs:
    print(f"    ({p1:3d}, \"{a1}\", {p2:3d}, \"{a2}\", {p3:3d}, \"{a3}\"),")
print("]")

# ── Validation ────────────────────────────────────────────────────
print("\n" + "=" * 65)
print("VALIDATION")
print("=" * 65)

# 1. Canonical Pro/Leu switch (Supp Table 1 row 106)
print("\n[1] Canonical Pro101(I)/Leu100(II) sensitivity marker:")
found_canonical = False
for col in range(aln_len):
    if aln["I"][col] == "P" and aln["II"][col] == "L":
        vc_pos = pos_maps["I"].get(col)
        cb_pos = pos_maps["II"].get(col)
        bv_aa  = aln["III"][col]
        sd_aa  = aln["IV"][col]
        print(f"  aln col {col}: vcEPSPS[{vc_pos}]=P  cbEPSPS[{cb_pos}]=L  "
              f"bvEPSPS={bv_aa}  sdEPSPS={sd_aa}")
        if vc_pos == 101:
            print("  ✓ vcEPSPS position 101 confirmed (matches Supp Table 1 row 106)")
            found_canonical = True
if not found_canonical:
    print("  ✗ WARNING: canonical position not found at vcEPSPS 101")
    print("    Inspect refs/aligned_tcoffee.fasta vs Supplementary Table 1")

# 2. Cross-check marker counts against Supp Table 1
# (Supp Table 1 has 464 alignment rows — some are gaps)
print(f"\n[2] Alignment length: {aln_len} columns")
print(f"    (Supplementary Table 1 has 464 rows — should be close)")

# 3. Self-classification check
print("\n[3] Self-classification (reference should classify as its own class):")
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.substitution_matrix = __import__('Bio.Align.substitution_matrices',
    fromlist=['load']).load("BLOSUM62")
aligner.open_gap_score   = -11
aligner.extend_gap_score = -1

def pct_id(aq, ar):
    m = a = 0
    for q, r in zip(aq, ar):
        if q == "-" or r == "-": continue
        a += 1
        if q == r: m += 1
    return (m/a*100) if a else 0

ref_seqs_raw = {}
for r in SeqIO.parse("refs/all_four_refs.fasta", "fasta"):
    cls = get_class(r.id)
    ref_seqs_raw[cls] = str(r.seq)

for query_cls in ["I","II","III","IV"]:
    q = ref_seqs_raw[query_cls]
    best_cls, best_pct = None, 0
    for ref_cls in ["I","II","III","IV"]:
        alns = list(aligner.align(q, ref_seqs_raw[ref_cls]))
        pct = pct_id(alns[0][0], alns[0][1])
        if pct > best_pct:
            best_pct, best_cls = pct, ref_cls
    ok = "✓" if best_cls == query_cls else "✗ FAIL"
    print(f"  Class {query_cls} ref: best match = Class {best_cls} ({best_pct:.1f}%)  {ok}")

print("\nDone. Copy the marker blocks above into epspsclass/classifier.py.")
PYEOF

# ── 5. Save alignment ─────────────────────────────────────────────
echo "[5/5] Saving files..."
echo "  refs/aligned_tcoffee.fasta — T-Coffee alignment"
echo "  Run this script with: ./derive_markers.sh 2>&1 | tee marker_derivation_output.txt"
echo ""
echo "=== Next steps ==="
echo "1. Copy the CLASS_I_MARKERS, CLASS_II_MARKERS, CLASS_IV_MARKERS, and"
echo "   CLASS_III_MOTIFS blocks from output into epspsclass/classifier.py"
echo "2. Update the aligner in classifier.py to use BLOSUM62 (see note below)"
echo "3. Run: python -m pytest tests/ -v"
echo "4. Commit and push to GitHub"
echo ""
echo "Note on BLOSUM62: also update classifier.py _ALIGNER to use:"
echo "  from Bio.Align import substitution_matrices"
echo "  _ALIGNER.substitution_matrix = substitution_matrices.load('BLOSUM62')"
echo "  _ALIGNER.open_gap_score   = -11  # standard BLOSUM62 gap penalties"
echo "  _ALIGNER.extend_gap_score = -1"
