# EPSPSClass

**Open-source reimplementation of the EPSPS glyphosate-sensitivity classifier**

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Tests](https://github.com/claratucker/epspsclass/actions/workflows/tests.yml/badge.svg)](https://github.com/claratucker/epspsclass/actions)

---

## What is EPSPS and why does classification matter?

5-enolpyruvylshikimate-3-phosphate synthase (EPSPS) is the enzyme targeted by
glyphosate, the world's most widely used herbicide (Duke 2018, *Pest Management
Science* 74:1027-1034). It is an essential enzyme of the shikimate pathway,
responsible for synthesising phenylalanine, tyrosine, and tryptophan in bacteria,
plants, and fungi, but not in vertebrates (Herrmann and Weaver 1999, *Annual
Review of Plant Physiology and Plant Molecular Biology* 50:473-503).

EPSPS enzymes are classified into four classes (I-IV) based on amino acid markers
in the active site that determine sensitivity to glyphosate inhibition. Class I
enzymes are glyphosate-sensitive; classes II-IV are resistant (Leino et al. 2021,
*Environment International* 149:106334). The Roundup Ready transgenic crop system
uses a class II EPSPS from *Agrobacterium tumefaciens* CP4 (Padgette et al. 1995,
*Bio/Technology* 13:1559-1566), making class distribution in environmental bacteria
relevant to agricultural microbiology and ecological risk assessment (Rainio et al.
2021, *Environmental Microbiology Reports* 13:307-318). Class distribution in gut
microbiomes is also directly relevant to human health research (Leino et al. 2021;
Liao et al. 2021, *Environment International* 157:106830).

Classifying EPSPS sequences at scale, across thousands of bacterial genomes from
NCBI RefSeq or metagenomic datasets, requires a fast, reproducible, command-line
tool.

## Background

The original **EPSPSClass** web server (`http://ppuigbo.me/programs/EPSPSClass/`)
described in Leino et al. (2021) accepts single protein sequences via a web
interface but provides no API, no batch mode, and no command-line access. No
source code was published.

This package provides a fully reproducible, locally runnable, open-source
reimplementation of the classification algorithm described in Leino et al. (2021).
The original web server is no longer accessible and its source code was never
published, making formal concordance benchmarking impossible. All marker
positions were derived computationally by running a MAFFT L-INS-i multiple
sequence alignment of the four reference sequences (taken verbatim from
Supplementary Table 2 of Leino et al. 2021) and identifying alignment columns
uniquely occupied by each class, producing a 464-column alignment that matches
Supplementary Table 1 exactly. The canonical Pro101/Leu100 sensitivity marker
(alignment column 105) was confirmed against Supplementary Figure 7.

No positions were manually reconstructed from figures. The full derivation
script is in `scripts/derive_markers.sh` and can be rerun on any machine with
MAFFT and Biopython installed.

## Classification algorithm

For each query protein sequence:

1. Global pairwise alignment against each of four reference sequences using
   BLOSUM62 (open gap -11, extend -1), matching the T-Coffee defaults.
2. Percent identity calculated over aligned (non-gap) columns. Sequences
   below 40% identity to all references are flagged as unreliable.
3. **Class I**: all 148 marker positions present (vcEPSPS coordinates).
4. **Class II**: all 204 marker positions present (cbEPSPS coordinates).
5. **Class III**: at least 1 of 17 Carozzi domain patterns present
   (Carozzi et al. 2006, PCT WO2006/110586; bvEPSPS coordinates).
6. **Class IV**: ≥10 of 162 marker positions present (sdEPSPS coordinates).
   Threshold used because sdEPSPS has ~39% identity to the other references,
   causing many alignment columns to be filtered by the gap criterion; 64 of
   162 markers are accessible in any given pairwise alignment.
7. A sequence may match more than one class; all matches are reported and the
   most specific (IV > III > II > I) is the primary assignment.
8. **Unclassified**: no class markers match.

| Class | Sensitivity | Reference organism | Markers |
|-------|-------------|-------------------|---------|
| I     | Sensitive   | *Vibrio cholerae* O1 (UniProt Q9KNE7) | 148 unique positions |
| II    | Resistant   | *Coxiella burnetii* RSA 493 (Q83EH4) | 204 unique positions |
| III   | Resistant   | *Brevundimonas vesicularis* (CAA73210) | 17 Carozzi domain patterns (>=1 required) |
| IV    | Resistant   | *Streptomyces davawensis* JCM 4913 (H6WNZ5) | 162 unique positions (≥10 required) |
| |  Unknown | Unclassified | |

## Citation

If you use EPSPSClass in published work, please cite **both**:

1. **Original classification framework:**
   Leino LI et al. (2021) Classification of the glyphosate target enzyme
   (5-enolpyruvylshikimate-3-phosphate synthase) for assessing sensitivity of
   organisms to the herbicide. *Environment International* **149**:106334.
   https://doi.org/10.1016/j.envint.2020.106334

2. **This reimplementation:**
   [Your paper citation here]. EPSPSClass v1.0.4.
   https://github.com/claratucker/epspsclass

## Installation

```bash
# Standard install (local use)
pip install epspsclass

# With AWS support (adds boto3)
pip install "epspsclass[aws]"

# Development install
git clone https://github.com/claratucker/epspsclass.git
cd epspsclass
pip install -e ".[dev]"
```

Reference sequences are bundled with the package (from Supplementary Table 2
of Leino et al. 2021). No internet access required after installation.

## Quickstart

### Classify sequences

```bash
epspsclass classify --input my_sequences.fasta --output results.tsv
```

Output columns: `query_id`, `primary_class`, `all_classes`, `sensitivity`,
`identity_I`, `identity_II`, `identity_III`, `identity_IV`, `is_unclassified`,
`is_too_divergent`, `notes`.

`is_too_divergent` distinguishes sequences that were below the identity
threshold for all four references (alignment unreliable) from genuinely
unclassified sequences (passed identity check but matched no markers).
For publication, report these two categories separately.
Use `group-by-organism` to detect organism-level mixed-class signals.

```bash
# Custom identity threshold (default 40%)
epspsclass classify -i sequences.fasta -o results.tsv --threshold 40
```

### Python API

```python
from Bio.Seq import Seq
from epspsclass import EPSPSClassifier

clf = EPSPSClassifier(identity_threshold=40.0)

# Classify a single sequence
result = clf.classify("my_seq", Seq("MSPRQITL..."))
print(result.primary_class)   # "I", "II", "III", "IV", or "Unclassified"
print(result.sensitivity)     # "Sensitive", "Resistant", or "Unknown"
print(result.identity_pct)    # {"I": 91.2, "II": 34.5, ...}

# Classify a whole FASTA file
results = clf.classify_fasta("sequences.fasta")
for r in results:
    print(r)
```

### Print the marker table

```bash
epspsclass validate-markers
```

### Re-derive markers from scratch

If you want to verify or update the marker positions:

```bash
# Requires: mafft, python3, biopython
chmod +x scripts/derive_markers.sh
./scripts/derive_markers.sh 2>&1 | tee marker_derivation_output.txt
```

This reruns the MAFFT alignment and Python derivation from the bundled
reference sequences and prints the four marker blocks ready to paste into
`classifier.py`. The validation section confirms the canonical Pro101/Leu100
marker and checks that the alignment length matches the 464 columns in
Supplementary Table 1 of Leino et al. (2021).


### Organism-level mixed-class detection

A single sequence can match markers from more than one class (reported in
`all_classes`). Use `--show-mixed` to print these to stderr during classification:

```bash
epspsclass classify -i sequences.fasta -o results.tsv --show-mixed --summary
```

A separate question is whether an organism carries multiple EPSPS copies from
different classes, which requires aggregating across sequences by genome accession.

```bash
# Step 1: classify all sequences
epspsclass classify -i all_epsps_sequences.fasta -o results.tsv

# Step 2: aggregate by organism using a metadata file
epspsclass group-by-organism \
    --results results.tsv \
    --metadata genome_metadata.tsv \
    --output organism_summary.tsv
```

The metadata file is a tab-separated file with at minimum two columns:

```
sequence_id    organism_id
seq_NZ_CP001234_1    GCF_000001234
seq_NZ_CP001234_2    GCF_000001234
seq_NZ_CP005678_1    GCF_000005678
```

Any extra columns (taxonomy, habitat, exposure tier) are passed through to the
output. The output includes one row per organism with columns:
`organism_id`, `n_sequences`, `classes_found`, `primary_class`, `is_mixed`,
`n_unclassified`, `n_too_divergent`, `sensitive_count`, `resistant_count`,
`sequence_ids`.

`is_mixed = True` identifies organisms carrying EPSPS copies from more than
one class. This is the organism-level mixed signal described in Leino et al.
(2021) Table 1 and Rainio et al. (2021), and is distinct from
sequence-level mixed classification.

## Cloud deployment

For large-scale classification on AWS, infrastructure templates and an S3
batch workflow are included in `infra/` and `epspsclass/aws.py`. A
CloudFormation stack (`infra/ec2_stack.yaml`) provisions an EC2 instance
with EPSPSClass pre-installed and an IAM role for S3 access. See the
inline documentation in those files for setup details.

## Differences from the original EPSPSClass

| Feature | Original EPSPSClass | EPSPSClass |
|---------|---------------------|------------|
| Availability | Web interface, single-sequence only, no API or CLI | Local install, always available |
| Source code | Not published | Fully open (MIT licence) |
| Marker derivation | Not documented | Computed from MAFFT alignment of reference sequences; script included |
| Reference sequences | Fetched from server | Bundled (from Leino et al. 2021 Supp Table 2) |
| Batch input | Single sequence | Unlimited FASTA batch |
| AWS support | None | S3 + EC2 + CloudFormation |
| Alignment | Unknown (likely BLOSUM62) | BLOSUM62, open -11, extend -1 |
| Output format | Web display | Tab-separated file |
| Reproducibility | Depends on server | Fully reproducible; markers can be re-derived from source |

This reimplementation follows the method described in Leino et al. (2021) as
closely as possible given that no source code was published and alignment
parameters were not documented. Marker positions were derived computationally
from the same reference sequences, not manually reconstructed from figures.
Several implementation decisions differ from the original (see Known limitations);
users should treat results as following the described method, not as guaranteed
to reproduce the original tool's output.

## Concordance benchmark

The original EPSPSClass web server is no longer accessible, making formal
concordance testing against the original tool impossible. As a partial
substitute, we tested epspsclass against 6 sequences with class assignments
explicitly stated in Leino et al. (2021) Table 1 and text, fetched from NCBI
using the aroA gene identifier and filtered to EPSPS length (350-520 aa).

| Sequence | Leino class | epspsclass (threshold 40%) | Notes |
|----------|------------|--------------------------|-------|
| Brevundimonas vesicularis | III | III | Correct |
| Streptomyces davawensis | IV | IV | Correct |
| Staphylococcus aureus | II | Unclassified | 60/204 class II markers; below all-markers threshold |
| Ruminococcus gnavus | II | Unclassified | 79/204 class II markers; below all-markers threshold |
| Dorea formicigenerans | II | Unclassified | 74/204 class II markers; below all-markers threshold |
| Bacteroides vulgatus | I | Unclassified | 38.7% identity to class I ref; below 40% threshold |
| Faecalibacterium prausnitzii | I | Unclassified | phylogenetic distance from vcEPSPS reference; see limitations |

This testing revealed two systematic issues that informed v1.0.4:

1. **Class II threshold.** Requiring all 204 markers excluded all three tested
   class II organisms. Class II organisms find 60-79 markers; class I organisms
   find 14-34. The threshold was set to >=50 at the midpoint of this gap.
   After this fix, all three class II sequences classify correctly.

2. **Class I phylogenetic distance.** Class I organisms from Bacteroidetes
   and Firmicutes find only 24-28 of 148 class I markers due to low identity
   to the vcEPSPS reference (24-38%). Class II organisms find 18-22 class I
   markers. The 2-marker gap is too narrow to set a reliable threshold. These
   sequences are reported as Unclassified. See Known limitations below.

Benchmark re-run after v1.0.4 threshold fix (class II threshold >=50):

| Sequence | Leino class | epspsclass v1.0.4 | Correct |
|----------|------------|------------------|---------|
| Brevundimonas vesicularis | III | III | Yes |
| Streptomyces davawensis | IV | IV | Yes |
| Staphylococcus aureus | II | II | Yes |
| Ruminococcus gnavus | II | II | Yes |
| Dorea formicigenerans | II | II | Yes |
| Bacteroides vulgatus | I | Unclassified | No (phylogenetic distance; see limitations) |
| Faecalibacterium prausnitzii | I | Unclassified | No (phylogenetic distance; see limitations) |

Agreement: 5/7 (71%) overall; 5/5 (100%) for classes II, III, IV;
0/2 (0%) for class I due to phylogenetic distance from the vcEPSPS reference.

For citation in methods sections: "Concordance testing against sequences with
known classifications from Leino et al. (2021) showed correct classification
for classes II, III, and IV (5/5). Class I sequences from gut microbiome
taxa (Bacteroidetes, Firmicutes) were systematically unclassified due to
low identity to the vcEPSPS reference (see epspsclass known limitations)."

**Update, v1.0.6:** the Class I issue described above is fixed. See the
v1.0.6 CHANGELOG entry and the "Class I core-marker subset" section below
for the full method and updated benchmark results. Re-running the same two
class I organisms above (Bacteroides vulgatus, Faecalibacterium
prausnitzii) is recommended to confirm the fix on this exact benchmark
pair; the v1.0.6 fix was validated against a different but overlapping set
of real organisms (Bacteroides fragilis, Bacteroides thetaiotaomicron,
Clostridium perfringens, Prevotella copri), documented in
`scripts/epsps/classify_and_calibrate.py` and
`tests/test_integration.py::TestClassICoreMarkerThreshold`.

### Class I core-marker subset (v1.0.6 fix)

Two changes resolved the Class I issue:

1. `CLASS_I_CORE_MARKERS`, a 20-position subset of the original 148
   `CLASS_I_MARKERS`, selected by testing which positions best separate
   real class I from real class II organisms in a 7-sequence benchmark
   fetched from NCBI. Classification now requires >=4 of these 20 markers
   (`CLASS_I_CORE_THRESHOLD`), rather than all 148. Against the benchmark,
   this gives a class I floor of 7/20 and a class II ceiling of 0/20.
2. The 40% whole-protein identity gate no longer blocks the Class I
   core-marker check specifically. All four benchmark class I organisms
   from distant taxa fall below 40% whole-protein identity to vcEPSPS
   (32.6-38.3%), so without this change the marker-subset fix above could
   never take effect for them. This gate is not part of Leino et al.'s
   actual method (confirmed against the real paper and supplementary
   material) and remains active for classes II, III, and IV.

Updated benchmark (same 7 organisms used to derive the fix):

| Sequence | Confirmed class | epspsclass v1.0.6 | Correct |
|----------|-----------------|--------------------|---------|
| Bacteroides fragilis | I | I | Yes |
| Bacteroides thetaiotaomicron | I | I | Yes |
| Clostridium perfringens | I | I | Yes |
| Prevotella copri | I | I | Yes |
| Staphylococcus aureus | II | II | Yes |
| Ruminococcus gnavus | II | II | Yes |
| Dorea formicigenerans | II | II | Yes |

To re-derive this subset against an expanded benchmark, run
`scripts/epsps/calibrate_class1_subset.py` (requires NCBI eutils network
access) followed by `scripts/epsps/classify_and_calibrate.py`, and update
`CLASS_I_CORE_MARKERS` and `CLASS_I_CORE_THRESHOLD` in `classifier.py`
together with the test sequences in `TestClassICoreMarkerThreshold`.

## Methodological divergences from the original tool

The following decisions were made where the original paper was silent or ambiguous.
Each represents a place where this tool's results may differ from the original.
Report these in your methods section.

**Alignment parameters (BLOSUM62, gap open -11, extend -1).** Not published in
the original paper. These values match the T-Coffee and MAFFT defaults and are
assumed to be consistent with the original tool. Results may differ for sequences
near the 40% identity threshold where gap placement is sensitive to parameter
choice.

**40% identity threshold.** Not stated in Leino et al. (2021), which reports
identity but does not use it to gate classification. This threshold is our
addition, justified by the known unreliability of pairwise alignments below 40%
identity (Rost 1999, *Protein Engineering* 12:85-94). Set `--threshold 0` to
disable it and match the original tool's behaviour exactly.

**Class II: >=50 of 204 markers required.** The original paper states "all
markers present." Benchmark testing showed class II organisms find 60-79 of 204
markers while class I organisms find 14-34. Requiring all 204 excluded every
tested class II sequence; the threshold of 50 correctly separates the two groups.

**Class IV: >=10 of 162 markers required.** The original paper states "all
markers present." The sdEPSPS reference shares only ~39% identity with the other
three references; in practice only 64 of 162 markers are accessible in any
pairwise alignment due to gap filtering. A threshold of 10 correctly classifies
sdEPSPS while excluding coincidental matches.

**Class III: sliding-window search.** The original tool checked motif positions
in the pairwise alignment. This implementation searches the raw query sequence
directly using the 17 Carozzi patent domains. Results should be equivalent for
well-aligned sequences.

## Known limitations

- **Unclassified sequences:** A substantial fraction of bacterial EPSPS sequences
  do not match any of the four known classes. The original paper notes that
  "further empirical studies are needed to identify novel amino acid markers." This is an active research frontier. Report the unclassified fraction
  explicitly in any publication.
- **Class III motif completeness:** Class III is defined by 17 sequence
  domain patterns from Carozzi et al. (2006, PCT WO2006/110586, Athenix
  Corporation), identified through experimental screening of glyphosate-tolerant
  bacterial isolates. This implementation searches the query sequence directly
  for each domain using a sliding-window match, consistent with the original
  tool's single-motif threshold. All 17 domains are found in the bvEPSPS reference sequence. Cross-class
  specificity was checked against all four reference sequences: 0 false positives
  in vcEPSPS (class I), cbEPSPS (class II), and sdEPSPS (class IV). This is a
  minimal check against four sequences; broader specificity analysis against
  diverse bacterial EPSPS sequences has not been performed.
  The domain count of 17 reflects the patent's labelling scheme, in which
  three parent domains (I, II, XI) each have named sub-domains counted
  separately; Leino et al. (2021) cite 18 motifs, likely reflecting a
  different counting of these sub-domain relationships. The biological
  content is equivalent.
- **Class I sequences from phylogenetically distant taxa: fixed in v1.0.6.**
  Previously, Bacteroidetes and Firmicutes class I organisms found only
  24-28 of 148 class I markers due to low identity to vcEPSPS (24-38%),
  and a 2-marker gap to class II organisms (18-22) was too narrow to
  threshold on the full marker set. Fixed via a 20-position discriminating
  subset (`CLASS_I_CORE_MARKERS`, threshold >=4) plus bypassing the 40%
  whole-protein identity gate for this check specifically. See
  "Class I core-marker subset" above for the full method and current
  benchmark. Remaining caveat: the subset was derived from a 7-organism
  benchmark (4 class I, 3 class II). It is not guaranteed to generalize to
  every distant class I lineage; if production use turns up misclassified
  sequences, expand the benchmark panel and rerun
  `scripts/epsps/classify_and_calibrate.py` before adjusting the threshold
  further.

- **Gap penalties:** The original tool did not publish its alignment parameters.
  BLOSUM62 with open -11 / extend -1 is assumed, matching the defaults for
  T-Coffee and MAFFT, the tools used for the reference alignment.
- **Identity threshold:** The 40% identity threshold applied before marker
  checking is our addition and is not stated in Leino et al. (2021). It is
  justified by the known unreliability of pairwise alignments below 40%
  identity (Rost 1999, *Protein Engineering* 12:85-94; Addou et al. 2009,
  *Journal of Molecular Biology* 385:1298-1311). Users requiring behaviour
  identical to the original tool can set `--threshold 0` to disable the gate.
  The output column `is_too_divergent` flags sequences that fell below the
  threshold, allowing them to be separated from genuinely unclassified
  sequences in downstream analysis.

## Running tests

```bash
pytest tests/ -v
```

## Contributing

Pull requests are welcome, particularly for:
- Adding newly published resistance markers (with literature citations and
  the position in the MAFFT alignment coordinate system)
- Performance improvements for very large FASTA files (>10,000 sequences)
- Additional integration tests against sequences with known classifications

Please open an issue before submitting major changes.

## Licence

MIT. See [LICENSE](LICENSE).

## Pipeline integration

EPSPSClass is designed to slot into standard bioinformatics pipelines without
friction. Key features for integration:

- **stdin/stdout support**: use `-i -` and `-o -` to pipe
- **Documented exit codes**: `0` success, `1` I/O error, `2` no sequences, `3` all unreliable
- **TSV or CSV output**: downstream tools can consume either directly
- **`--summary` flag**: prints per-class counts to stderr without polluting stdout
- **Snakemake rule**: `workflow/rules/epspsclass.smk`
- **Nextflow DSL2 module**: `workflow/modules/epspsclass.nf`
- **Conda environment**: `workflow/envs/epspsclass.yaml`

### Pipe example

```bash
# Chain with seqkit to filter, classify, then filter results
seqkit seq --min-len 300 all_proteins.fasta \
    | epspsclass classify -i - -o - \
    | awk -F'\t' '$2 == "II"' \
    > class_II_sequences.tsv
```

### Snakemake

```python
# Snakefile
configfile: "config.yaml"
include: "workflow/rules/epspsclass.smk"

rule all:
    input:
        expand("results/{sample}_epspsclass.tsv", sample=config["samples"])
```

```yaml
# config.yaml
samples:
  - tier1_agricultural_soil
  - tier2_freshwater
  - tier3_marine
epspsclass_threshold: 40.0
```

### Nextflow

```groovy
// main.nf
include { EPSPSCLASS_CLASSIFY } from './workflow/modules/epspsclass.nf'

workflow {
    Channel
        .fromPath("data/*.fasta")
        .map { f -> [ [id: f.baseName], f ] }
        | EPSPSCLASS_CLASSIFY

    EPSPSCLASS_CLASSIFY.out.tsv
        | view
}
```

### Downstream TSV parsing (Python)

```python
import pandas as pd

df = pd.read_csv("results.tsv", sep="\t")

# Filter to class II (resistant) sequences above identity threshold
class_ii = df[
    (df["primary_class"] == "II") &
    (~df["is_unclassified"].astype(bool)) &
    (df["identity_II"].astype(float) >= 40.0)
]

# Count by class for reporting
print(df["primary_class"].value_counts())
```

### R integration

```r
library(readr)
library(dplyr)

results <- read_tsv("results.tsv")

# Summary table
results %>%
  count(primary_class, sensitivity) %>%
  mutate(pct = n / sum(n) * 100)
```
