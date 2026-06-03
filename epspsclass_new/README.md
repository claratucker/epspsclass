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
described in Leino et al. (2021) is no longer reliably available, having been
hosted on a personal academic domain that is no longer maintained. No source
code was published alongside the tool.

This package provides a fully reproducible, locally runnable, open-source
reimplementation of the identical classification algorithm. All marker
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
   (bvEPSPS coordinates). Threshold of 2 rather than 1 reduces false positives
   from conserved tripeptides.
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
| III   | Resistant   | *Brevundimonas vesicularis* (CAA73210) | 21 exclusive triplet motifs |
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
   [Your paper citation here]. EPSPSClass v1.0.3.
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
`all_classes`). A separate question is whether an organism carries multiple
EPSPS copies from different classes, which requires aggregating across sequences
by genome accession.

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

## AWS Deployment

For large-scale analysis (thousands of bacterial genomes), we provide AWS
infrastructure support.

### EC2 instance setup

Provision a reproducible compute environment using the included CloudFormation
template. A **t3.micro** is sufficient for classification alone; use larger
instances for the full EPSPS evolutionary genomics pipeline.

```bash
aws cloudformation deploy \
  --template-file infra/ec2_stack.yaml \
  --stack-name epspsclass \
  --parameter-overrides \
      InstanceType=c6i.2xlarge \
      S3BucketName=my-epsps-bucket \
      KeyName=my-keypair \
  --capabilities CAPABILITY_NAMED_IAM
```

This creates an EC2 instance with EPSPSClass pre-installed and an IAM role
granting the instance read/write access to your S3 bucket only (no
hard-coded credentials).

Connect via SSM Session Manager (no SSH key or open port 22 required):

```bash
aws ssm start-session --target <instance-id>
```

### S3 batch workflow

```python
from epspsclass.aws import run_batch_from_s3

run_batch_from_s3(
    input_s3_uri  = "s3://my-epsps-bucket/input/epsps_sequences.fasta",
    output_s3_uri = "s3://my-epsps-bucket/results/classification.tsv",
)
```

### Recommended instance sizes

| Sequences | Recommended instance | Approx. runtime |
|-----------|---------------------|----------------|
| Marker derivation only | t3.micro | < 1 min |
| ≤ 500     | t3.medium           | < 5 min        |
| 500–2,000 | c6i.xlarge          | ~15 min        |
| 2,000–5,000 | c6i.2xlarge       | ~45 min        |
| 5,000+    | c6i.4xlarge         | ~90 min        |

## Differences from the original EPSPSClass

| Feature | Original EPSPSClass | EPSPSClass |
|---------|---------------------|------------|
| Availability | Web server (currently down) | Local install, always available |
| Source code | Not published | Fully open (MIT licence) |
| Marker derivation | Not documented | Computed from MAFFT alignment of reference sequences; script included |
| Reference sequences | Fetched from server | Bundled (from Leino et al. 2021 Supp Table 2) |
| Batch input | Single sequence | Unlimited FASTA batch |
| AWS support | None | S3 + EC2 + CloudFormation |
| Alignment | Unknown (likely BLOSUM62) | BLOSUM62, open −11, extend −1 |
| Output format | Web display | Tab-separated file |
| Reproducibility | Depends on server | Fully reproducible; markers re-derivable |

The classification **logic is identical** to Leino et al. (2021). The marker
positions were derived computationally from the same reference sequences, not
manually reconstructed from figures.

## Known limitations

- **Unclassified sequences:** A substantial fraction of bacterial EPSPS sequences
  do not match any of the four known classes. The original paper notes that
  "further empirical studies are needed to identify novel amino acid markers"
  ; this is an active research frontier. Report the unclassified fraction
  explicitly in any publication.
- **Class III motif completeness:** Class III is defined by 17 sequence
  domain patterns from Carozzi et al. (2006, PCT WO2006/110586, Athenix
  Corporation), identified through experimental screening of glyphosate-tolerant
  bacterial isolates. This implementation searches the query sequence directly
  for each domain using a sliding-window match, consistent with the original
  tool's single-motif threshold. All 17 domains are found in the bvEPSPS
  reference sequence and none produce false positives in vcEPSPS (class I).
  The domain count of 17 reflects the patent's labelling scheme, in which
  three parent domains (I, II, XI) each have named sub-domains counted
  separately; Leino et al. (2021) cite 18 motifs, likely reflecting a
  different counting of these sub-domain relationships. The biological
  content is equivalent.
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
