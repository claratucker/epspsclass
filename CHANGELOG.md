# Changelog

## v1.0.2 (2026)

### Fixed

- Class III classification now uses the 17 experimentally validated domain
  patterns from Carozzi et al. (2006, PCT WO2006/110586, Athenix Corporation)
  rather than computationally derived triplets from the MAFFT alignment.
  Domains are matched by sliding-window search over the raw query sequence,
  consistent with the original EPSPSClass tool approach.
- Threshold restored to >=1 domain match (original definition), replacing
  the provisional >=2 triplet match used in v1.0.0 and v1.0.1. Specificity
  is maintained by the length and degenerate character of the Carozzi domains
  (Domain I is 17 positions; Domain VII is 16 positions) rather than by
  requiring multiple shorter matches.
- All 17 domains validated: bvEPSPS reference correctly self-classifies as
  class III with 17/17 domains found; vcEPSPS (class I) yields 0 false
  positives across all 17 domains.

## v1.0.1 (2026)

### Added

- `--show-mixed` flag on `epspsclass classify`: prints sequences matching
  more than one class simultaneously to stderr, with their full class list,
  primary class, and per-reference identity percentages.
- Mixed-class count added to `--summary` output.
- README section clarifying that mixed-class assignments are preserved in the
  `all_classes` TSV column and that `primary_class` applies the priority rule
  (IV > III > II > I). The Ia/Ib subdivision visible in Leino et al. (2021)
  Supplementary Figure 5 is descriptive (PCA-derived) and not a separate
  classification output of EPSPSClass — both Ia and Ib use vcEPSPS markers
  and are reported as Class I.

### Fixed

- CI workflow: corrected setuptools build backend (`setuptools.build_meta`)
  and updated GitHub Actions to Node.js 24 compatible versions.

## v1.0.0 (2026)

### Initial release

- Open-source reimplementation of the EPSPSClass algorithm (Leino et al. 2021,
  *Environment International* 149:106334)
- Marker positions derived computationally from MAFFT L-INS-i alignment of the
  four reference sequences (Supplementary Table 2, Leino et al. 2021), producing
  a 464-column alignment matching Supplementary Table 1 exactly
- Canonical Pro101/Leu100 sensitivity marker confirmed at alignment column 105
- Classification: Class I (148 markers), Class II (204 markers), Class III
  (21 exclusive triplet motifs, ≥2 required), Class IV (162 markers, ≥10 required)
- BLOSUM62 substitution matrix, gap open −11, extend −1
- Reference sequences bundled from Supplementary Table 2 — fully offline operation
- CLI with stdin/stdout support, TSV/CSV output, documented exit codes
- Python API: `EPSPSClassifier.classify()`, `classify_fasta()`
- AWS support: S3 batch workflow, CloudFormation EC2 stack
- Snakemake rule: `workflow/rules/epspsclass.smk`
- Nextflow DSL2 module: `workflow/modules/epspsclass.nf`
- Conda environment: `workflow/envs/epspsclass.yaml`
- Marker re-derivation script: `scripts/derive_markers.sh`
- 17 unit tests, GitHub Actions CI (Python 3.9–3.12)
