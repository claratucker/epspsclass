# Changelog

## v1.0.3 (2026)

### Fixed and improved for publication quality

- Corrected classifier.py docstring: alignment has 464 columns not 426
  (vcEPSPS has 426 residues but the alignment with gaps is 464 columns)
- Documented 40% identity threshold as our addition, not from Leino et al.
  (2021), with citations: Rost (1999) Protein Engineering 12:85-94;
  Addou et al. (2009) J Mol Biol 385:1298-1311
- Documented class IV >=10 marker threshold with explicit calibration
  justification (64 of 162 markers accessible in practice due to gap filtering)
- Documented Carozzi sub-domain hierarchy in classifier.py comments and
  validate-markers output (Ia/Ib/Ic are sub-domains of I; IIa/IIb of II)
- Added non-standard amino acid detection: sequences containing IUPAC
  ambiguity codes (B, Z, X, U etc.) now receive a warning note in output
- Added error isolation in classify_fasta(): a malformed sequence no longer
  aborts the entire run; errors are caught, logged, and returned as
  unclassified with an error note
- Added is_too_divergent property to ClassificationResult: distinguishes
  sequences that were below the identity threshold from genuinely unclassified
  sequences. Added as a column in TSV output. For publication, these two
  categories should be reported separately
- Added `group-by-organism` subcommand: aggregates per-sequence classifications
  by genome accession to detect organism-level mixed-class signals (organisms
  carrying EPSPS copies from more than one class). Accepts a metadata TSV
  mapping sequence IDs to organism IDs with optional extra columns
  (taxonomy, habitat, exposure tier) passed through to output
- Added is_too_divergent column to TSV/CSV output
- Fixed Nextflow example in README: was using paired-end FASTQ glob
  (fromFilePairs), corrected to FASTA input (fromPath)
- Created tests/test_integration.py: 15 integration tests covering
  self-classification of all four reference sequences, cross-classification
  checks, identity value sanity checks, and new feature tests. Test count
  increased from 17 to 32

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
  classification output of EPSPSClass; both Ia and Ib use vcEPSPS markers
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
- Reference sequences bundled from Supplementary Table 2, with fully offline operation
- CLI with stdin/stdout support, TSV/CSV output, documented exit codes
- Python API: `EPSPSClassifier.classify()`, `classify_fasta()`
- AWS support: S3 batch workflow, CloudFormation EC2 stack
- Snakemake rule: `workflow/rules/epspsclass.smk`
- Nextflow DSL2 module: `workflow/modules/epspsclass.nf`
- Conda environment: `workflow/envs/epspsclass.yaml`
- Marker re-derivation script: `scripts/derive_markers.sh`
- 17 unit tests, GitHub Actions CI (Python 3.9–3.12)
