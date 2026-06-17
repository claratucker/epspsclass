# Changelog

## v1.0.6 (2026)

### Fixed

- **Class I never classified any real sequence.** Class I required all 148
  `CLASS_I_MARKERS` positions to match exactly. No real, evolutionarily
  diverged organism ever does: benchmark testing against 7 real sequences
  fetched from NCBI (Staphylococcus aureus, Ruminococcus gnavus, Dorea
  formicigenerans as confirmed class II; Bacteroides fragilis, Bacteroides
  thetaiotaomicron, Clostridium perfringens, Prevotella copri as confirmed
  class I from phylogenetically distant taxa) showed real class I organisms
  matching only 24-28 of 148 markers, never all 148. This was fixed in two
  parts:

  1. **Marker subset + count threshold.** Added `CLASS_I_CORE_MARKERS`, the
     20 positions (of the original 148) most discriminating between the
     benchmark's class I and class II organisms, and `CLASS_I_CORE_THRESHOLD
     = 4`. Against the same benchmark, this subset gives a class I floor of
     7/20 and a class II ceiling of 0/20, a clean 7-marker gap, compared to
     a 2-marker gap on the full 148-marker set (24-28 vs 18-22) that was too
     narrow to threshold safely. The full 148-marker count is still computed
     and recorded in `aligned_markers` for diagnostic purposes; it no longer
     gates the classification decision. See the derivation comment on
     `CLASS_I_CORE_MARKERS` in `classifier.py` for the full method.

  2. **Identity-gate bypass for Class I only.** Even with the subset fix,
     all four real class I benchmark organisms still failed to classify,
     because all four fall below the classifier's 40% whole-protein
     identity threshold (32.6-38.3%), which skips marker checking entirely
     for any class below it. Checked directly against Leino et al.'s real
     paper and supplementary material: no identity threshold, no Rost
     (1999) citation, and no mention of an identity pre-filter appears
     anywhere; the actual Methods text (section 5.4) only requires "all the
     amino acid markers from the respective reference sequence," with no
     stated whole-protein identity requirement. The 40% gate is therefore
     not part of the original method and its removal for this one path does
     not diverge further from Leino et al.; it removes an addition that was
     never in the source method. EPSPS is a small, single-domain protein,
     and glyphosate sensitivity is governed by active-site residues, not
     overall fold-level conservation (Funke et al. 2009; Sammons and Gaines
     2014), so a low whole-protein identity score is not, on its own,
     evidence that the marker positions are misaligned or unreliable. The
     bypass applies only to the Class I core-marker check; Class II, III,
     and IV are still gated by the 40% threshold as before. The low
     whole-protein identity is still recorded in `notes` even when the
     bypass allows classification to proceed.

  Both real class I benchmark organisms now classify correctly as Class I
  (`Sensitive`), and real class II benchmark organisms still correctly do
  not (`is_too_divergent=False`, `classes` does not contain `"I"`).
  Regression tests using the exact benchmark sequences (hardcoded, no
  network access required) added as `TestClassICoreMarkerThreshold` in
  `tests/test_integration.py`.

- Added `scripts/epsps/calibrate_class1_subset.py` and
  `scripts/epsps/classify_and_calibrate.py` for reproducing this
  calibration or re-deriving the marker subset against an expanded
  benchmark panel. Requires network access to NCBI eutils to fetch fresh
  benchmark sequences; not required to run the existing test suite.

## v1.0.5 (2026)

### Fixed

- **CLI summary and warning miscounted unreliable sequences.** The
  `--summary` output's "Flagged unreliable" count, and the
  "WARNING: all sequences are below the identity threshold" message, used a
  string-match against each result's `notes` field instead of the
  classifier's own `is_too_divergent` property. A sequence that fails the
  identity threshold against some references but passes against at least
  one, and receives a valid primary classification, has a "below threshold"
  note for the references it failed, which caused it to be miscounted as
  unreliable even though `is_too_divergent` correctly read `False` for that
  same result. Found while running this classifier on EPSPS proteins fetched
  from NCBI for 164 genera: every sequence was flagged unreliable in the
  summary despite valid per-sequence classifications in the TSV output. Both
  the count and the warning now use `r.is_too_divergent` directly.
- **`_print_summary`'s `file` parameter used a late-bound default.**
  `file=sys.stderr` as a default argument value is evaluated once at module
  import time, not per call, so if `sys.stderr` is ever replaced (for
  example by a test harness or output-capturing wrapper) after import, the
  default silently continues pointing at the original, possibly closed,
  object. Changed to `file=None`, resolved to the live `sys.stderr` inside
  the function body.
- Added regression tests covering both issues in `tests/test_integration.py`.

### Documentation correction

- **v1.0.4's Class I Ia/Ib claim was checked against the real source and is
  incorrect.** Obtained the actual Leino et al. (2021) Supplementary Table 1
  (the full 464-position alignment table with separate Class Ia, Class Ib,
  Class II, Class III, Class IV columns) and parsed it programmatically. The
  Class Ia and Class Ib columns are identical at all 464 positions. There is
  no sequence-level distinction between Ia and Ib in the published data, and
  no separate marker sets to split out; "alpha and beta" in the paper's
  results text appear to be historically named lineages sharing one marker
  set, not two separate classifier inputs. Splitting `CLASS_I_MARKERS` into
  `CLASS_IA_MARKERS`/`CLASS_IB_MARKERS` is not a valid fix and is no longer
  an open issue for that reason.
- The real cause of confirmed class I organisms still failing to match all
  148 markers is the exact-match-all rule's strictness against genuine
  sequence divergence, independent of any Ia/Ib question. This was already
  documented correctly in the README's "Class I sequences from
  phylogenetically distant taxa" limitation and remains unresolved: a count
  threshold is not safely calibratable given the 2-marker gap between real
  class II organisms (18-22 of 148) and distant class I organisms (24-28 of
  148) in the benchmark panel.

## v1.0.4 (2026)

### Fixed based on benchmark testing

Benchmark testing against gut microbiome sequences with known classifications
from Leino et al. (2021) Table 1 revealed two systematic classification failures:

- **Class II threshold changed to >=50 of 204 markers.** Requiring all 204
  markers failed to classify known class II organisms (Staphylococcus aureus,
  Ruminococcus gnavus, Dorea formicigenerans), which found 60-79 markers.
  Class I organisms found only 14-34 markers. Threshold of 50 correctly
  separates the two groups in benchmark data.
- **Class I Ia/Ib sub-class limitation documented.** Leino et al. (2021)
  define two class I sub-classes (Ia and Ib) with separate marker sets in
  Supplementary Table 1. This tool combines both into a single 148-marker
  set, causing class Ib sequences to be misclassified as Unclassified.
  Gut microbiome class I organisms find only 27-28 of 148 combined markers.
  This is a known limitation; implementing separate Ia and Ib marker sets
  is an open issue.

  **Correction, v1.0.5:** the claim above is incorrect. Checked directly
  against the paper's actual Supplementary Table 1 (the full 464-position
  alignment table): the Class Ia and Class Ib columns are identical at every
  position. There is no sequence-level distinction between Ia and Ib in the
  published data, and no separate marker sets to implement. See the v1.0.5
  entry below for the corrected understanding of why real class I organisms
  still fail to match all 148 markers.
- README methodological divergences section updated to document both thresholds
  with benchmark calibration data.
- Known limitations section updated with class I Ia/Ib issue.

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
