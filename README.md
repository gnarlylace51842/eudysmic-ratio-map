# Carried but Conditionally Used — v1.0 release package

Frozen submission package for: *Carried but Conditionally Used: A Frozen Held-Out Audit of Whether Chirality-Aware Molecular Fingerprints Predict Measured Same-Assay Enantioselectivity* (D. Ashraf).

## What this is
A same-assay, target-resolved audit separating whether radius-3 chirality-aware Morgan fingerprints **carry** stereochemistry (a representation property) from whether an exactly antisymmetric model **uses** it to predict measured enantiomer potency differences under scaffold extrapolation. Analysis executed under a protocol frozen and hashed **before** held-out modeling; development targets (BACE1, BTK, JAK3) are excluded from the held-out correction family.

## Headline results
- **Carry:** near-complete at radius 3; radius-2 collisions are a neighborhood-radius effect (config-invariant across 1,024–8,192 bits, count, and unhashed identifiers), not hashing.
- **Use (held-out, primary):** **5 of 11 held-out targets** FDR-positive (JAK1, JAK2, ERK2, ROR-γ, HDAC6).
- **Use (14-target characterization):** 8 of 14 positive (adds BACE1, BTK, EGFR — EGFR is a correction-family boundary case).
- **Model class:** ridge ≈ RF (Wilcoxon p=0.76) → largely additive signal.
- **Signal-to-variability:** frozen median-based R_t suggestive but **not confirmed** on held-out (ρ=0.42, p=0.20); an exploratory tail-based statistic is stronger (held-out ρ=0.77) but **post-hoc**.

## Folder layout
```
manuscript/   MANUSCRIPT_DRAFT.md, SUPPLEMENT.md
figures/      Fig1–5 (main) + SFig_* (supplement), PNG, generated from frozen data
tables/       Table1–2 (main) + STable_* (supplement) + STable_targets (ChEMBL metadata)
code/         all analysis scripts (Python)
results/      PROTOCOL_FROZEN, manifests, cohort_results, FINAL_RESULTS, per-target checkpoints,
              null_arrays, REFERENCES_VERIFIED, FROZEN_ARTIFACTS (hashes)
data/         cohort_pairs.csv (processed same-assay enantiomer pairs; ChEMBL-derived)
README.md  LICENSE  CITATION.cff  requirements.txt
REPRODUCTION.md  DATA_DICTIONARY.md  REFERENCES_VERIFIED.md  HASH_MANIFEST.txt
```

## Reproduce
See `REPRODUCTION.md` for clean-room instructions beginning from ChEMBL 37 and a fresh environment, with exact commands, seeds, versions, runtimes, and expected hashes. `HASH_MANIFEST.txt` lists sha256 for every file.

## Licensing
Code: MIT (see LICENSE). Processed data in `data/` and `results/` derives from ChEMBL 37, released by EMBL-EBI under **CC BY-SA 3.0**; redistribution and reuse of the derived data must comply with that license and cite ChEMBL (ref 9). See LICENSE for details.

## Status
Not yet publicly deposited; a DOI/URL will be added on deposit. This package is an author draft prepared for external review; see the manuscript Limitations and the human-verification list in the accompanying report.
