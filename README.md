# Carried but Conditionally Used — v1.0 reproducibility package

Reproducibility package for the study *Carried but Conditionally Used: A Frozen Held-Out Audit of Whether Chirality-Aware Molecular Fingerprints Predict Measured Same-Assay Enantioselectivity* (D. Ashraf). This repository contains the data, code, and result artifacts only; the manuscript itself is not included here.

## What this is
A same-assay, target-resolved audit separating whether radius-3 chirality-aware Morgan fingerprints **carry** stereochemistry (a representation property) from whether an exactly antisymmetric model **uses** it to predict measured enantiomer potency differences under scaffold extrapolation. Analysis executed under a protocol frozen and hashed **before** held-out modeling; development targets (BACE1, BTK, JAK3) are excluded from the held-out correction family.

## Headline results
- **Carry:** Most radius-2 stereochemical collisions in the audited development targets were radius-limited rather than hash-induced, and stereochemical carry was approximately complete at radius 3 (collision rate configuration-invariant across 1,024–8,192 bits, count, and unhashed identifiers).
- **Use (held-out, primary):** **5 of 11 held-out targets** were FDR-positive under scaffold extrapolation (JAK1, JAK2, ERK2/MAPK1, ROR-γ/RORC, HDAC6).
- **Use (14-target characterization):** 8 of 14 positive (adds BACE1, BTK, EGFR); EGFR is positive under the combined 14-target correction family but **not** the held-out-only family.
- **KCNH2/hERG:** slightly negative S_Δ; **not** a positive predictive target.
- **Model class:** Ridge and random forest showed comparable observed performance. This suggests that much of the accessible signal may be additive in signed-fingerprint space, but does not establish formal equivalence.
- **Signal-to-variability (frozen):** the median-based R_t relationship is suggestive but **not confirmed** on the held-out targets (ρ=0.42, p=0.20).

## Exploratory, post-hoc analysis (not a frozen headline)
A tail-based signal statistic (upper-quartile |Δy| / 0.45) correlates more strongly with skill across the map (held-out ρ=0.77). It was chosen *after* observing the HDAC6 counterexample and is therefore hypothesis-generating and post-hoc — **not** a frozen confirmatory result. It is reported only as a concrete hypothesis for future, larger-scale confirmatory work.

## Folder layout
```
code/      all analysis scripts (Python)
data/      cohort_pairs.csv (processed same-assay enantiomer pairs; ChEMBL-derived)
results/   PROTOCOL_FROZEN, confirmatory_manifest (+ .sha256), cohort_results, FINAL_RESULTS,
           per-target checkpoints, null_arrays, S_interp_vs_extrap, target_metadata_14,
           FROZEN_ARTIFACTS (hashes)
tables/    Table1–2 (main) + STable_* + STable_targets (ChEMBL target metadata)
figures/   Fig1–5 (main) + SFig_*, PNG, generated from frozen data
README.md  LICENSE  CITATION.cff  requirements.txt
REPRODUCTION.md  DATA_DICTIONARY.md  TARGET_IDENTITY_AND_CHECKSUM_AUDIT.md  HASH_MANIFEST.txt
```

## Reproduce
See `REPRODUCTION.md` for clean-room instructions beginning from ChEMBL 37 and a fresh environment, with exact commands, seeds, versions, runtimes, and expected hashes. `DATA_DICTIONARY.md` documents every field of every data/result file. `HASH_MANIFEST.txt` lists a sha256 for every file; verify with `shasum -a 256 -c HASH_MANIFEST.txt`.

## Licensing
Code: MIT (see LICENSE). Processed data in `data/` and `results/` derives from ChEMBL 37, released by EMBL-EBI under **CC BY-SA 3.0**; redistribution and reuse of the derived data must comply with that license and cite ChEMBL (Zdrazil et al., *Nucleic Acids Res* 2024;52(D1):D1180–D1192; doi:10.1093/nar/gkad1004). See LICENSE for details.

## Integrity & status
Publicly available in this repository (https://github.com/gnarlylace51842/eudysmic-ratio-map). Every file's SHA-256 is listed in `HASH_MANIFEST.txt`; the prospectively frozen protocol and artifact hashes are in `results/PROTOCOL_FROZEN.json` and `results/FROZEN_ARTIFACTS.json`; target-identity resolution and checksum provenance are documented in `TARGET_IDENTITY_AND_CHECKSUM_AUDIT.md`. Persistent archive: Zenodo **https://doi.org/10.5281/zenodo.21436571**.
