# Changelog

## v1.0.1 — 2026-07-20 (documentation / packaging only; frozen scientific artifacts unchanged)
- Removed the manuscript draft and supplement from the public package; the repository now contains data, code, results, tables, figures, and reproducibility documentation only.
- Corrected README "Headline results" wording for claim discipline:
  - collision result stated as "radius-limited rather than hash-induced" (not "not hashing");
  - ridge vs random forest stated as "comparable observed performance ... does not establish formal equivalence" (removed "largely additive signal" shorthand);
  - the tail-based signal statistic was moved out of the headlines into a clearly labeled "Exploratory, post-hoc analysis" section and marked hypothesis-generating, not a frozen confirmatory result.
- Added the Zenodo DOI to README and CITATION.cff; bumped package version to 1.0.1.
- Removed bibliography-verification files (REFERENCES_VERIFIED.md/.json) from the public package.

Frozen scientific artifacts (PROTOCOL_FROZEN.json, cohort_pairs.csv, cohort_results.json, FINAL_RESULTS.json,
confirmatory_manifest.json + .sha256, per-target checkpoints, null_arrays.json, S_interp_vs_extrap.json,
target_metadata_14.csv/.json, figures, and result tables) are BYTE-IDENTICAL to v1.0. No analysis, p-value,
q-value, positive-target list, figure, or numerical result changed.

Why a new Zenodo version: the v1.0 archive (DOI 10.5281/zenodo.21436571) was cut from the initial commit, which
still contained the manuscript draft; v1.0.1 provides a manuscript-free archive that matches the intended
reproducibility package and the corrected documentation wording. Under Zenodo's immutable-versioning model the v1.0
record is preserved in the version history and v1.0.1 supersedes it as the latest version; the concept
(all-versions) DOI resolves to the latest version.

## v1.0 — 2026-07-18
- Initial public reproducibility package (frozen carry-vs-use enantioselectivity audit).
