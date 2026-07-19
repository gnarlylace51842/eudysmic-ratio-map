# Clean-room reproduction guide

Reproduce the analysis from ChEMBL 37 and a fresh environment. Two paths: **(A) fast** — reproduce all figures/tables from the shipped processed data (`data/cohort_pairs.csv`), no ChEMBL download; **(B) full** — regenerate the processed data from ChEMBL 37 as well.

## 0. Environment
```
python 3.12
pip install rdkit==2026.03.3 scikit-learn==1.7.2 numpy scipy matplotlib
# (see requirements.txt)
```
All randomness is seeded: **orientation-symmetrization = 7; primary/secondary nulls = 202; supplement null regen = 303; stereo-audit sample = 20260711.** Models are frozen: RF (n_estimators=80, min_samples_leaf=3), ridge (alpha=10); representation = Morgan radius 3, 4096 bits, binary, useChirality=True.

## A. Fast reproduction (no ChEMBL; ~10–20 min total)
Working dir = the release root. Scripts read `data/cohort_pairs.csv` (or expect it in the script directory — copy it next to the scripts, or adjust the `D=` path).
```
python code/make_figs.py     # Fig 1–5 + Table1_characterization_14.csv + Table2_heldout_11.csv     (~2–4 min; re-derives OOF + scaffold-bootstrap CIs)
python code/gen_supp.py      # STable_funnel/provenance/collision/Rt/target_full + SFig_collision/interp_extrap/Rt (~3–5 min)
python code/null_regen.py    # null_arrays.json (300 sign-flip nulls/target)                            (~40–50 min)
# then regenerate SFig_nulls.png from null_arrays.json (panel code in the paper's methods)
```
**Must reproduce (expected, ±rounding):** 14-family positives = {BACE1,BTK,EGFR,ERK2,HDAC6,JAK1,JAK2,ROR-γ}; held-out (m=11) positives = {JAK1,JAK2,ERK2,ROR-γ,HDAC6}; ridge−RF median +0.005, Wilcoxon p=0.76; R_t(median) vs S_Δ all ρ=0.62 / held-out ρ=0.42; collision config-invariant = True for all 14. `null_regen.py` self-validates its regenerated null means against the stored 2,000-rep summaries (max deviation observed: 0.014).

## B. Full reproduction from ChEMBL 37 (adds ~6–9 h of compute)
1. Download `chembl_37_sqlite.tar.gz` from EMBL-EBI (ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/), extract, and place/symlink the SQLite DB at `data/chembl_37.db` (~28 GB extracted).
2. Enumerate the confirmatory manifest (data-only eligibility): `python code/manifest_enumerate.py` → `confirmatory_manifest.json` (~3–8 min; expects the DB path).
3. Freeze the cohort pairs: `python code/build_cohort.py` → `cohort_pairs.csv` (should match the shipped file; sha256 in HASH_MANIFEST).
4. Run the confirmatory cohort (resumable, per-target checkpoints): `python code/cohort_run.py` → `cohort_results.json` + `cohort_ckpt/*.json` (**~3–5 h**; RF is the cost; 2,000 sign-flip + 1,000 permutation nulls × 2 models × 14 targets).
5. Recompute the two BH families → `FINAL_RESULTS.json` (BH code embedded in the results step / make_figs; ~seconds).
6. Regenerate figures/tables and the supplement as in Path A.

**Notes.** The 28 GB DB is not shipped. ChEMBL's `pref_name` string for tid 165 ("…inwardly rectifying…") is a misleading label artifact; the authoritative identity is **KCNH2 (hERG/Kv11.1)** — ChEMBL CHEMBL240, gene KCNH2, UniProt Q12809, *Homo sapiens* (a Kv11.1 voltage-gated channel, not a Kir/KCNJ channel) — see `results/target_metadata_14.csv` and `tables/STable_targets.csv`. Exact floating-point p-values depend on the null draw sequence (seed 202); the shipped `cohort_results.json` is canonical and figures/tables are derived from it, so the fast path reproduces the published numbers exactly, while a full re-run reproduces them within Monte-Carlo error of the empirical nulls.

## Hashes
`HASH_MANIFEST.txt` lists sha256 for every file in this release; `results/FROZEN_ARTIFACTS.json` records the canonical analysis-artifact hashes and seeds.
