# Target-Identity and Release-Integrity Audit

**Verdict: PASS — metadata corrected; scientific results unchanged.**

Date: 2026-07-17. Scope: metadata, naming, and reproducibility-package correction only. No analysis was rerun; no p-value, q-value, positive-target list, development/held-out assignment, figure result, or conclusion was changed. The integrity gate (target-row alignment + no interchange, §6) was run **first**; it passed, so the correction proceeded.

## 1. Method
Target identities were resolved **from the ChEMBL 37 database** (`target_dictionary` for ChEMBL ID / preferred name / organism / target type; `target_components` → `component_sequences` for UniProt accession; `component_synonyms` `GENE_SYMBOL` for gene symbol), keyed on the internal `tid`. Identities were **not** inferred from display strings. Counts (pairs/scaffolds/docs) were re-derived from the frozen `cohort_pairs.csv` and cross-checked against `confirmatory_manifest.json` and `cohort_results.json`.

## 2. Files inspected
Entire package under `pilot_stereo/` and `pilot_stereo/release_v1.0/`: manuscript (`MANUSCRIPT_DRAFT.md`, `MANUSCRIPT_SKELETON.md`), supplement (`SUPPLEMENT.md`), all figures (`figures/*.png`), all tables (`*.csv`, `tables/*.csv`), all results/JSON (`cohort_results.json`, `FINAL_RESULTS.json`, `confirmatory_manifest.json`, `S_interp_vs_extrap.json`, `null_arrays.json`, per-target checkpoints, `FROZEN_ARTIFACTS.json`, `PROTOCOL_FROZEN.json`), all code (`*.py`), and release docs (`README.md`, `REPRODUCTION.md`, `DATA_DICTIONARY.md`, `REFERENCES_VERIFIED.md`, `HASH_MANIFEST.txt`).

## 3. Resolved target identities (14 targets, from ChEMBL 37)
| tid | ChEMBL ID | gene | UniProt | organism | type | display name | split | pairs | scaf | docs |
|---|---|---|---|---|---|---|---|---|---|---|
| 100097 | CHEMBL5251 | BTK | Q06187 | Homo sapiens | SINGLE PROTEIN | BTK | development | 266 | 114 | 56 |
| 10919 | CHEMBL2835 | JAK1 | P23458 | Homo sapiens | SINGLE PROTEIN | JAK1 | held-out | 207 | 119 | 46 |
| 103982 | CHEMBL1741186 | RORC | P51449 | Homo sapiens | SINGLE PROTEIN | ROR-γ (RORC) | held-out | 201 | 119 | 27 |
| 10938 | CHEMBL2971 | JAK2 | O60674 | Homo sapiens | SINGLE PROTEIN | JAK2 | held-out | 199 | 117 | 48 |
| 19639 | CHEMBL3778 | IRAK4 | Q9NWZ3 | Homo sapiens | SINGLE PROTEIN | IRAK4 | held-out | 172 | 82 | 29 |
| 10849 | CHEMBL2148 | JAK3 | P52333 | Homo sapiens | SINGLE PROTEIN | JAK3 | development | 139 | 78 | 45 |
| 165 | CHEMBL240 | KCNH2 | Q12809 | Homo sapiens | SINGLE PROTEIN | KCNH2 (hERG/Kv11.1) | held-out | 129 | 95 | 65 |
| 12252 | CHEMBL4822 | BACE1 | P56817 | Homo sapiens | SINGLE PROTEIN | BACE1 | development | 112 | 59 | 32 |
| 11638 | CHEMBL4040 | MAPK1 | P28482 | Homo sapiens | SINGLE PROTEIN | ERK2 (MAPK1) | held-out | 102 | 64 | 18 |
| 11177 | CHEMBL3130 | PIK3CD | O00329 | Homo sapiens | SINGLE PROTEIN | PI3Kδ (PIK3CD) | held-out | 98 | 57 | 37 |
| 12694 | CHEMBL3553 | TYK2 | P29597 | Homo sapiens | SINGLE PROTEIN | TYK2 | held-out | 95 | 51 | 28 |
| 11307 | CHEMBL1865 | HDAC6 | Q9UBN7 | Homo sapiens | SINGLE PROTEIN | HDAC6 | held-out | 83 | 53 | 27 |
| 9 | CHEMBL203 | EGFR | P00533 | Homo sapiens | SINGLE PROTEIN | EGFR | held-out | 82 | 64 | 51 |
| 10906 | CHEMBL2599 | SYK | P43405 | Homo sapiens | SINGLE PROTEIN | SYK | held-out | 82 | 52 | 21 |

Written to `target_metadata_14.csv` and `target_metadata_14.json` (fields: tid, target_chembl_id, official_preferred_name, gene_symbol, uniprot_accession, organism, target_type, manuscript_display_name, development_or_heldout, n_pairs, n_scaffolds, n_documents).

## 4. tid 165 — the KCNH2/hERG correction
- **tid 165 = ChEMBL `CHEMBL240` = gene `KCNH2` = UniProt `Q12809` = *Homo sapiens*, SINGLE PROTEIN.** Display name: **KCNH2 (hERG/Kv11.1)** at first mention, **KCNH2/hERG** thereafter.
- ChEMBL's stored `pref_name` for this target is the literal string **"Voltage-gated inwardly rectifying potassium channel KCNH2"**. The "inwardly rectifying" descriptor is a **ChEMBL database-label artifact** (hERG/Kv11.1 is a voltage-gated *delayed-rectifier*, not an inwardly-rectifying Kir/KCNJ channel). Authoritative identity is fixed by the ChEMBL ID, gene symbol, and UniProt accession — all of which resolve unambiguously to KCNH2/hERG. The raw ChEMBL string is retained verbatim in the data/result files as faithful provenance and in `official_preferred_name`; the manuscript-facing display name is the standardized **KCNH2 (hERG/Kv11.1)**.
- **KCNH2/hERG is not a positive target.** Its held-out predictive skill is slightly **negative** (S_Δ = -0.007); it is not described as a successful predictive target anywhere.

## 5. Stale `KCNJ` occurrences found and corrected
| file | before | after |
|---|---|---|
| `S_interp_vs_extrap.json` (field `ab`, tid 165) | `"KCNJ"` | `"hERG"` |
| `figures/SFig_interp_extrap.png` (panel label, derived from the JSON) | `KCNJ` | `hERG` (regenerated via `gen_supp.py`) |
| `null_regen.log` (transient run log; not shipped, not in FROZEN_ARTIFACTS) | line `KCNJ real=-0.007 …` | file removed (regenerable) |

No file anywhere claimed "hERG was corrected to a KCNJ/Kir target"; no "not hERG" statement existed. `FINAL_RESULTS.json` had already been relabeled KCNJ→hERG in an earlier pass (numbers identical; see §7).

## 6. Remaining `KCNJ`/`Kir` occurrences (intentional — NOT stale mislabels)
These are correct and were deliberately kept:
- `SUPPLEMENT.md`, `DATA_DICTIONARY.md`, `REPRODUCTION.md`: disambiguation sentences stating KCNH2/hERG **is not** a Kir/KCNJ channel (they explain the misleading ChEMBL pref_name for reviewers).
- `FROZEN_ARTIFACTS.json`: the post-freeze change-ledger recording the `KCNJ→hERG` label correction.

## 7. Numerical integrity — no analysis changed
**7a. Frozen scientific artifacts are byte-identical to their pre-registration hashes:**
| artifact | frozen sha256 (12) | current sha256 (12) | status |
|---|---|---|---|
| PROTOCOL_FROZEN.json | d758fc5e9dd3 | d758fc5e9dd3 | UNCHANGED |
| cohort_pairs.csv | de0f79afc2bf | de0f79afc2bf | UNCHANGED |
| cohort_results.json | 31168ab5ed10 | 31168ab5ed10 | UNCHANGED |
| manifest_enumerate.py | e515459d1b71 | e515459d1b71 | UNCHANGED |
| build_cohort.py | 8c5e66dc5968 | 8c5e66dc5968 | UNCHANGED |
| cohort_run.py | 0177bc0b8822 | 0177bc0b8822 | UNCHANGED |
| audit_reps.py | 3e352fa2777a | 3e352fa2777a | UNCHANGED |
| build_frozen.py | b20b1a832082 | b20b1a832082 | UNCHANGED |
| cohort_ckpt/ (14, combined) | 645d07be77dd | 645d07be77dd | UNCHANGED |

**7b. Target-row alignment (no interchange):** all 14 targets have identical pairs/scaffolds across `cohort_pairs.csv`, `confirmatory_manifest.json`, and `cohort_results.json`; total = 1,967 pairs. Spot-checks confirm the two 82-pair targets are **not** interchanged:
| target | tid | pairs | scaffolds |
|---|---|---|---|
| EGFR | 9 | 82 | 64 |
| SYK | 10906 | 82 | 52 |
| KCNH2 (hERG/Kv11.1) | 165 | 129 | 95 |
| BACE1 | 12252 | 112 | 59 |
| ERK2 (MAPK1) | 11638 | 102 | 64 |
**7c.** Positive-target lists unchanged: 8/14 = BACE1, BTK, EGFR, ERK2, HDAC6, JAK1, JAK2, ROR-γ; 5/11 held-out = ERK2, HDAC6, JAK1, JAK2, ROR-γ. EGFR remains a combined-family (not held-out-family) positive. KCNH2/hERG remains negative (S_Δ = -0.007), not positive.

## 8. Checksum procedure (`confirmatory_manifest.json`)
The manifest was regenerated with immutable identifiers; the old embedded/internal hash was **not** copied into it. A fresh SHA-256 was computed over exactly the saved bytes and stored in a **detached** companion file, then verified by re-hashing.
```
# compute detached checksum over the exact saved bytes
shasum -a 256 confirmatory_manifest.json > confirmatory_manifest.json.sha256
# verify (re-hash and compare)
shasum -a 256 -c confirmatory_manifest.json.sha256      # -> confirmatory_manifest.json: OK
```
- Pre-registration manifest sha256 (retained in `FROZEN_ARTIFACTS.json` as `pre_registration_sha256`): `7b510150b6161f2171748c5a8a29809b2fdf25f92355714ceb7f215392294ea6`
- New manifest sha256 (detached, `confirmatory_manifest.json.sha256`): `9c2d03f917f4284ce3168cb89cf123eefc4f7977c4ac1c20cde56771eae3206b`
- Recomputed now over shipped bytes: `9c2d03f917f4284ce3168cb89cf123eefc4f7977c4ac1c20cde56771eae3206b`  → match: **True**

## 9. Files modified / added
**Modified:** `MANUSCRIPT_DRAFT.md` (ROR-γ→ROR-γ (RORC); stale manifest-hash ref→detached checksum; nomenclature note), `MANUSCRIPT_SKELETON.md` (manifest-hash ref), `SUPPLEMENT.md` (nomenclature + identifier note), `confirmatory_manifest.json` (immutable IDs, no embedded hash), `S_interp_vs_extrap.json` (KCNJ→hERG), `figures/SFig_interp_extrap.png` (regenerated), `FROZEN_ARTIFACTS.json` (documented reconciliation of manifest + FINAL_RESULTS hashes), `release_v1.0/DATA_DICTIONARY.md`, `release_v1.0/REPRODUCTION.md`, `release_v1.0/HASH_MANIFEST.txt` (regenerated, 60 files). Release copies of the above synced.
**Added:** `target_metadata_14.csv`, `target_metadata_14.json`, `confirmatory_manifest.json.sha256` (+ release copies under `results/`).
**Removed:** `null_regen.log` (transient).

## 10. Diff summary of substantive text changes
- Manuscript §Methods: "…internal hash `21e44d41`" → "…records immutable target identifiers; its bytes are covered by the detached checksum `confirmatory_manifest.json.sha256`" + target-nomenclature sentence naming ERK2 (MAPK1), ROR-γ (RORC), KCNH2 (hERG/Kv11.1), PI3Kδ (PIK3CD).
- Manuscript §3.3: "…ROR-γ, and HDAC6" → "…ROR-γ (RORC), and HDAC6".
- Supplement S1/S12: added target-nomenclature + immutable-identifier + detached-checksum notes; explicit KCNH2/hERG disambiguation (not Kir/KCNJ; S_Δ negative).
- Release DATA_DICTIONARY / REPRODUCTION: tid-165 note restandardized to KCNH2 (hERG/Kv11.1), CHEMBL240, KCNH2, Q12809, with metadata-table pointer.
- No sentence expressing a scientific claim, number, or conclusion was altered.

**Verdict: PASS — metadata corrected; scientific results unchanged.** No genuine numerical or target-row mismatch was discovered.
