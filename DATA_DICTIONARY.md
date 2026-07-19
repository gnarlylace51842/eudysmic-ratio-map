# Data dictionary

## data/cohort_pairs.csv — one row per unique structural enantiomer pair (n=1,967)
| column | type | description |
|---|---|---|
| tid | int | ChEMBL target id (see tables/STable_targets.csv) |
| target | str | ChEMBL preferred name |
| gkey | str | connectivity+relative-config key (InChI with /m,/s stripped); groups the two enantiomers |
| smiles_A, smiles_B | str | parent (salt-stripped) canonical SMILES of the two enantiomers |
| pchembl_A, pchembl_B | float | per-hand pChEMBL (=−log10 molar IC50), median-aggregated within assay |
| dy | float | Δy label = median over the pair's assays of pChEMBL(A)−pChEMBL(B) (log10 units) |
| scaffold | str | Bemis–Murcko scaffold SMILES (achiral); the extrapolation grouping unit |
| n_stereocenters | int | tetrahedral stereocentres |
| n_assays | int | number of assays contributing to the median Δy |
| n_docs | int | number of source documents contributing |

## results/cohort_results.json — per target (list under "targets")
`tid, name, n, n_scaffold, median_absdy, collision_r3, tie_frac, Rt`; and per model (`rf`,`ridge`): `real{S,spear,bal}`, `S_null_mean`, `S_null_q` (2.5/50/97.5 of the 2,000-rep sign-flip null), `S_p_primary` (empirical p vs sign-flip null), `S_real_minus_null`, `S_p_secondary` (vs 1,000-rep pair-permutation null). S = S_Δ = 1 − MAE(Δ̂,Δy)/MAE(0,Δy).

## results/FINAL_RESULTS.json
`development_targets`; `family_14` (each: name, dev, n, nsc, S, spear, p, q14, Rt); `family_11_heldout` (name, n, S, spear, p, q11); `positives_14`; `positives_11_heldout`. q14/q11 = Benjamini–Hochberg q within the 14-target and 11-held-out-target families respectively.

## results/null_arrays.json — per target
`ab` (label), `real_S`, `null` (300 regenerated sign-flip null S_Δ), `regen_mean`/`stored_mean` and `regen_q975`/`stored_q975` (validation), `stored_p` (from the frozen 2,000-rep run).

## tables/STable_targets.csv
`tid, abbr, pref_name_chembl, gene, uniprot, organism, chembl_id, status`. Note: tid 165 is **KCNH2 (hERG/Kv11.1)** — ChEMBL target **CHEMBL240**, gene **KCNH2**, UniProt **Q12809**, *Homo sapiens* single protein (a voltage-gated Kv11.1 channel, **not** a Kir/KCNJ inwardly-rectifying channel). The "inwardly rectifying" phrase inside ChEMBL's `pref_name` string is a database-label artifact; identity is fixed by the ChEMBL ID / gene / UniProt accession. Full immutable identifiers for all 14 targets are in `results/target_metadata_14.csv`/`.json`.

## results/target_metadata_14.csv / .json
Immutable target identity table (one row per target): `tid, target_chembl_id, official_preferred_name, gene_symbol, uniprot_accession, organism, target_type, manuscript_display_name, development_or_heldout, n_pairs, n_scaffolds, n_documents`. Identities are taken from the ChEMBL target dictionary / component tables (not inferred from display names). Compact gene-symbol labels used in figures/tables map to these rows.

## results/confirmatory_manifest.json (+ .sha256)
Data-only eligibility manifest with immutable identifiers per qualifying target (ChEMBL ID, gene, UniProt, organism, target type, display name) plus `pairs/scaffolds/docs`; eligibility `criteria`; 14 targets, 1,967 pairs, 3 development (BACE1/BTK/JAK3) + 11 held-out. Its exact bytes are covered by the **detached** checksum `confirmatory_manifest.json.sha256` (verify: `shasum -a 256 -c confirmatory_manifest.json.sha256`).

## Units & conventions
pChEMBL and Δy are in log10 potency units. σ(Δy) ≈ 0.45 is an approximate secondary reference (see manuscript §Methods). Positive = q ≤ 0.05 AND S_Δ > 0. "held-out" = not among the three development targets (BACE1, BTK, JAK3).
