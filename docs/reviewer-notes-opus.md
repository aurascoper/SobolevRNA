# Reviewer Notes — Whitepaper Deviations from Spec

Every place the whitepaper departs from the prompt spec, with the code that justifies the deviation. Cell/line references are to `hungarian-chain-mapping-hws(1).ipynb`.

## Resolved discrepancies (spec → code)

### 1. Sobolev α

- **Spec:** α = 10.0 throughout SHR.
- **Code:** `sobolev_h1_smooth(gradient, alpha=10.0)` has default 10.0, but:
  - `shr_refine_single` (cell 24) sets `'alpha': 10.0` — main refinement.
  - `shr_polish` (cell 24) sets `'alpha': 5.0` — seam polish after stitching.
- **Paper (§5.3):** documents both regimes; polish uses stronger smoothing because it runs after Kabsch stitching and targets high-wavenumber seam texture.

### 2. HWS windowing vs 3D-chunk windowing (spec conflation)

- **Spec:** "Window decomposition... W = 1022, stride s = 768" applied uniformly.
- **Code:** Two independent windowing systems exist.
  - *HWS (1D, RNA-FM embeddings)*: cell 19, `HWS_WINDOW_SIZE = 1022`, `HWS_STRIDE = 768`, overlap 254 nt. ✓ matches spec.
  - *3D stitcher (Protenix/Boltz chunks)*: cell 17, `stitch_megascale_rna(window_size=512, overlap=128)`, stride 384. Spec did not mention this.
- **Paper (§6.1):** explicitly separates the two: HWS is 1D embedding stitching, `stitch_megascale_rna` is 3D coordinate stitching, and their window parameters are set by different hardware constraints (RNA-FM positional encoding vs Boltz-1 VRAM).

### 3. Contact map threshold: raw cosine vs pseudo-probability

- **Spec:** threshold θ = 0.85 applied to cosine similarity.
- **Code (cell 19, `embeddings_to_contact_map`):** applies `contact_prob = (sim + 1.0) / 2.0` first, then thresholds at 0.85. This means the effective cosine threshold is 2·0.85 − 1 = 0.70.
- **Paper (§3.5):** writes the threshold correctly on the rescaled pseudo-probability and explicitly notes the equivalent raw-cosine threshold of 0.70.

### 4. Steric regularization ε

- **Spec:** ε = 1e-2 (confirmed correct); ε = 1e-8 is the wrong alternative.
- **Code (cell 24, `energy_steric`):** uses `+ 1e-2`. ✓ matches spec.
- **Paper (§4.3):** reproduces the full gradient-bound argument from the notebook's own comment.

### 5. `BOLTZ_MAX_LEN` vs `MEGASCALE_THRESH`

- **Spec Section 7:** implies "Boltz-1 + Protenix for L ≤ 512 nt where GPU memory permits."
- **Code (cell 26):** `MEGASCALE_THRESH = 512` (routing cutoff to SHR) but `BOLTZ_MAX_LEN = 800` (hard OOM guard; Boltz can still run as an ensemble contributor between 512 and 800 nt even when the main path is SHR).
- **Paper (§7):** documents both constants and explains their different roles.

### 6. Radius-of-gyration target for 9MME

- **Spec:** "For 9MME (N = 4640): R_g^target = 3.5 · 4640^0.45 ≈ 104 Å."
- **Arithmetic:** 3.5 · 4640^0.45 = **156.3 Å**, not 104 Å. Spec's numerical answer is wrong for the stated formula.
- **Code (cell 24, `shr_refine_single`):** `rg_target = 3.5 * (N ** 0.45) if N > 200 else 0.0` — matches the formula, so the computed target for 9MME is 156 Å.
- **Code inline comment (`get_shr_params`):** "Compaction gravity needed to reach Rg ~104 A for 9MME." This comment is inconsistent with the formula the same function uses — likely drift from an earlier coefficient choice.
- **Paper (§4.5):** reports the correct value (156 Å) for the formula as implemented, notes the code comment's discrepancy, discusses the over-target regime, and notes cryo-EM-compatible range (90–130 Å) for reference.

### 7. Bond stiffness `k_bond`

- **Spec:** "k_bond = 100.0."
- **Code:** Module constant `K_BOND = 10.0` (cell 24, line 21) exists but is unused in the live path. `energy_bond(coords, d0=5.95, k_bond=10.0)` has 10.0 as its *default*, but `total_energy` always calls `energy_bond(coords, d0=5.95, k_bond=params['k_bond'])` where `params['k_bond']` comes from `get_shr_params` and is always 100.0 across all three tiers.
- **Paper (§4.2):** reports the running value (100), flags the unused module constant.

### 8. Contact map is blended, not RNA-FM alone

- **Spec:** treats the contact map as coming from HWS only.
- **Code (cell 26, `get_or_build_contact_map`):** performs `cmap = 0.5 * rna_fm_cmap + 0.5 * ribonanza_cmap` when RibonanzaNet-2 is available.
- **Paper (§3.6):** documents the 50/50 blend and explains the complementary information sources (coevolutionary vs chemical-probing-trained).

### 9. Template contact cutoff

- **Spec:** implies `d_contact = 8.0 Å` throughout.
- **Code:** `d_contact = 8.0 Å` in `energy_contacts` (cell 24) is correct for the DL restraint, but `extract_contacts_from_template` (cell 19) uses a different cutoff `d_cutoff = 12.0 Å` for pulling contacts from a reference structure when RNA-FM is unavailable.
- **Paper:** does not dwell on this (the template-contact path is a fallback), but §4.4 is careful to scope "d_contact = 8.0" to the DL restraint only.

### 10. A-form helix fallback parameters

- **Spec:** not mentioned.
- **Code (cell 17, `_aform_helix`):** rise = 2.81 Å, twist = 32.7°, radius = 9.0 Å.
- **Paper (§6.3):** cites the exact values from code, frames the fallback as a stability mechanism rather than a prediction.

## Points where the spec was followed (no deviation)

- `HWS_WINDOW_SIZE = 1022`, `HWS_STRIDE = 768`, `HWS_EMB_DIM = 640` ✓
- `d_0 = 5.95` Å (bond equilibrium) ✓
- `d_contact = 8.0` Å (DL restraint cutoff) ✓
- `min_seq_sep = 6` in contact map ✓
- `chunk_size = 500` for large-N cosine matrix ✓
- Adaptive TBM thresholds 50% / 35% / 20% at L ≤ 1000 / 1000–3000 / > 3000 ✓
- Flory exponent ν = 0.45 ✓
- Kabsch with chirality guard (D[2,2] = sgn(det)) ✓
- Linear taper blending in stitcher ✓

## Numbers the spec stated without justification — I preserved them

- Public score 0.38650, private score 0.50934 — taken as given per user direction (use retrospective final result).
- Oracle ceiling 0.554 — cited in §9; traces to the user's prior gap-analysis memory.
