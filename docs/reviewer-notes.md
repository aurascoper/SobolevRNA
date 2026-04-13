# Reviewer Notes — Spec vs. Code Discrepancy Log

Ground truth: `hungarian-chain-mapping-hws(1).ipynb` (28 cells, verified 2026-04-12)
Spec source: `docs/research-prompt.md` (drafted from code reading, pre-verification)

All 10 discrepancies are resolved and baked into `docs/whitepaper.md` via the updated research prompt.

---

| # | Spec says | Code actually does | Resolution |
|---|-----------|-------------------|------------|
| 1 | Sobolev α = 10.0 throughout | α=10.0 in `shr_refine_single`; α=5.0 in `shr_polish` | Document both. Exploration pass uses α=10.0 (aggressive low-freq bias, coarser moves); polish pass uses α=5.0 (less filtering, finer local refinement). Physically coherent: relax structure at coarse scale first, then refine. |
| 2 | HWS stride s=768, W=1022 | `HWS_STRIDE=768`, `HWS_WINDOW_SIZE=1022`, overlap=254 nt | **Spec correct.** Match confirmed. |
| 3 | Stitcher uses W=1022, s=768 | `stitch_megascale_rna` defaults: `window_size=512`, `overlap=128`, `stride=384` | **Two distinct systems.** Embedding-HWS (for RNA-FM): W=1022, stride=768, overlap=254. 3D-stitcher (for Protenix/Boltz coordinates): W=512, overlap=128, stride=384. Whitepaper separates into §3 (Embedding-HWS) and §3.5 (3D Stitcher). |
| 4 | Contact threshold at raw cosine > 0.85 | Code: `contact_prob = (cos + 1) / 2 ∈ [0,1]`; threshold applied at 0.85 on this rescaled value → equivalent raw cosine threshold = **0.70** | Document the (cos+1)/2 linear rescaling. Threshold 0.85 on contact_prob = raw cosine ≥ 0.70. State both values in whitepaper. |
| 5 | ε = 1e-2 in steric regularization | `dist = sqrt(sum(diff²) + 1e-2)` — **Confirmed** | Match. |
| 6 | BOLTZ_MAX_LEN implied ≤ 512 | `BOLTZ_MAX_LEN = 800`; `MEGASCALE_THRESH = 512` | Routing: L ≤ 512 → Protenix + Boltz-1 ensemble (Boltz can run up to 800 but gated by MEGASCALE_THRESH for ensemble path); L > 512 → SHR megascale path. Use 800/512 split throughout. |
| 7 | R_g(4640) ≈ 104 Å | Formula `3.5 * N^0.45` yields **156 Å** for N=4640 (code comment says ~104 Å — internal comment error) | Use formula value: $R_g^{\text{target}} = 3.5 \times 4640^{0.45} \approx 156$ Å. Flag code comment discrepancy in footnote. Verify against cryo-EM data for large ribosomal subunits (expected R_g for 60S subunit ~100–150 Å — 156 Å is at the upper end but plausible for a 4,640 nt fragment). |
| 8 | k_bond units: kcal/mol/Å² | Pure dimensionless arbitrary energy units. Module-level `K_BOND = 10.0` is **unused**; running value is `k_bond = 100.0` from `get_shr_params` | Clarify: SHR energies are dimensionless. k_bond=100.0 (from dynamic schedule). Module constant K_BOND=10.0 is a dead variable. Do not claim physical units. |
| 9 | TBM thresholds 50/35/20 at L≤1000 / 1000<L≤3000 / L>3000 | `get_adaptive_identity`: 50.0 / 35.0 / 20.0 at same breakpoints — **Confirmed** | Match. |
| 10 | Contact map is RNA-FM (HWS) only | Pipeline blends **50% RNA-FM + 50% RibonanzaNet-2** into the consensus contact map that feeds E_DL | Document the blend in §4.4. E_DL restraints come from a two-model consensus: $\mathcal{C} = 0.5 \cdot C^{\text{HWS}} + 0.5 \cdot C^{\text{RNet2}}$, thresholded to binary. RibonanzaNet-2 provides sequence-based contact prediction as a complementary signal to RNA-FM's embedding-derived contacts. |

---

## Impact on Whitepaper Sections

| Section | Affected by # |
|---------|--------------|
| §3 Embedding-HWS | 2, 4, 10 |
| §3.5 3D Stitcher (new section) | 3 |
| §4.3 Steric repulsion | 5 |
| §4.4 Contact restraint | 4, 6, 10 |
| §4.5 Radius of gyration | 7 |
| §4.6 Bond energy | 8 |
| §5 Sobolev preconditioning | 1 |
| §7 Routing | 6 |

---

## Notes on R_g Discrepancy (#7)

The code comment `# Compaction gravity needed to reach Rg ~104 Å for 9MME` appears to be a stale comment from an earlier version of the formula (possibly when the exponent was ~0.38 or the coefficient was ~2.4). The formula as written — $3.5 \times N^{0.45}$ — is physically defensible (Flory scaling, Hyeon & Thirumalai 2011 empirical exponent) and should be treated as authoritative. The 156 Å value at N=4640 is at the upper bound of reported R_g values for large ribosomal subunit fragments in solution but not implausible. The code enforces this as a *lower-bound* constraint (penalizes only over-extension), so any structure more compact than 156 Å is unpenalized.

---

## Notes on Contact Map Blend (#10)

The actual `get_or_build_contact_map` function implements a three-tier strategy with 50/50 blending when both RNA-FM and RibonanzaNet-2 are available. The blended contact map is then passed as the `contact_map` argument to `total_energy` → `energy_contacts`. This is a more robust design than the spec implied — two orthogonal models (sequence embedding vs. pairwise sequence prediction) averaging out each other's false positives.
