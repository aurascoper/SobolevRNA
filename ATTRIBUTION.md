# Attribution

This pipeline integrates several open-source models and datasets. The following components are not authored by Hunter Kinder and are used under their respective licenses.

---

## RNA-FM

**Authors:** Chen, J., Hu, Z., Sun, S., Tan, Q., Wang, Y., Yu, Q., Zong, L., Hong, L., Xiao, J., King, I., Wang, B., Yin, B. (2022)
**Paper:** "Interpretable RNA Foundation Model from Unannotated Data for Highly Accurate RNA Structure and Function Predictions" — arXiv:2204.00300
**Kaggle dataset:** aarryyaannnnnn/rna-fm-weights
**Role in pipeline:** Sequence embeddings for contact map generation (HWS input)
**License:** MIT

---

## RibonanzaNet-2

**Author:** shujun717 (Kaggle)
**Kaggle model:** shujun717/ribonanzanet2
**Role in pipeline:** Pairwise contact map prediction (E_DL restraint fallback)

---

## Boltz-1

**Authors:** Wohlwend, J. et al. (2024) — MIT Computational and Systems Biology
**Paper:** "Boltz-1: Democratizing Biomolecular Structure Prediction" — bioRxiv 2024
**Kaggle dataset:** odat1248/boltz-0511
**Role in pipeline:** 3D structure seed generation for targets ≤ 800 nt
**License:** MIT

---

## Protenix v1

**Source:** ByteDance Research (adapted by qiweiyin/Kaggle)
**Kaggle dataset:** qiweiyin/protenix-v1-adjusted
**Role in pipeline:** AlphaFold3-class 3D generation for targets ≤ 512 nt
**License:** Apache 2.0

---

## USalign

**Authors:** Zhang, C., Shine, M., Pyle, A.M., Zhang, Y. (2022)
**Paper:** "US-align: universal structure alignments of proteins, nucleic acids, and macromolecular complexes" — Nature Methods 2022
**Role in pipeline:** TM-score computation for evaluation
**License:** Academic/non-commercial (see Zhang Lab license)

---

## Competition Data

**Source:** Stanford RNA 3D Folding Part 2, Kaggle
**URL:** https://www.kaggle.com/competitions/stanford-rna-3d-folding-2
**License:** Competition-specific terms (non-commercial, academic use)

---

## Original Contributions

The following components were designed and implemented by Hunter Kinder:

- **HWS (Hierarchical Windowed Sensor):** Sliding window embedding extraction and taper-weighted accumulation for sequences exceeding RNA-FM's 1,022 nt limit
- **SHR (Stochastic Hamiltonian Relaxation):** Four-term physics Hamiltonian (E_bond + E_rep + E_DL + E_Rg) with Sobolev H¹ gradient preconditioning via DCT-II spectral filtering, implemented in JAX with 64-bit precision
- **Dynamic parameter schedule:** Length-dependent force constant and step schedule across three sequence-length tiers
- **overlap_stitcher:** Kabsch-based 3D chunk stitching with linear taper blending and A-form helix fallback
- **Adaptive TBM routing:** Identity-threshold routing between template-based modeling and deep learning predictors
- **Multi-GPU device routing:** Deterministic allocation of RNA-FM + RibonanzaNet-2 to cuda:0, Boltz-1 + Protenix to cuda:1

The `hungarian-chain-mapping-hws.ipynb` base notebook was forked from a public Kaggle notebook by amirrezaaleyasin, which provided the Protenix integration scaffolding. All physics engine code, HWS, stitching, and pipeline orchestration logic are original.
