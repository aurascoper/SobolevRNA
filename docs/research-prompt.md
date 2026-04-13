# Deep Research Prompt — RNA 3D Megascale Whitepaper

**Model:** Claude Opus 4.6
**Task:** Analyze the pipeline codebase and produce a rigorous technical whitepaper with full mathematical derivations, GitHub-renderable LaTeX, and biological context.

---

## Instructions for Claude Opus

You are writing a technical whitepaper for a Kaggle competition submission: **RNA 3D Megascale Structure Prediction Pipeline**, which achieved 0.509 private TM-score on the Stanford RNA 3D Folding Part 2 competition.

Read the following source files in full before writing:
- `hungarian-chain-mapping-hws(1).ipynb` — the full pipeline notebook (28 cells)
- `README.md` — architecture overview and mathematical summaries
- `ATTRIBUTION.md` — model credits and original contribution boundaries

**CRITICAL — Read cells selectively.** Infrastructure cells (installs, imports, path setup, model loading) are boilerplate. The original scientific content is concentrated in:

| Cell | Content | Instruction |
|------|---------|-------------|
| 17 | `overlap_stitcher`, `kabsch_align`, `_aform_helix`, `_linear_taper_blend` | **Read in full. Cite function names and variable names.** |
| 19 | `get_rna_fm_embeddings_safe` (HWS), `embeddings_to_contact_map`, `get_contact_map_for_sequence` | **Read in full. Cite function names and variable names.** |
| 24 | `energy_bond`, `energy_steric`, `energy_contacts`, `energy_rg`, `total_energy`, `sobolev_h1_smooth`, `get_shr_params`, `shr_refine_single`, `shr_polish` | **Read in full. Every constant must be verified against code.** |
| 26 | `MEGASCALE_THRESH`, `BOLTZ_MAX_LEN`, routing logic in `get_or_build_contact_map` | Skim for constants and routing structure. |
| 0–16, 18, 20–23, 25, 27 | Installs, imports, model loading, markdown headers | Skip — infrastructure only. |

Your output is a whitepaper in GitHub-flavored Markdown with LaTeX math (use `$...$` for inline math and `$$...$$` for display math — these render on GitHub via MathJax). Do NOT use `\[...\]` or `\(...\)` — GitHub does not render those.

---

## Ground-Truth Corrections (verified against code 2026-04-12)

**These override anything in README.md or the spec sections below. Do not reconcile silently — cite both values where a discrepancy exists.**

| # | Topic | Correct value | Note |
|---|-------|--------------|------|
| 1 | Sobolev α | `shr_refine_single`: α=10.0; `shr_polish`: α=5.0 | Document both. Exploration uses α=10.0 (aggressive low-freq bias); polish uses α=5.0 (finer moves). Physically: coarse relaxation then local refinement. |
| 2 | HWS parameters | `HWS_WINDOW_SIZE=1022`, `HWS_STRIDE=768`, overlap=254 nt | Confirmed correct. |
| 3 | Two distinct systems | **Embedding-HWS** (RNA-FM): W=1022, stride=768, overlap=254. **3D-Stitcher** (coordinates): W=512, overlap=128, stride=384 | Whitepaper must separate these into distinct subsections. They are not the same system. |
| 4 | Contact threshold | Code rescales cosine: `contact_prob = (cos + 1) / 2 ∈ [0,1]`, then thresholds at 0.85. Equivalent raw cosine threshold = **0.70** | State both: threshold 0.85 on contact_prob = raw cosine ≥ 0.70. |
| 5 | Steric ε | `dist = sqrt(sum(diff²) + 1e-2)` | Confirmed. ε=1e-2. |
| 6 | Routing split | `MEGASCALE_THRESH = 512`, `BOLTZ_MAX_LEN = 800` | L ≤ 512: Protenix + Boltz-1 ensemble (Boltz capacity up to 800 nt); L > 512: SHR megascale. Use 800/512 split. |
| 7 | R_g at N=4640 | Formula `3.5 * N^0.45` gives **≈156 Å**, not 104 Å. Code comment `~104 Å` is a stale internal error. | Use 156 Å. Add footnote flagging the code comment discrepancy. |
| 8 | k_bond units | Dimensionless arbitrary energy. Module constant `K_BOND=10.0` is **unused dead code**. Running value from `get_shr_params`: `k_bond=100.0`. | Do not claim kcal/mol/Å². State dimensionless. k_bond=100.0. |
| 9 | TBM thresholds | 50.0/35.0/20.0 at L≤1000 / 1000<L≤3000 / L>3000 | Confirmed correct. |
| 10 | Contact map source | Pipeline blends **50% RNA-FM + 50% RibonanzaNet-2** into consensus contact map for E_DL | Not RNA-FM only. Two-model consensus: $\mathcal{C} = 0.5 C^{\text{HWS}} + 0.5 C^{\text{RNet2}}$, thresholded to binary. |

---

## Whitepaper Structure

Write each section in full. Do not summarize or placeholder. The output will be uploaded directly to `docs/whitepaper.md`.

---

### Section 1: Abstract (250–350 words)

Summarize: the problem (RNA 3D structure prediction at megascale), the key contributions (HWS, SHR, Kabsch stitching, multi-model ensemble), the result (0.387 public → 0.509 private TM-score), and why the inversion is scientifically meaningful.

---

### Section 2: Background and Motivation

Cover:
- Why RNA 3D structure prediction is harder than protein folding (flexibility, non-canonical base pairs, pseudoknots, scarcity of experimental data)
- What TM-score measures and why it is the appropriate metric for this task (cite Zhang & Skolnick 2004)
- The megascale challenge: why models like RNA-FM, Boltz-1, and Protenix truncate at 512–1022 nt, and why this is a fundamental limitation for ribosomal and viral RNA targets
- Reference 9MME specifically (4,640 nt, large ribosomal subunit fragment) as the canonical stress test

---

### Section 3: HWS — Hierarchical Windowed Sensor

Derive the full mathematical formulation from first principles:

1. **The truncation problem:** RNA-FM processes sequences as token strings with positional encodings bounded at length $W = 1022$. For $L > W$, naive truncation loses long-range covariation signal. Quantify the information loss.

2. **Window decomposition:** Partition $[0, L)$ into overlapping windows $\{[a_w, b_w)\}$ with stride $s = 768$ and window size $W = 1022$, yielding $n_w = \lceil (L - W) / s \rceil + 1$ windows.

3. **Taper weight derivation:** Explain why uniform averaging at boundaries causes force discontinuities in the downstream SHR integrator (derivative discontinuity → high-frequency gradient noise). Derive the linear taper as the minimal-smoothness boundary condition.

4. **Weighted accumulation formula** (expand fully):
$$\tilde{e}(i) = \frac{\sum_{w:\, i \in [a_w, b_w)} \omega_w(i)\, e_w(i - a_w)}{\sum_{w:\, i \in [a_w, b_w)} \omega_w(i)}$$

5. **Contact map construction:** Explain the choice of cosine similarity as a proxy for base-pairing propensity, and justify the $\theta = 0.85$ threshold. Discuss minimum sequence separation $|i - j| \geq 6$ (avoids trivial stack contacts).

6. **Complexity:** Analyze memory and compute costs. Why chunked cosine similarity ($chunk\_size = 500$) is necessary for $N > 2000$.

---

### Section 4: SHR — Stochastic Hamiltonian Relaxation

This is the core mathematical contribution. Derive every term from physical principles.

#### 4.1 Coordinate System

Explain why C1′ (ribose sugar C1-prime) atoms are the canonical reduced representation for RNA backbone modeling — one atom per nucleotide, captures backbone geometry without all-atom overhead.

#### 4.2 Bond Energy

$$E_{\text{bond}}(\mathbf{x}) = k_{\text{bond}} \sum_{i=1}^{N-1} \left(\|\mathbf{x}_{i+1} - \mathbf{x}_i\|_2 - d_0\right)^2$$

Derive from a harmonic approximation to the C1′–C1′ virtual bond potential. Justify $d_0 = 5.95$ Å (measured from PDB RNA structures — average C1′–C1′ distance across dinucleotide steps). Discuss the choice $k_{\text{bond}} = 100.0$ (stiff bond prevents chain breakage during compaction; units: kcal/mol/Å² when combined with the learning rate schedule).

#### 4.3 Steric Repulsion

$$E_{\text{rep}}(\mathbf{x}) = \sum_{\substack{i < j \\ |i-j| > 1}} \left(\max\!\left(0,\; \sigma_{\text{clash}} - \|\mathbf{x}_i - \mathbf{x}_j\|_2\right)\right)^2$$

Derive as a soft-sphere excluded volume potential. Explain why the $|i - j| > 1$ mask is necessary (bonded neighbors are governed by $E_{\text{bond}}$, not repulsion). Discuss the gradient explosion problem at $\|\mathbf{x}_i - \mathbf{x}_j\| \to 0$:

The gradient of $E_{\text{rep}}$ w.r.t. $\mathbf{x}_i$ scales as:
$$\left\|\nabla_{\mathbf{x}_i} E_{\text{rep}}\right\| \sim \frac{2\sigma_{\text{clash}}}{\sqrt{\epsilon}}$$

For $\epsilon = 10^{-8}$: max gradient $\approx 60{,}000$ — causes explosion. For $\epsilon = 10^{-2}$: max gradient $\approx 60$ — stable. Prove this bound analytically, and justify why physical RNA structures never have $\|\mathbf{x}_i - \mathbf{x}_j\| < 1$ Å (minimum van der Waals contact distance).

#### 4.4 Deep Learning Contact Restraint

$$E_{\text{DL}}(\mathbf{x}) = w_{\text{DL}} \sum_{(i,j) \in \mathcal{C}} \left(\max\!\left(0,\; \|\mathbf{x}_i - \mathbf{x}_j\|_2 - d_{\text{contact}}\right)\right)^2$$

Explain the contact set $\mathcal{C}$ as predicted by HWS (primary) or RibonanzaNet-2 (fallback). Justify $d_{\text{contact}} = 8.0$ Å (C1′–C1′ distance at which Watson-Crick and non-canonical base pairs are formed — derived from PDB statistics). The hinge loss $\max(0, d - d_c)^2$ imposes an upper-bound constraint, not a fixed-distance constraint, allowing atoms to be closer than $d_c$ without penalty.

#### 4.5 Radius of Gyration (Flory Scaling)

$$E_{R_g}(\mathbf{x}) = k_{R_g} \left(\max\!\left(0,\; R_g^{\text{target}} - R_g(\mathbf{x})\right)\right)^2$$

$$R_g(\mathbf{x}) = \sqrt{\frac{1}{N} \sum_{i=1}^{N} \|\mathbf{x}_i - \bar{\mathbf{x}}\|_2^2}$$

Derive the Flory scaling law for RNA:
$$R_g^{\text{target}} = 3.5 \cdot N^{0.45}$$

Connect to polymer physics: for an ideal polymer (Gaussian chain), $R_g \sim N^{0.5}$. RNA is more compact due to base stacking and tertiary contacts — empirical exponent $\nu \approx 0.45$ from Hyeon & Thirumalai (2011). Explain why this term uses $\max(0, R_g^{\text{target}} - R_g)$ (penalizes over-extension only — already-compact structures are not penalized for being compact).

For 9MME ($N = 4640$): $R_g^{\text{target}} = 3.5 \times 4640^{0.45} \approx 104$ Å. Verify this is consistent with experimental cryo-EM data for large ribosomal subunits.

---

### Section 5: Sobolev H¹ Gradient Preconditioning

This is the most mathematically novel section. Derive from functional analysis.

#### 5.1 The Problem with $L^2$ Gradients

Standard gradient descent uses:
$$\mathbf{x}^{(t+1)} = \mathbf{x}^{(t)} - \eta \nabla H(\mathbf{x}^{(t)})$$

For a chain of $N$ atoms, the gradient $\nabla H \in \mathbb{R}^{N \times 3}$ is a sequence indexed by residue position $i$. The $E_{\text{DL}}$ and $E_{\text{rep}}$ terms generate gradients with high spatial frequency components (oscillating sign, residue-by-residue). These high-frequency modes:
- Compete with and mask the low-frequency compaction signal from $E_{R_g}$
- Cause the Kabsch-aligned chain to develop local kinks during SHR polish
- Require very small learning rates, slowing convergence

#### 5.2 Sobolev Seminorm and H¹ Space

Define the discrete H¹ Sobolev seminorm on sequences $f: \{0, \ldots, N-1\} \to \mathbb{R}$:

$$|f|_{H^1}^2 = \sum_{k=0}^{N-1} (1 + \alpha k^2) |\hat{f}_k|^2$$

where $\hat{f}_k$ are the DCT-II coefficients of $f$ and $\alpha > 0$ is a smoothing parameter. The H¹ inner product induces a preconditioned gradient descent in which high-wavenumber components are penalized by $(1 + \alpha k^2)$.

#### 5.3 The Preconditioned Update

The H¹-preconditioned gradient of $H$ is:

$$\widetilde{\nabla H}(i) = \text{IDCT-II}\!\left(\frac{(\text{DCT-II}\; \nabla H)_k}{1 + \alpha k^2}\right)$$

This is equivalent to solving the linear system $(I - \alpha \Delta)\widetilde{\nabla H} = \nabla H$ where $\Delta$ is the discrete Laplacian on the chain — i.e., it applies an implicit smoothing step via the resolvent of the Laplacian.

The full JAX update step:
$$\mathbf{x}^{(t+1)} = \mathbf{x}^{(t)} - \eta(t) \cdot \text{clip}\!\left(\widetilde{\nabla H}^{(t)},\; -c,\; c\right)$$

with $c \in \{2.0, 5.0\}$ (tier-dependent clip threshold) and linear learning rate decay $\eta(t) = \eta_0(1 - t/T)$.

#### 5.4 Why DCT-II (Not FFT)

The chain has Neumann boundary conditions (free ends, $\nabla \cdot \mathbf{x}|_{\partial} = 0$). DCT-II is the natural spectral basis for Neumann BCs on $[0, N)$, while DFT assumes periodic BCs. Using DFT would introduce Gibbs ringing at chain termini — DCT-II eliminates this artifact.

---

### Section 6: Kabsch Alignment and Overlap Stitching

Provide the full derivation of the Kabsch algorithm as implemented in `overlap_stitcher`.

#### 6.1 Problem Statement

Given predicted chunk $\mathbf{M} \in \mathbb{R}^{m \times 3}$ (mobile) and the assembled global array $\mathbf{T} \in \mathbb{R}^{m \times 3}$ in the overlap region (target), find the rigid body transform $(\mathbf{R}, \mathbf{t})$ minimizing:

$$\underset{\mathbf{R} \in SO(3),\, \mathbf{t} \in \mathbb{R}^3}{\arg\min} \sum_{i=1}^{m} \|\mathbf{R}\mathbf{m}_i + \mathbf{t} - \mathbf{T}_i\|^2$$

#### 6.2 Closed-Form Solution

1. Center both arrays: $\mathbf{M}' = \mathbf{M} - \boldsymbol{\mu}_M$, $\mathbf{T}' = \mathbf{T} - \boldsymbol{\mu}_T$
2. Compute cross-covariance: $\mathbf{H} = \mathbf{M}'^{\top} \mathbf{T}' \in \mathbb{R}^{3 \times 3}$
3. SVD: $\mathbf{H} = \mathbf{U} \boldsymbol{\Sigma} \mathbf{V}^{\top}$
4. Chirality guard: $\mathbf{D} = \text{diag}(1, 1, \text{sgn}(\det(\mathbf{V}^{\top}\mathbf{U}^{\top})))$
5. Optimal rotation: $\mathbf{R}^* = \mathbf{V} \mathbf{D} \mathbf{U}^{\top}$
6. Translation: $\mathbf{t}^* = \boldsymbol{\mu}_T - \mathbf{R}^* \boldsymbol{\mu}_M$

The chirality guard in step 4 ensures $\det(\mathbf{R}^*) = +1$ (proper rotation). Without it, the SVD may return an improper rotation (reflection) when the overlap region is nearly planar — which would invert the chain handedness and produce an unphysical mirror structure.

#### 6.3 Seam Blending

After alignment, the overlap region (positions $[a_w, a_w + m)$) contains two estimates: $\mathbf{T}_{\text{prev}}$ (from prior chunks) and $\mathbf{M}_{\text{new}}^*$ (Kabsch-aligned). They are blended with a linear ramp:

$$\mathbf{x}^{\text{blend}}(i) = \left(1 - \frac{i - a_w}{m}\right) \mathbf{T}_{\text{prev}}(i) + \frac{i - a_w}{m}\, \mathbf{M}_{\text{new}}^*(i)$$

Explain why this is necessary: even after Kabsch alignment, residual RMSD in the overlap region creates a positional discontinuity. Without blending, E_bond at the seam boundary experiences a sudden jump, causing gradient spikes in the SHR polish.

---

### Section 7: Adaptive TBM and Multi-Model Routing

Describe the three-tier routing logic:

1. **TBM (Template-Based Modeling):** For sequences with $\geq 50\%$ identity to known PDB structures, template coordinates are adapted via pairwise sequence alignment (Bio.Align.PairwiseAligner). The identity threshold is adaptive: $\tau(L) = 50\%$ for $L \leq 1000$, $35\%$ for $L \in (1000, 3000]$, $20\%$ for $L > 3000$ (data-sparse region).

2. **Deep learning ensemble:** Boltz-1 + Protenix (5 diverse samples each) for $L \leq 512$ nt where GPU memory permits. Explain why ensemble diversity matters for TM-score optimization vs. RMSD minimization.

3. **SHR megascale:** For $L > 512$ nt, SHR with HWS contact restraints is the primary path. Explain why DL models fail at this scale (OOM on T4 × 2, and positional encoding generalization failures).

---

### Section 8: Results and Analysis

#### 8.1 Score Table

| Set | TM-score |
|-----|---------|
| Public leaderboard | 0.38650 |
| **Private leaderboard** | **0.50934** |

#### 8.2 Why Private > Public (The Inversion)

Analyze the score gap. Hypotheses to develop:
1. **Benchmark composition:** The public test set may over-represent short, well-templated RNA sequences where TBM baselines score well. The private set likely includes more megascale targets (ribosomal, viral) where the HWS + SHR path dominates.
2. **Physics generalization:** $E_{\text{DL}}$ is conditioned on RNA-FM embeddings, which capture sequence coevolution. For novel RNA families with no template, physics-based compaction (E_Rg + E_bond) provides a constraint that pure DL cannot.
3. **Overfitting avoidance:** No task-specific fine-tuning was performed. The pipeline's parameters ($d_0$, $\sigma_{\text{clash}}$, $\alpha$) are derived from physical constants and RNA structural statistics, not optimized on the public test set.

---

### Section 9: Limitations and Future Work

Be honest. Cover:
- HWS contact maps are binary (threshold $\theta = 0.85$) — continuous contact probabilities would improve E_DL
- SHR is a CPU/JAX optimizer, not end-to-end differentiable with the DL predictors — joint training is unexplored
- Kabsch stitching accumulates error across many chunks for very long sequences (9MME at 4,640 nt has ~6 chunks) — a global bundle adjustment would reduce drift
- The Flory exponent $\nu = 0.45$ is from solution NMR data; under cryo-EM conditions the effective compaction may differ
- No multi-sequence alignment (MSA) was used for Protenix/Boltz-1 calls (disabled for speed) — enabling MSA for mid-length targets would likely improve ensemble quality

---

### Formatting Requirements

- All math: GitHub-flavored markdown LaTeX (`$...$` inline, `$$...$$` display). No `\[...\]`.
- All equations numbered: use `(1)`, `(2)`, etc. inline after the equation.
- Each section should be 400–900 words minimum (technical depth, not padding).
- Code snippets where helpful: use triple-backtick Python blocks.
- No bullet lists for mathematical content — use full sentences with proper notation.
- References formatted as: Author (Year) — *Title* — Venue.

---

### Output

Write the complete whitepaper to `docs/whitepaper.md`. Begin immediately after reading the source files. Do not ask clarifying questions — use your best scientific judgment on anything underspecified. If a derivation requires an assumption, state it explicitly and defend it.
