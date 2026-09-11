# Protein-Templated DNA Synthesis — Simulation Code

**An Information Dynamics Model of Protein-Templated DNA Synthesis: From Qualitative Observations to Quantitative Predictions**

---

## Overview

This repository contains the simulation code supporting the manuscript's
quantitative model of Drt3b protein-templated DNA synthesis. The model
recasts the Drt3b mechanism within the information-dynamics framework:

- **Constraint space** (formerly called "virtual space") — a two-state
  automaton {S_A, S_C} with deterministic alternation S_A → S_C → S_A,
  encoding the geometric constraints of the active pocket.
- **Real space** — the cellular dNTP pool, with state-dependent matching
  energies set by the active-site geometry.
- **Coupling** — a Boltzmann selection rule, P(n | s) ∝ exp(−E(n, s)/kT),
  where the energy function E(n, s) maps the substrate n to the
  constraint-space state s via the explicit state→nucleotide map
  m(S_A) = A, m(S_C) = C.

The framework reproduces wild-type strict alternation, the E26A partial
fidelity loss, and the E26Q dG misincorporation, and it generates
testable predictions for mutants not yet experimentally characterized.

---

## Repository structure

```
protein-templated-dna-sim/
├── README.md
└── Drt3b/
    ├── drt3b_deterministic_simulator.py
    ├── drt3b_deterministic_simulator_log.txt
    ├── drt3b_mutation_effects.py
    ├── drt3b_mutation_effects_log.txt
    ├── drt3b_dinucleotide_complexity.py
    ├── drt3b_dinucleotide_complexity_log.txt
    ├── drt3b_noise_robustness.py
    └── drt3b_noise_robustness_log.txt
```

| Script | Purpose |
|---|---|
| `drt3b_deterministic_simulator.py` | Deterministic greedy simulator in the large-barrier limit (Δ → ∞). Reproduces strict AC alternation. |
| `drt3b_mutation_effects.py` | Two-layer mutation ensemble (A/C-only and G/T-extended) with empirical standard errors over 50 independent seeds. |
| `drt3b_dinucleotide_complexity.py` | Dinucleotide frequencies, three distinct entropy measures, LZ complexity, and the (AC)₆ ↔ error-rate mapping. |
| `drt3b_noise_robustness.py` | Temperature robustness scan; regenerates manuscript Figure 2. |

---

## Model layers

Two model layers are reported **separately**, each with its own metric.
Every output table labels the layer and metric explicitly.

### Layer (A) — A/C-only sub-model

- G and T excluded (E = ∞).
- Per-step mismatch probability: `p = e^(−Δ) / (1 + e^(−Δ))`.
- **Metric**: ε = fraction of adjacent identical bases.
- **Closed form**: `ε(Δ) = 2 e^(−Δ) / (1 + e^(−Δ))²`.

### Layer (B) — G/T-extended sub-model

- G allowed only in the S_A state; T excluded everywhere.
- `Δ_G,A` calibrated by Boltzmann inversion of the E26Q product-level
  dG fraction (10.16%): `Δ_G,A = ln[(1 − p_G)/p_G] ≈ 1.386 kT`
  with `p_G = 0.20`.
- **Metric**: product-level dG fraction = 0.5 × p_G.
- The adjacent-identical metric ε is identically zero in this layer
  (the S_A and S_C emission alphabets are disjoint).

---

## Key results

### Mutation effects (Layer A — A/C-only)

| Mutant | Δ (kT) | ε (mean ± SE) | ε analytic |
|---|---|---|---|
| Wild-type | 100 | 0.0000 ± 0.0000 | 0.0000 |
| E26A | 2.5 | 0.1406 ± 0.0024 | 0.1402 |
| E26Q | 2.5 | 0.1406 ± 0.0024 | 0.1402 |
| R253A (pred.) | 2.5 | 0.1406 ± 0.0024 | 0.1402 |
| E26A_R253A (pred.) | 1.0 | 0.3956 ± 0.0042 | 0.3932 |
| Random (A/C-only) | 0.0 | 0.5017 ± 0.0035 | 0.5000 |
| Random (G/T allowed) | — | 0.2503 ± 0.0028 | — |

### Mutation effects (Layer B — G/T-extended)

| Mutant | Δ_G,A (kT) | dG fraction (mean ± SE) |
|---|---|---|
| E26A | 1.386 | 0.1000 ± 0.0015 |
| E26Q | 1.386 | 0.1000 ± 0.0015 |

Experimental anchor: E26Q product-level dG fraction = 10.16%
(Deng et al. 2026, Fig. S11D/S11F). This is a **calibration** of the
model parameter to the experimental value, not an independent prediction.

### Entropy measures (three distinct quantities)

| Sequence | H₁ (bit) | H₂ (bit) | H(X₂\|X₁) (bit) |
|---|---|---|---|
| WT (strict AC) | 1.0000 | 1.0000 | 0.0000 |
| E26A (Δ = 2.5) | 0.9996 | 1.5814 | 0.5818 |
| E26A_R253A (pred.) | 0.9988 | 1.9625 | 0.9637 |
| Random (A/C-only) | 0.9984 | 1.9956 | 0.9971 |
| poly(AAC), phase unknown | 0.9185 | 1.5850 | 0.6665 |

For a deterministic poly(AAC) repeat: H₁ = 0.918 bit,
H₂ = log₂ 3 = 1.585 bit, H(X₂\|X₁) = 0.667 bit (phase unknown) or
0 bit (phase known). The value 1.58 bit is the **dinucleotide** entropy
H₂, not the marginal entropy H₁.

### Temperature robustness (manuscript Figure 2)

| T (kT) | ε(Δ=100) | ε(Δ=1.0) |
|---|---|---|
| 0.1 | 0.0000 | 0.0001 |
| 0.5 | 0.0000 | 0.2100 |
| 1.0 | 0.0000 | 0.3932 |
| 2.0 | 0.0000 | 0.4700 |
| 4.0 | 0.0000 | 0.4923 |
| 6.0 | 0.0000 | 0.4965 |
| 8.0 | 0.0000 | 0.4981 |
| 10.0 | 0.0001 | 0.4988 |

The wild-type conclusion (ε ≈ 0) is unchanged for any Δ ≳ 10.
Δ = 100 kT is retained only as a numerically safe surrogate for the
hard-exclusion limit (Δ → ∞).

---

## How to run

```bash
cd Drt3b

python3 drt3b_deterministic_simulator.py
python3 drt3b_mutation_effects.py
python3 drt3b_dinucleotide_complexity.py
python3 drt3b_noise_robustness.py
```

**Requirements**: Python 3.11+, `numpy`.

All scripts use 50 independent random seeds (1000–1049) for ensemble
statistics, with a fixed seed (42) reserved for single-run
demonstration. All reported values are reproducible from the
provided code.

---

## Data sources

Experimental anchors used for calibration:

- Deng, P., Lee, H., Armijo, C., Wang, H., Gao, A. (2026).
  Protein-templated synthesis of dinucleotide repeat DNA by an
  antiphage reverse transcriptase. *Science* **392**, 1274–1281.
  Figs. 4H, S11A, S11D, S11F.
- Kiran, S. (2026). A substrate recursion principle for biological
  information, with empirical anchoring through a templating-mode
  taxonomy. Preprint (June 2026).

**Note on R253A**: this mutant has not been experimentally
characterized; the original study did not purify the R253A mutant
protein. R253A parameters in this code are a **model prediction**
based on the structural symmetry between Glu26 (A-state gate) and
Arg253 (C-state gate), and remain to be tested.

---

## Terminology

- **Constraint space** — the abstract space of geometric and chemical
  constraints imposed by the Drt3b active pocket (the two-state
  automaton {S_A, S_C}). Corresponds to the "virtual space" of the
  broader information-dynamics framework.
- **Real space** — the physical space of the dNTP pool.
- **State s ∈ {S_A, S_C}** — a conformational variable of the active
  pocket.
- **Nucleotide n ∈ {A, C, G, T}** — the substrate.
- **m(·)** — the state→nucleotide map m(S_A) = A, m(S_C) = C.

---

## Repository

<https://github.com/hkaiopen/protein-templated-dna-sim>

## License

Code is released for academic and non-commercial use. Please cite the
manuscript if you use this code in published work.

## Contact

Kai Huang — hkaiopen@foxmail.com — ORCID 0009-0007-0898-180X
