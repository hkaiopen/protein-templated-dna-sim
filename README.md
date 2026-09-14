# Protein-Templated DNA Synthesis — Information-Dynamics Simulations

The model treats the Drt3b active pocket as a deterministic two-state
automaton in **constraint space** (`S_A → S_C → S_A`), with nucleotide
emission governed by a **Boltzmann selection rule**. All energies are in
**kT units (kT = 1)**.

The state chain is deterministic (period 2); the nucleotide emission is
stochastic. The chain admits a unique **stationary occupation
distribution** π(S_A) = π(S_C) = 1/2, not a static fixed point.

---

## Repository structure

| File | Purpose |
|---|---|
| `config.py` | Energy barriers, seeds, analytical helpers |
| `drt3b_two_state_simulator.py` | Per-mutant single-sequence simulator with analytical cross-check (fixed seed 42) |
| `drt3b_stochastic_simulator.py` | 50-replicate ensemble for mean ± std |
| `drt3b_gt_extension.py` | G/T-extended model (E26Q only; two-state A/G approximation) |
| `drt3b_dinucleotide_complexity.py` | Analytical entropy metrics |
| `drt3b_delta_scan.py` | ε vs Δ at T = 1 (analytical + ensemble) |
| `drt3b_temperature_scan.py` | ε vs T at fixed Δ (exploratory sanity check) |
| `validate_manuscript.py` | One-command validation of all numerical claims in the manuscript |

Requires Python 3.8+. No external dependencies beyond the standard
library.

---

## Scientific status of reported values

Values in this repository reflect three distinct epistemic categories:
**calibration**, **prediction**, and **model limit**. They should not be
conflated.

| Result | Status |
|---|---|
| WT strict AC alternation | Model consequence (hard-exclusion limit, Δ → ∞) |
| E26A ε ≈ 0.144 | **Calibrated** to the experimental alternation defect |
| R253A ε ≈ 0.144 | **Symmetry-based prediction** (not experimentally tested) |
| E26A_R253A ε ≈ 0.393 | **Residual-barrier prediction** (Δ = 1.0 working assumption) |
| Random (4-base) ε = 0.25 | **Theoretical limit** (Δ = 0) |
| E26Q product-level dG fraction = 0.10 | **Calibrated reproduction** of the reported 10.16% |
| Δ at ε = 0.05 | **Analytical threshold** (Δ ≈ 3.64 kT) |

---

## Key parameter values

| Parameter | Value | Basis |
|---|---|---|
| `DELTA_E26A` | 2.47 kT | Calibrated from experimental ε = 0.144 |
| `DELTA_R253A` | 2.47 kT | Symmetry with E26A (model prediction) |
| `DELTA_DOUBLE` | 1.0 kT | Working assumption for both-barrier-loss |
| `DELTA_RANDOM` | 0.0 kT | Thermodynamic limit |
| `P_G_AT_A_STATE_E26Q` | 0.20 | Structural ~80/20 dA/dG preference |
| `P_G_AT_A_STATE_E26A` | 0.0 | Gly248 gate still excludes G/T |
| `GAMMA_G_E26Q` | ≈1.386 kT | ln(4) |

---

## Model layers

- **A/C-only layer** — G and T are sterically excluded (`E = ∞`).
  Error rate is the fraction of adjacent identical A/C bases, with the
  closed form

  ```
  ε(Δ) = 2·exp(−Δ) / (1 + exp(−Δ))²  =  ½·sech²(Δ/2)
  ```

- **G/T-extended layer** — applies to **E26Q only**. At the A-selecting
  state, A and G compete under a two-state approximation; C and T are
  either sterically disfavored or absorbed into an effective barrier.
  This approximation is valid when the A/C discrimination barrier at
  S_A is large compared to γ_G. The E26A mutant does **not** enter this
  layer: the Gly248 steric gate still excludes G and T, so
  P(G | S_A) = 0 for E26A.

---

## Entropy metrics (analytical)

Computed from the merged-phase dinucleotide distribution
`P(AA) = P(CC) = ε/2` and `P(AC) = P(CA) = (1 − ε)/2`, giving
`H₂ = 1 + h₂(ε)` and `H(X₂|X₁) = h₂(ε)`.

| Mutant | Δ (kT) | H₁ (bit) | H₂ (bit) | H(X₂\|X₁) (bit) |
|---|---:|---:|---:|---:|
| Wild-type | → ∞ | 1.000 | 1.000 | 0.000 |
| E26A (calibrated) | 2.47 | 1.000 | 1.594 | 0.594 |
| R253A (prediction) | 2.47 | 1.000 | 1.594 | 0.594 |
| E26A_R253A (prediction) | 1.0 | 1.000 | 1.967 | 0.967 |
| Random (4-base) | 0.0 | 2.000 | 4.000 | 2.000 |

---

## Notes on tabulated values

- **E26A**: analytical value is ε = 0.1440 (from Δ calibration); the
  reported standard deviation 0.0244 is the finite-length stochastic
  ensemble spread over 50 independent sequences, not experimental or
  parameter uncertainty.
- **R253A**: analytical prediction is ε = 0.1440 under the symmetry
  assumption; the finite-length ensemble gives ε ≈ 0.141 with standard
  error ≈ 0.024. The manuscript reports the analytical value.
- **E26A_R253A**: analytical value at Δ = 1.0 is ε = 0.3930; the
  finite-length ensemble gives 0.3956 ± 0.0042.
- **Entropy values**: all analytical, computed from `H₂ = 1 + h₂(ε)`
  and `H(X₂|X₁) = h₂(ε)`. The random 4-base limit gives marginal
  Shannon entropy `H₁ = 2.000 bit`.
- **Temperature scan** (`drt3b_temperature_scan.py`): an exploratory
  numerical check at fixed Δ (kT units). Because energies are expressed
  in kT units, this scan does not constitute an independent test of
  temperature dependence; a physical temperature study would require Δ
  expressed in kcal/mol.

---

## Quick start

```bash
git clone https://github.com/hkaiopen/protein-templated-dna-sim.git
cd protein-templated-dna-sim/Drt3b

# One-command validation of all manuscript numerical claims
python validate_manuscript.py

# Individual scripts
python drt3b_stochastic_simulator.py       # Table 1: mean ± std per mutant
python drt3b_delta_scan.py                 # Table 2 / Fig. 2: ε(Δ) curve
python drt3b_gt_extension.py               # Table 3: E26Q dG fraction
python drt3b_dinucleotide_complexity.py    # Table 4: entropy metrics
python drt3b_two_state_simulator.py        # Quick single-sequence check
python drt3b_temperature_scan.py           # Exploratory ε–T scan
```

`validate_manuscript.py` prints PASS/FAIL for each numerical claim and
exits 0 only if all checks pass. Expected output ends with:

```
ALL MANUSCRIPT VALUES VERIFIED
```

---

## Reproducibility

- Ensemble statistics use 50 independent seeds (`1000`–`1049`).
- The single-sequence demonstration uses seed `42`.
- All entropy values reported in the manuscript are analytical; the
  finite-sample simulation gives `H₁ = 0.9996 ± 0.0004` for
  A/C-symmetric systems, consistent with the analytical value
  `H₁ = 1.000`.

---

## Predictions not implemented in code

Predictions involving a single-state automaton (poly(A)/poly(C)),
asymmetric barriers (AA/CC ratio = `exp(δ_low − δ_high)`), and R253A
as an open prediction are derived analytically in the manuscript and
can be evaluated in closed form. They are not part of this simulation
suite.

---

## Data and code availability

The repository provides the computational implementation of the model:
analytical evaluators for the closed-form master curve and entropy
metrics, stochastic sequence simulators for the A/C-only and
G/T-extended layers, and a Δ-scan script. A one-command validation
script (`validate_manuscript.py`) checks all numerical claims reported
in the manuscript.

The code is intended as a **computational companion** to the theoretical
model; experimental values (e.g., the 10.16% product-level dG fraction,
the E26A error rate ε = 0.144) are not independently regenerated from
raw experimental data.

---

## License

Code is released for academic and non-commercial use. Please cite the
manuscript if you use this code in published work.

## Contact

Kai Huang — hkaiopen@foxmail.com
