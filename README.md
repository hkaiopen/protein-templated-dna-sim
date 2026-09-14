# Protein-Templated DNA Synthesis — Information-Dynamics Simulations

The model treats the Drt3b active pocket as a deterministic two-state
automaton in **constraint space**, with nucleotide emission governed by
a Boltzmann selection rule. All energies are in **kT units (kT = 1)**.

## Repository structure

| File | Purpose |
|---|---|
| `config.py` | Energy barriers, seeds, analytical helpers |
| `drt3b_deterministic_simulator.py` | Single-sequence simulator (fixed seed 42) |
| `drt3b_stochastic_simulator.py` | 50-replicate ensemble for mean ± std |
| `drt3b_gt_extension.py` | G/T-extended model (E26Q only; two-state A/G approximation) |
| `drt3b_dinucleotide_complexity.py` | Analytical entropy metrics |
| `drt3b_delta_scan.py` | Reproduces Fig. 2 (ε vs Δ at T = 1) |
| `drt3b_temperature_scan.py` | ε vs T at fixed Δ (supplementary check) |

## Key parameter values

| Parameter | Value | Basis |
|---|---|---|
| `DELTA_E26A` | **2.47 kT** | Calibrated from experimental ε = 0.144 |
| `DELTA_R253A` | 2.47 kT | Symmetry with E26A (model prediction) |
| `DELTA_DOUBLE` | 1.0 kT | Working assumption |
| `DELTA_RANDOM` | 0.0 kT | Thermodynamic limit |
| `P_G_AT_A_STATE_E26Q` | 0.20 | Structural ~80/20 dA/dG preference |
| `P_G_AT_A_STATE_E26A` | **0.0** | Gly248 gate still excludes G/T |
| `GAMMA_G_E26Q` | ≈1.386 kT | `ln(4)` |

## Model layers

- **A/C-only layer**: G and T excluded (`E = ∞`). Error rate is the
  fraction of adjacent identical A/C bases. Closed form:
  `ε(Δ) = 2·exp(-Δ)/(1+exp(-Δ))²`.
- **G/T-extended layer**: applies to **E26Q only**. At the A-selecting
  state, A and G compete under the two-state approximation described in
  the manuscript (Sec. 3.2). E26A does not enter this layer because the
  Gly248 steric gate still excludes G and T in the E26A mutant.

## Entropy metrics (analytical)

Computed from the merged-phase dinucleotide distribution
`P(AA) = P(CC) = ε/2`, `P(AC) = P(CA) = (1 - ε)/2`, giving
`H₂ = 1 + h₂(ε)` and `H(X₂|X₁) = h₂(ε)`.

| Mutant | Δ (kT) | H₁ (bit) | H₂ (bit) | H(X₂\|X₁) (bit) |
|---|---:|---:|---:|---:|
| Wild-type | → ∞ | 1.000 | 1.000 | 0.000 |
| E26A (calibrated) | 2.47 | 1.000 | 1.594 | 0.594 |
| R253A (prediction) | 2.47 | 1.000 | 1.594 | 0.594 |
| E26A_R253A (prediction) | 1.0 | 1.000 | 1.967 | 0.967 |
| Random (4-base) | 0.0 | 2.000 | 4.000 | 2.000 |

## Notes on tabulated values

- All entropy values above are **analytical**; finite-sample
  simulations give H₁ = 0.9996 ± 0.0004 for A/C-symmetric systems,
  consistent with the analytical value H₁ = 1.000.
- The R253A row in Table 1 of the manuscript is the analytical exact
  value under the symmetry assumption, with no standard deviation
  reported. Ensemble simulations give a standard error of ~0.0024.
- The E26A_R253A row in Table 1 of the manuscript reports the
  analytical value 0.3930. Ensemble simulations give 0.3956 ± 0.0042.

## Predictions not implemented in code

Predictions A (single-state automaton, poly(A)/poly(C)), C (asymmetric
barriers, AA/CC ratio = exp(δ_low − δ_high)), and E (R253A as open
prediction) are derived analytically in the manuscript and can be
evaluated in closed form.

## Reproducing the tables and figures

```bash
python drt3b_stochastic_simulator.py       # Table 1
python drt3b_delta_scan.py                 # Table 2 and Fig. 2
python drt3b_gt_extension.py               # Table 3
python drt3b_dinucleotide_complexity.py    # Table 4
```

Ensemble statistics use 50 independent seeds (1000–1049); the
single-sequence demonstration uses seed 42.

## License

Code is released for academic and non-commercial use. Please cite the
manuscript if you use this code in published work.

## Contact

Kai Huang — hkaiopen@foxmail.com
