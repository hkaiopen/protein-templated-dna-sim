#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
drt3b_noise_robustness.py
=========================

Noise robustness of the Drt3b information-dynamics model.
Reproduces manuscript Figure 2: error rate as a function of temperature
for two representative A/C mismatch barriers (wild-type Δ = 100 kT and
the predicted double mutant E26A_R253A Δ = 1.0 kT).

Model
-----
Two-state automaton S_A ↔ S_C with deterministic alternation and
Boltzmann nucleotide selection at finite temperature T. The effective
per-step mismatch probability is

    p_wrong(Δ, T) = e^(−Δ/T) / (1 + e^(−Δ/T)),

giving the exact closed-form error rate for the A/C-only layer

    ε(Δ, T) = 2 e^(−Δ/T) / (1 + e^(−Δ/T))^2.

The Monte Carlo simulation must reproduce this closed form to within
the empirical standard error.

Aligned with
------------
- drt3b_mutation_effects.py    (same constants, same empirical SE)
- drt3b_dinucleotide_complexity.py (same metric definition)
- Reviewer 1, Point 1 : empirical SE from 50 independent replicas
- Reviewer 1, Point 4 : Δ = 100 kT is a numerical proxy for the
                        hard-exclusion limit; ε < 0.001 for any
                        Δ ≳ 10 at T = 1.
- Reviewer 2, Point 3 : terminology "constraint space".

Authors: Kai Huang, Hongkui Liu, Ziwei Huang
Repository: https://github.com/hkaiopen/protein-templated-dna-sim
"""

import numpy as np


# ------------------------------------------------------------------
# Shared constants (must match drt3b_mutation_effects.py)
# ------------------------------------------------------------------
RNG_SEEDS  = list(range(1000, 1050))   # 50 independent seeds
SEQ_LENGTH = 500

# Representative barriers for Figure 2
DELTA_WT       = 100.0   # numerical proxy for the hard-exclusion limit
DELTA_DOUBLE   = 1.0     # predicted double mutant E26A_R253A


# ==================================================================
# 1. Simulator (A/C-only layer, temperature-scaled Boltzmann)
# ==================================================================

def simulate_ac_at_T(delta, T, L, rng):
    """
    A/C-only two-state automaton at temperature T (kT units).
    Per-step mismatch probability:
        p_wrong = e^(−Δ/T) / (1 + e^(−Δ/T)).
    """
    x = delta / T
    p_wrong = np.exp(-x) / (1.0 + np.exp(-x))
    state = 'S_A' if rng.random() < 0.5 else 'S_C'
    seq = []
    for _ in range(L):
        if state == 'S_A':
            seq.append('C' if rng.random() < p_wrong else 'A')
        else:
            seq.append('A' if rng.random() < p_wrong else 'C')
        state = 'S_C' if state == 'S_A' else 'S_A'
    return ''.join(seq)


# ==================================================================
# 2. Metrics and closed form
# ==================================================================

def error_rate_adjacent(seq):
    """Fraction of adjacent identical bases (Layer A metric)."""
    if len(seq) < 2:
        return 0.0
    return sum(1 for i in range(len(seq) - 1)
               if seq[i] == seq[i + 1]) / (len(seq) - 1)


def closed_form_eps(delta, T):
    """Exact ε(Δ, T) = 2 e^(−Δ/T) / (1 + e^(−Δ/T))^2."""
    x = delta / T
    return 2.0 * np.exp(-x) / (1.0 + np.exp(-x)) ** 2


# ==================================================================
# 3. Ensemble runner (empirical SE, Reviewer 1 P1)
# ==================================================================

def run_ensemble(delta, T):
    eps = []
    for s in RNG_SEEDS:
        seq = simulate_ac_at_T(delta, T, SEQ_LENGTH,
                               np.random.default_rng(s))
        eps.append(error_rate_adjacent(seq))
    eps = np.asarray(eps, dtype=float)
    mean = float(eps.mean())
    se   = float(eps.std(ddof=1) / np.sqrt(len(eps)))
    return mean, se


# ==================================================================
# 4. Main driver — manuscript Figure 2
# ==================================================================

def main():
    sep = '=' * 82
    print(sep)
    print("drt3b_noise_robustness.py — Drt3b information-dynamics model")
    print(f"50 seeds × {SEQ_LENGTH} bp | empirical SEs (Reviewer 1, P1)")
    print("Reproduces manuscript Figure 2")
    print(sep)

    temps = [0.1, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    cases = [
        ("Wild-type (Δ = 100)",     DELTA_WT),
        ("E26A_R253A (Δ = 1.0)",    DELTA_DOUBLE),
    ]

    for label, delta in cases:
        print(f"\n[{label}]")
        print(f"  {'T (kT)':>8}  {'ε (mean ± SE)':>22}  "
              f"{'ε analytic':>12}  {'|Δε|':>8}")
        print("  " + "-" * 62)
        for T in temps:
            m, se = run_ensemble(delta, T)
            ana   = closed_form_eps(delta, T)
            print(f"  {T:>8.2f}  {m:>10.4f} ± {se:>6.4f}  "
                  f"{ana:>12.4f}  {abs(m - ana):>8.4f}")

    # ------------------------------------------------
    # Summary for manuscript Figure 2
    # ------------------------------------------------
    print()
    print(sep)
    print("Figure 2 data (for direct plotting):")
    print(sep)
    print(f"  {'T (kT)':>8}  {'ε_WT(Δ=100)':>14}  {'ε_pred(Δ=1.0)':>16}")
    print("  " + "-" * 44)
    for T in temps:
        e_wt  = closed_form_eps(DELTA_WT, T)
        e_dbl = closed_form_eps(DELTA_DOUBLE, T)
        print(f"  {T:>8.2f}  {e_wt:>14.4f}  {e_dbl:>16.4f}")

    print()
    print("Interpretation:")
    print("  * Wild-type Δ = 100: error rate remains < 10⁻³ for all T ≤ 10.")
    print("  * Δ = 1 (model prediction): error rate rises with T; at T = 1")
    print("    it equals ε(1) = 0.3932 and at T = 10 it approaches the")
    print("    A/C-only random limit ε → 0.50.")
    print("  * Δ = 100 is a numerical proxy for the hard-exclusion limit;")
    print("    ε < 0.001 for any Δ ≳ 10 at T = 1 (Reviewer 1, P4).")
    print(sep)


if __name__ == "__main__":
    main()
