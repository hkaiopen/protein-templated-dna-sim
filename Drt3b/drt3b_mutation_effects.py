#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
drt3b_mutation_effects.py
=========================

Mutation effects on Drt3b protein-templated DNA synthesis under the
information-dynamics model (constraint space + Boltzmann selection).

Data sources: Deng et al., Science 392, 1274 (2026), Figs. 4H, S11A,
              S11D, S11F.
Note:         R253A has not been experimentally characterized.

"""

import numpy as np
from collections import Counter


# ------------------------------------------------------------------
# Shared constants (must match dinucleotide_complexity.py)
# ------------------------------------------------------------------
RNG_SEEDS  = list(range(1000, 1050))   # 50 independent seeds
SEQ_LENGTH = 500
DEMO_SEED  = 42

# Δ_AC defaults (A/C-only layer)
DELTA_AC_DEFAULTS = {
    'wildtype'   : 100.0,
    'E26A'       : 2.5,
    'E26Q'       : 2.5,     # treated as symmetric to E26A in Layer A
    'R253A'      : 2.5,
    'E26A_R253A' : 1.0,
}

# Δ_G,A (G/T-extended layer); p_G = 0.20 → Δ ≈ 1.3863 kT
DELTA_GA_DEFAULTS = {
    'E26A' : float(np.log(0.8 / 0.2)),
    'E26Q' : float(np.log(0.8 / 0.2)),
}


# ==================================================================
# 1. Calibration utilities
# ==================================================================

def calibrate_delta_AC_from_AC6(f_AC6):
    """
    Return the A/C mismatch barrier Δ_AC consistent with an experimental
    (AC)6 motif frequency f_AC6 (both phases). Small-p closed form:
        p ≈ 1 − f_AC6^(1/6);  Δ_AC = −ln(p/(1−p))
    For large p, f_AC6 is dominated by q^6 + p^6 and the approximation
    loses accuracy; use the bisection-based inverse for precision.
    """
    if f_AC6 <= 0 or f_AC6 > 1:
        raise ValueError("f_AC6 must lie in (0, 1].")
    p = 1.0 - f_AC6 ** (1.0 / 6.0)
    if p <= 0:
        return float('inf')
    return -np.log(p / (1.0 - p))


def calibrate_delta_G_A_from_A_state_fraction(p_G):
    """Δ_G,A = ln[(1 − p_G) / p_G] from A-state G-incorporation prob."""
    if p_G <= 0:
        return float('inf')
    if p_G >= 1:
        return 0.0
    return float(np.log((1.0 - p_G) / p_G))


def calibrate_delta_G_A_from_product_fraction(f_G_product):
    """
    Δ_G,A from product-level dG fraction, assuming A-state occupies half
    of the cycle (so per-step A-state probability p_G = 2 f_G_product).
    """
    return calibrate_delta_G_A_from_A_state_fraction(2.0 * f_G_product)


# ==================================================================
# 2. Simulators (one per layer)
# ==================================================================

def simulate_ac_only(mutation, rng, L=SEQ_LENGTH):
    """Layer (A): A/C-only two-state automaton with mismatch barrier Δ."""
    delta = DELTA_AC_DEFAULTS[mutation]
    p_wrong = np.exp(-delta) / (1.0 + np.exp(-delta))
    state = 'S_A' if rng.random() < 0.5 else 'S_C'
    seq = []
    for _ in range(L):
        if state == 'S_A':
            seq.append('C' if rng.random() < p_wrong else 'A')
        else:
            seq.append('A' if rng.random() < p_wrong else 'C')
        state = 'S_C' if state == 'S_A' else 'S_A'
    return ''.join(seq)


def simulate_gt_ext(mutation, rng, L=SEQ_LENGTH):
    """Layer (B): G/T-extended. G only in S_A, T excluded everywhere."""
    delta_GA = DELTA_GA_DEFAULTS[mutation]
    p_G = np.exp(-delta_GA) / (1.0 + np.exp(-delta_GA))
    state = 'S_A' if rng.random() < 0.5 else 'S_C'
    seq = []
    for _ in range(L):
        if state == 'S_A':
            seq.append('G' if rng.random() < p_G else 'A')
        else:
            seq.append('C')
        state = 'S_C' if state == 'S_A' else 'S_A'
    return ''.join(seq)


def simulate_random_gt(rng, L=SEQ_LENGTH):
    """Fully random model with G/T allowed (all four bases uniform)."""
    state = 'S_A' if rng.random() < 0.5 else 'S_C'
    seq = []
    for _ in range(L):
        seq.append(rng.choice(['A', 'C', 'G', 'T']))
        state = 'S_C' if state == 'S_A' else 'S_A'
    return ''.join(seq)


# ==================================================================
# 3. Metrics
# ==================================================================

def error_rate_adjacent(seq):
    """Fraction of adjacent identical bases (Layer A metric)."""
    if len(seq) < 2:
        return 0.0
    return sum(1 for i in range(len(seq) - 1)
               if seq[i] == seq[i + 1]) / (len(seq) - 1)


def dG_fraction(seq):
    """Product-level dG fraction (Layer B metric)."""
    return seq.count('G') / len(seq) if seq else 0.0


def composition(seq):
    N = len(seq)
    return {b: seq.count(b) / N for b in 'ACGT'} if N else {b: 0.0 for b in 'ACGT'}


# ==================================================================
# 4. Ensemble runners
# ==================================================================

def _stats(arr):
    arr = np.asarray(arr, dtype=float)
    return float(arr.mean()), float(arr.std(ddof=1) / np.sqrt(len(arr)))


def run_ensemble_ac(mutation):
    eps, comp = [], []
    for s in RNG_SEEDS:
        seq = simulate_ac_only(mutation, np.random.default_rng(s))
        eps.append(error_rate_adjacent(seq))
        comp.append(composition(seq))
    eps_m, eps_se = _stats(eps)
    comp_mean = {b: float(np.mean([c[b] for c in comp])) for b in 'ACGT'}
    return {'eps_mean': eps_m, 'eps_se': eps_se, 'composition': comp_mean}


def run_ensemble_gt(mutation):
    dG, comp = [], []
    for s in RNG_SEEDS:
        seq = simulate_gt_ext(mutation, np.random.default_rng(s))
        dG.append(dG_fraction(seq))
        comp.append(composition(seq))
    dG_m, dG_se = _stats(dG)
    comp_mean = {b: float(np.mean([c[b] for c in comp])) for b in 'ACGT'}
    return {'dG_mean': dG_m, 'dG_se': dG_se, 'composition': comp_mean}


def run_ensemble_random_gt():
    eps, dG, comp = [], [], []
    for s in RNG_SEEDS:
        seq = simulate_random_gt(np.random.default_rng(s))
        eps.append(error_rate_adjacent(seq))
        dG.append(dG_fraction(seq))
        comp.append(composition(seq))
    eps_m, eps_se = _stats(eps)
    dG_m, dG_se   = _stats(dG)
    comp_mean = {b: float(np.mean([c[b] for c in comp])) for b in 'ACGT'}
    return {'eps_mean': eps_m, 'eps_se': eps_se,
            'dG_mean': dG_m, 'dG_se': dG_se,
            'composition': comp_mean}


# ==================================================================
# 5. Main
# ==================================================================

def main():
    sep = '=' * 94
    print(sep)
    print("drt3b_mutation_effects.py — Drt3b information-dynamics model")
    print(f"50 seeds × {SEQ_LENGTH} bp | empirical SEs (Reviewer 1, P1)")
    print("Data: Deng et al., Science 392, 1274 (2026), Figs. 4H, S11A, "
          "S11D, S11F")
    print("Note: R253A has not been experimentally characterized.")
    print(sep)

    # --------------------------------------------------------------
    # Layer (A) — A/C-only sub-model
    # --------------------------------------------------------------
    print("\n[Layer (A)] A/C-only sub-model (E = ∞ for G, T)")
    print("            Metric ε = fraction of adjacent identical bases")
    print("            Closed form: ε(Δ) = 2 e^(−Δ) / (1 + e^(−Δ))^2")
    print("-" * 94)
    print(f"{'Mutation':<14} {'Δ_AC':>6} {'ε (mean ± SE)':>22} "
          f"{'ε analytic':>12} {'|Δε|':>8} "
          f"{'A':>6} {'C':>6} {'G':>6} {'T':>6}")
    print("-" * 94)
    ac_muts = ['wildtype', 'E26A', 'E26Q', 'R253A', 'E26A_R253A']
    res_ac = {}
    for mut in ac_muts:
        r = run_ensemble_ac(mut)
        res_ac[mut] = r
        d = DELTA_AC_DEFAULTS[mut]
        eps_a = 2.0 * np.exp(-d) / (1.0 + np.exp(-d)) ** 2
        comp = r['composition']
        print(f"{mut:<14} {d:>6.2f} "
              f"{r['eps_mean']:>10.4f} ± {r['eps_se']:>6.4f} "
              f"{eps_a:>12.4f} {abs(r['eps_mean'] - eps_a):>8.4f} "
              f"{comp['A']:>6.3f} {comp['C']:>6.3f} "
              f"{comp['G']:>6.3f} {comp['T']:>6.3f}")

    # Random (G/T allowed) — Layer A-adjacent metric but 4-letter alphabet
    r_rand = run_ensemble_random_gt()
    comp = r_rand['composition']
    print(f"{'random(G/T)':<14} {'—':>6} "
          f"{r_rand['eps_mean']:>10.4f} ± {r_rand['eps_se']:>6.4f} "
          f"{'—':>12} {'—':>8} "
          f"{comp['A']:>6.3f} {comp['C']:>6.3f} "
          f"{comp['G']:>6.3f} {comp['T']:>6.3f}")

    # --------------------------------------------------------------
    # Layer (B) — G/T-extended sub-model
    # --------------------------------------------------------------
    print("\n[Layer (B)] G/T-extended model — G only in S_A; T excluded")
    print(f"            Δ_G,A = ln[(1 − p_G)/p_G] = "
          f"{DELTA_GA_DEFAULTS['E26Q']:.4f} kT  (p_G = 0.20)")
    print("            Metric = product-level dG fraction = 0.5 · p_G")
    print("            NOTE: adjacent-identical ε is identically zero in "
          "this layer")
    print("-" * 94)
    print(f"{'Mutation':<14} {'Δ_G,A':>8} {'dG frac ± SE':>20} "
          f"{'A':>6} {'C':>6} {'G':>6} {'T':>6}")
    print("-" * 94)
    gt_muts = ['E26A', 'E26Q']
    res_gt = {}
    for mut in gt_muts:
        r = run_ensemble_gt(mut)
        res_gt[mut] = r
        dGA = DELTA_GA_DEFAULTS[mut]
        comp = r['composition']
        print(f"{mut:<14} {dGA:>8.4f} "
              f"{r['dG_mean']:>10.4f} ± {r['dG_se']:>6.4f} "
              f"{comp['A']:>6.3f} {comp['C']:>6.3f} "
              f"{comp['G']:>6.3f} {comp['T']:>6.3f}")
    print("-" * 94)
    print("  Experimental anchor: E26Q dG fraction ≈ 10.16% "
          "(Deng et al., Fig. S11D/S11F).")
    print("  Model reproduces this by construction (calibration, not "
          "independent prediction).")
    print("  E26A shares the same per-step Δ_G,A; its lower overall "
          "yield is")
    print("  attributed to a processivity defect (Fig. 4H), outside "
          "this model's scope.")

    # --------------------------------------------------------------
    # Representative sequences
    # --------------------------------------------------------------
    print("\n" + "=" * 94)
    print(f"Representative sequences (first trial, seed {DEMO_SEED}, "
          "80 nt)")
    print("=" * 94)
    print("\n  Layer (A) — A/C-only sub-model:")
    for mut in ac_muts:
        seq = simulate_ac_only(mut, np.random.default_rng(DEMO_SEED), L=80)
        print(f"    {mut:<14} : {seq}")
    print("\n  Layer (B) — G/T-extended sub-model:")
    for mut in gt_muts:
        seq = simulate_gt_ext(mut, np.random.default_rng(DEMO_SEED), L=80)
        print(f"    {mut:<14} : {seq}")

    # --------------------------------------------------------------
    # Cross-check against the reviewer response letter
    # --------------------------------------------------------------
    print("\n" + "=" * 94)
    print("Cross-check against the reviewer response letter")
    print("=" * 94)
    print(f"  E26A Layer A ε(2.5)      = "
          f"{res_ac['E26A']['eps_mean']:.4f}  "
          f"(resp. 0.144, closed form 0.1402)")
    print(f"  E26A_R253A Δ=1.0 ε(1.0)  = "
          f"{res_ac['E26A_R253A']['eps_mean']:.4f}  "
          f"(resp. 0.393, closed form 0.3932)")
    print(f"  Δ=3.0 closed form ε(3.0) = "
          f"{2*np.exp(-3.0)/(1+np.exp(-3.0))**2:.4f}  "
          f"(resp. 0.090)")
    print(f"  Random (G/T allowed) ε   = "
          f"{r_rand['eps_mean']:.4f}  "
          f"(manuscript 0.244)")
    print(f"  E26Q Layer B dG fraction = "
          f"{res_gt['E26Q']['dG_mean']:.4f}  "
          f"(experimental anchor 0.1016)")

    print("\n" + sep)
    print("DONE. Two layers kept separate; no mixed-model numbers.")
    print(sep)


if __name__ == "__main__":
    main()
