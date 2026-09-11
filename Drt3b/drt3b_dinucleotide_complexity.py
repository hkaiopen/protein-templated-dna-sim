#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dinucleotide_complexity.py
==========================

Dinucleotide frequencies, entropy measures, and sequence complexity
for the information-dynamics model of Drt3b protein-templated DNA synthesis.

Model layers (labelled in every output table):

  Layer (i)  A/C-only sub-model (E = ∞ for G, T)
      ε(Δ) = 2 e^(−Δ) / (1 + e^(−Δ))^2

  Layer (ii) G/T-extended sub-model (only G allowed, only in S_A)
      Δ_G,A = ln[(1 − p_G) / p_G], p_G ≈ 0.20 → Δ_G,A ≈ 1.386 kT
      Adjacent-identical metric is identically zero here; the reported
      quantity is the product-level dG fraction = 0.5 × p_G.

"""

import numpy as np
from collections import Counter


# ------------------------------------------------------------------
# Global simulation constants (shared with drt3b_mutation_effects.py)
# ------------------------------------------------------------------
RNG_SEEDS  = list(range(1000, 1050))    # 50 independent seeds
SEQ_LENGTH = 500                        # bp per sequence
DEMO_SEED  = 42                         # reproducible single-run demo


# ==================================================================
# 1. Emission laws (Boltzmann selection on the constraint space)
# ==================================================================

def boltzmann_ac(state, delta):
    """A/C-only layer: constraint space {S_A, S_C}, mismatch penalty Δ."""
    p_wrong = np.exp(-delta) / (1.0 + np.exp(-delta))
    p_correct = 1.0 - p_wrong
    return ({'A': p_correct, 'C': p_wrong} if state == 'S_A'
            else {'C': p_correct, 'A': p_wrong})


def boltzmann_gt(state, delta_GA):
    """G/T-extended layer: G allowed only in S_A. T excluded everywhere."""
    if state == 'S_A':
        p_G = np.exp(-delta_GA) / (1.0 + np.exp(-delta_GA))
        return {'A': 1.0 - p_G, 'G': p_G}
    return {'C': 1.0}


# ==================================================================
# 2. Sequence generators
# ==================================================================

def simulate_ac(delta, L, rng):
    """A/C-only two-state automaton; initial phase uniform."""
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


def simulate_gt(delta_GA, L, rng):
    """G/T-extended two-state automaton; G only in S_A; T excluded."""
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


def simulate_three_state_aac(L, rng):
    """Deterministic 3-state cycle S_A → S_A → S_C → ... (poly(AAC))."""
    cycle = ['A', 'A', 'C']
    start = int(rng.integers(0, 3))
    return ''.join(cycle[(start + i) % 3] for i in range(L))


# ==================================================================
# 3. Entropy measures (Reviewer 1, Point 5)
# ==================================================================

def marginal_entropy(seq):
    """H1 = −Σ p(n) log2 p(n) over single nucleotides."""
    N = len(seq)
    counts = Counter(seq)
    return -sum((c / N) * np.log2(c / N) for c in counts.values())


def dinucleotide_entropy(seq):
    """H2 = −Σ p(n1n2) log2 p(n1n2) over adjacent pairs."""
    pairs = [seq[i:i + 2] for i in range(len(seq) - 1)]
    N = sum(1 for _ in pairs)
    counts = Counter(pairs)
    return -sum((c / N) * np.log2(c / N) for c in counts.values())


def entropy_rate(seq):
    """H(X2|X1) = H2 − H1. For deterministic alternation this is zero;
    for poly(AAC) with unknown phase it is 0.667 bit."""
    return dinucleotide_entropy(seq) - marginal_entropy(seq)


# ==================================================================
# 4. Error metric and (AC)6 mapping (Reviewer 1, Point 3)
# ==================================================================

def error_rate_adjacent(seq):
    """Fraction of adjacent identical bases (A/C-only layer)."""
    if len(seq) < 2:
        return 0.0
    matches = sum(1 for i in range(len(seq) - 1) if seq[i] == seq[i + 1])
    return matches / (len(seq) - 1)


def ac6_frequency(seq):
    """Fraction of length-6 windows that are perfect AC alternations,
    counting BOTH phases (ACACAC and CACACA)."""
    if len(seq) < 6:
        return 0.0
    hits = sum(1 for i in range(len(seq) - 5)
               if seq[i:i + 6] in ('ACACAC', 'CACACA'))
    return hits / (len(seq) - 5)


def f_ac6_from_p(p):
    """Analytic forward map for BOTH phases: f_AC6 = q^6 + p^6."""
    q = 1.0 - p
    return q ** 6 + p ** 6


def p_from_f_ac6(f):
    """Inverse map by bisection on [0, 0.5]; monotone decreasing."""
    f = max(f, 0.5 ** 6)  # f can't go below 2·(1/2)^6
    if f >= 1.0:
        return 0.0
    lo, hi = 0.0, 0.5
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if f_ac6_from_p(mid) > f:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def eps_from_f_ac6(f):
    """Bridge from experimental f_AC6 to adjacent-identical error rate."""
    p = p_from_f_ac6(f)
    return 2.0 * p * (1.0 - p)


# ==================================================================
# 5. LZ76 complexity (normalised by n / log2 n)
# ==================================================================

def lz76_complexity(seq):
    """Normalised LZ76 complexity C_LZ = c(n) / (n / log2 n)."""
    n = len(seq)
    if n <= 1:
        return 0.0
    i, c = 0, 0
    while i < n:
        length = 1
        while i + length <= n:
            if seq[i:i + length] not in seq[:i]:
                break
            length += 1
        c += 1
        i += length
    return c / (n / np.log2(n))


# ==================================================================
# 6. Ensemble analysis with empirical standard errors
# ==================================================================

def analyze_ensemble(seqs):
    def _stats(arr):
        arr = np.asarray(arr, dtype=float)
        return float(arr.mean()), float(arr.std(ddof=1) / np.sqrt(len(arr)))
    return {
        'error_rate': _stats([error_rate_adjacent(s) for s in seqs]),
        'H1':         _stats([marginal_entropy(s)     for s in seqs]),
        'H2':         _stats([dinucleotide_entropy(s) for s in seqs]),
        'H_rate':     _stats([entropy_rate(s)         for s in seqs]),
        'f_AC6':      _stats([ac6_frequency(s)        for s in seqs]),
        'LZ':         _stats([lz76_complexity(s)      for s in seqs]),
        'n_runs':     len(seqs),
    }


def run_ensemble(simulator, **kwargs):
    seqs = [simulator(rng=np.random.default_rng(s), L=SEQ_LENGTH, **kwargs)
            for s in RNG_SEEDS]
    return analyze_ensemble(seqs)


# ==================================================================
# 7. Model registry
# ==================================================================

AC_MODELS = {
    'Wild-type (Δ=100)'          : 100.0,
    'E26A (Δ=2.5)'               : 2.5,
    'R253A (Δ=2.5, pred.)'       : 2.5,
    'E26A_R253A (Δ=1.0, pred.)'  : 1.0,
    'Random (Δ=0, A/C-only)'     : 0.0,
}

DELTA_GA_CALIBRATED = float(np.log(0.8 / 0.2))   # ≈ 1.3863 kT

GT_MODELS = {
    'E26A (Δ_G,A calib.)': DELTA_GA_CALIBRATED,
    'E26Q (Δ_G,A calib.)': DELTA_GA_CALIBRATED,
}


# ==================================================================
# 8. Main driver
# ==================================================================

def main():
    sep = '=' * 82
    print(sep)
    print("dinucleotide_complexity.py — Drt3b information-dynamics model")
    print(f"50 seeds × {SEQ_LENGTH} bp | empirical SEs (Reviewer 1, P1)")
    print(sep)

    # --------------------------------------------------------------
    # Part 1 — A/C-only sub-model
    # --------------------------------------------------------------
    print("\n[Part 1] A/C-only sub-model (G, T excluded; E = ∞)")
    print("         Metric ε = fraction of adjacent identical bases")
    print("         LZ = LZ76 phrases / (n / log2 n)")
    print("-" * 82)
    print(f"{'Mutant':<28} {'ε':>10} {'H1':>9} {'H2':>9} "
          f"{'H(X2|X1)':>10} {'f_AC6':>9} {'LZ':>8}")
    print("-" * 82)

    ac_results = {}
    for name, delta in AC_MODELS.items():
        res = run_ensemble(simulate_ac, delta=delta)
        ac_results[name] = res
        print(f"{name:<28} "
              f"{res['error_rate'][0]:>10.4f} "
              f"{res['H1'][0]:>9.4f} "
              f"{res['H2'][0]:>9.4f} "
              f"{res['H_rate'][0]:>10.4f} "
              f"{res['f_AC6'][0]:>9.4f} "
              f"{res['LZ'][0]:>8.4f}")

    # -- Analytic cross-check of ε(Δ) --
    print("\n  Analytic cross-check: ε(Δ) = 2 e^(−Δ) / (1 + e^(−Δ))^2")
    print(f"  {'Mutant':<28} {'Δ':>6} {'ε analytic':>12} "
          f"{'ε simulated':>14} {'Δε':>10}")
    for name, delta in AC_MODELS.items():
        eps_a = 2.0 * np.exp(-delta) / (1.0 + np.exp(-delta)) ** 2
        eps_s = ac_results[name]['error_rate'][0]
        print(f"  {name:<28} {delta:>6.2f} {eps_a:>12.4f} "
              f"{eps_s:>14.4f} {eps_s - eps_a:>+10.4f}")

    # -- (AC)6 mapping self-consistency --
    print("\n  (AC)6 mapping self-consistency (Reviewer 1, P3):")
    print("  Forward: f_AC6 = q^6 + p^6   (both phases, exact)")
    print("  Inverse: p from bisection;  ε = 2 p (1 − p)")
    print(f"  {'Mutant':<28} {'p':>8} {'f_AC6 analytic':>16} "
          f"{'f_AC6 simulated':>17} {'ε from f_AC6':>15}")
    for name, delta in AC_MODELS.items():
        p       = np.exp(-delta) / (1.0 + np.exp(-delta))
        f_an    = f_ac6_from_p(p)
        f_sim   = ac_results[name]['f_AC6'][0]
        eps_map = eps_from_f_ac6(f_sim)
        print(f"  {name:<28} {p:>8.4f} {f_an:>16.4f} "
              f"{f_sim:>17.4f} {eps_map:>15.4f}")

    # --------------------------------------------------------------
    # Part 2 — G/T-extended model (E26Q, E26A)
    # --------------------------------------------------------------
    print("\n[Part 2] G/T-extended model — G allowed only in S_A; "
          "T excluded")
    print(f"         Δ_G,A = ln[(1 − p_G)/p_G] = "
          f"{DELTA_GA_CALIBRATED:.4f} kT  (p_G = 0.20)")
    print("-" * 82)
    print(f"{'Mutant':<28} {'f_G':>9} {'dG frac':>10} {'H1':>9} "
          f"{'H2':>9} {'H(X2|X1)':>10}")
    for name, dGA in GT_MODELS.items():
        seqs = [simulate_gt(dGA, SEQ_LENGTH, np.random.default_rng(s))
                for s in RNG_SEEDS]
        f_G   = float(np.mean([s.count('G') / len(s) for s in seqs]))
        h1    = float(np.mean([marginal_entropy(s)     for s in seqs]))
        h2    = float(np.mean([dinucleotide_entropy(s) for s in seqs]))
        dG_prod = 0.5 * (np.exp(-dGA) / (1.0 + np.exp(-dGA)))
        print(f"{name:<28} {f_G:>9.4f} {dG_prod:>10.4f} "
              f"{h1:>9.4f} {h2:>9.4f} {h2 - h1:>10.4f}")
    print("         Note: state alternation is deterministic in this "
          "layer, so the")
    print("         adjacent-identical metric ε is identically zero "
          "and uninformative.")

    # --------------------------------------------------------------
    # Part 3 — 3-state AAC model (Reviewer 1, Point 5)
    # --------------------------------------------------------------
    print("\n[Part 3] 3-state AAC model — VERIFICATION of Reviewer 1, P5")
    print("-" * 82)
    seqs = [simulate_three_state_aac(SEQ_LENGTH, np.random.default_rng(s))
            for s in RNG_SEEDS]
    res = analyze_ensemble(seqs)
    print(f"  Deterministic poly(AAC) repeat, length = {SEQ_LENGTH} bp, "
          f"{len(seqs)} replicas")
    print(f"    Marginal entropy  H1        = {res['H1'][0]:.4f} bit "
          f"  [corrected value: 0.918 bit]")
    print(f"    Dinucleotide       H2        = {res['H2'][0]:.4f} bit "
          f"  [theory: log2(3) = 1.585 bit]")
    print(f"    Entropy rate       H(X2|X1)  = {res['H_rate'][0]:.4f} bit")
    print( "      [phase known:   0.0000 bit]")
    print( "      [phase unknown: 0.6667 bit = (2/3)*1 + (1/3)*0]")
    h1_exact = -(2/3) * np.log2(2/3) - (1/3) * np.log2(1/3)
    print(f"    Closed form H1 = −(2/3)log2(2/3) − (1/3)log2(1/3) "
          f"= {h1_exact:.4f} bit")
    print("    NOTE: 1.585 bit is the *dinucleotide* entropy H2, not "
          "the marginal H1.")
    print("          The original manuscript conflated the two.")

    # --------------------------------------------------------------
    # Part 4 — Summary cross-check
    # --------------------------------------------------------------
    print("\n[Part 4] Summary — manuscript cross-check (Reviewer 1, P1)")
    print("-" * 82)
    print(f"  {'Quantity':<38} {'Value':>12} {'Source':>22}")
    print(f"  {'-'*38} {'-'*12} {'-'*22}")
    rows = [
        ('E26A error rate ε(2.5)',
         f"{ac_results['E26A (Δ=2.5)']['error_rate'][0]:.4f}",
         'Table 2 / Resp. R1-P1'),
        ('E26A_R253A ε(1.0)',
         f"{ac_results['E26A_R253A (Δ=1.0, pred.)']['error_rate'][0]:.4f}",
         'Table 2 / Resp. R1-P1'),
        ('Random (A/C-only) ε',
         f"{ac_results['Random (Δ=0, A/C-only)']['error_rate'][0]:.4f}",
         'Table 2'),
        ('Analytic ε(1.0)',
         f"{2*np.exp(-1.0)/(1+np.exp(-1.0))**2:.4f}",
         'Closed form'),
        ('Analytic ε(3.0)',
         f"{2*np.exp(-3.0)/(1+np.exp(-3.0))**2:.4f}",
         'Closed form'),
        ('E26Q product-level dG',
         f"{0.5*(np.exp(-DELTA_GA_CALIBRATED)/(1+np.exp(-DELTA_GA_CALIBRATED))):.4f}",
         'Exp. 10.16% → calib.'),
    ]
    for label, val, src in rows:
        print(f"  {label:<38} {val:>12} {src:>22}")

    print("\n" + sep)
    print("DONE. All quantities labelled by (model layer, metric).")
    print(sep)


if __name__ == "__main__":
    main()
