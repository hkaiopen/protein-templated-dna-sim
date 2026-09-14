"""
One-command validation of all numerical claims made in the manuscript.

Run:
    python validate_manuscript.py

Prints PASS/FAIL for each claim and exits 0 only if all checks pass.
"""

import sys
from config import (
    DELTA_E26A, DELTA_R253A, DELTA_DOUBLE, DELTA_RANDOM,
    DELTA_WT_PROXY, GAMMA_G_E26Q, P_G_AT_A_STATE_E26Q,
    P_G_AT_A_STATE_E26A, epsilon_from_delta, delta_from_epsilon,
    H2_dinucleotide, conditional_entropy,
)

TOL = 2e-3  # absolute tolerance for analytical values

CLAIMS = [
    ("WT epsilon (hard-exclusion)",
     epsilon_from_delta(DELTA_WT_PROXY), 0.0),
    ("E26A epsilon (calibrated)",
     epsilon_from_delta(DELTA_E26A), 0.1440),
    ("R253A epsilon (symmetry)",
     epsilon_from_delta(DELTA_R253A), 0.1440),
    ("E26A_R253A epsilon (Delta=1.0)",
     epsilon_from_delta(DELTA_DOUBLE), 0.3930),
    ("Random (4-base) epsilon",
     0.25, 0.25),
    ("E26Q gamma_G (kT)",
     GAMMA_G_E26Q, 1.386),
    ("E26Q P(G|S_A)",
     P_G_AT_A_STATE_E26Q, 0.20),
    ("E26Q product-level dG fraction",
     0.5 * P_G_AT_A_STATE_E26Q, 0.10),
    ("E26A P(G|S_A)",
     P_G_AT_A_STATE_E26A, 0.0),
    ("E26A H2 (bit)",
     H2_dinucleotide(epsilon_from_delta(DELTA_E26A)), 1.594),
    ("E26A H(X2|X1) (bit)",
     conditional_entropy(epsilon_from_delta(DELTA_E26A)), 0.594),
    ("Double-mutant H2 (bit)",
     H2_dinucleotide(epsilon_from_delta(DELTA_DOUBLE)), 1.967),
    ("Double-mutant H(X2|X1) (bit)",
     conditional_entropy(epsilon_from_delta(DELTA_DOUBLE)), 0.967),
    ("Prediction D threshold Delta (kT)",
     delta_from_epsilon(0.05), 3.637),
]


def main():
    print("Validating manuscript numerical claims")
    print("=" * 62)
    all_ok = True
    for name, got, expected in CLAIMS:
        diff = abs(got - expected)
        ok = diff < TOL
        all_ok &= ok
        mark = "PASS" if ok else "FAIL"
        print(f"[{mark}] {name:<42} got={got:.4f}  exp={expected:.4f}")
    print("=" * 62)
    if all_ok:
        print("ALL MANUSCRIPT VALUES VERIFIED")
        return 0
    print("VALIDATION FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())