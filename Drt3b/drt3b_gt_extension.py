"""
G/T-extended Drt3b model — applies to **E26Q only**.

At the A-selecting state, A and G compete under the two-state
approximation stated in the manuscript (Sec. 3.2). The E26A mutant
does **not** enter this layer: the Gly248 steric gate still excludes
G and T, so P(G | S_A) = 0 for E26A.
"""

from config import (
    GAMMA_G_E26Q,
    P_G_AT_A_STATE_E26Q,
    P_G_AT_A_STATE_E26A,
    DELTA_E26A,
    q_from_delta,
)

import random


def simulate_e26q(length=400, seed=42, p_g=None):
    """Simulate the E26Q product under the two-state approximation at S_A."""
    if p_g is None:
        p_g = P_G_AT_A_STATE_E26Q
    rng = random.Random(seed)
    state_is_A = True
    seq = []
    for _ in range(length):
        if state_is_A:
            seq.append("G" if rng.random() < p_g else "A")
        else:
            seq.append("C")
        state_is_A = not state_is_A
    return "".join(seq)


def product_level_dG_fraction(seq):
    """Fraction of dG in the product (all positions)."""
    if not seq:
        return 0.0
    return seq.count("G") / len(seq)


def main():
    print("G/T-extended Drt3b model (two-state A/G approximation at S_A)")
    print("-" * 68)
    print(f"gamma_G (E26Q) = {GAMMA_G_E26Q:.4f} kT")
    print(f"P(G|S_A)       = {P_G_AT_A_STATE_E26Q:.4f}")
    print()

    # E26Q: product-level dG fraction
    seq_q = simulate_e26q()
    f_g_q = product_level_dG_fraction(seq_q)
    print(f"E26Q: simulated product-level dG fraction = {f_g_q:.4f}")

    # E26A: Gly248 still excludes G/T
    seq_a = "".join(
        "A" if i % 2 == 0 else "C" for i in range(400)
    )
    f_g_a = product_level_dG_fraction(seq_a)
    print(f"E26A: product-level dG fraction = {f_g_a:.4f}  "
          f"(P(G|S_A) = {P_G_AT_A_STATE_E26A:.2f})")


if __name__ == "__main__":
    main()
