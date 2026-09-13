"""
Dinucleotide entropy and complexity metrics for the A/C-only Drt3b model.

All tabulated entropy values are analytical, computed from the closed
form H_2 = 1 + h2(epsilon) and H(X2|X1) = h2(epsilon), where epsilon is
the error rate of the state (Eq. 8 in the manuscript).

The merged-phase dinucleotide distribution is P(AA) = P(CC) = eps/2 and
P(AC) = P(CA) = (1-eps)/2.
"""

from config import (
    DELTA_AC_DEFAULTS,
    epsilon_from_delta,
    H2_dinucleotide,
    conditional_entropy,
)


def main():
    print("Analytical entropy metrics for the A/C-only Drt3b model")
    print("-" * 72)
    print(f"{'Mutant':<12}{'Delta':>8}{'eps':>10}"
          f"{'H1':>8}{'H2':>8}{'H(X2|X1)':>12}")

    for name, delta in DELTA_AC_DEFAULTS.items():
        if delta >= 50:
            eps = 0.0
        else:
            eps = epsilon_from_delta(delta)
        h1 = 1.0                                # A/C-symmetric
        h2 = H2_dinucleotide(eps)
        hc = conditional_entropy(eps)
        print(f"{name:<12}{delta:>8.2f}{eps:>10.4f}"
              f"{h1:>8.3f}{h2:>8.3f}{hc:>12.3f}")

    print()
    print("Random four-base limit (reference):")
    print(f"{'Random':<12}{0.0:>8.2f}{0.25:>10.4f}"
          f"{2.000:>8.3f}{4.000:>8.3f}{2.000:>12.3f}")


if __name__ == "__main__":
    main()
