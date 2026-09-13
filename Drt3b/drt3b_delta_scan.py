"""
Scan of the analytical error rate epsilon(Delta) at T = 1 (G/T excluded).

Reproduces Fig. 2 and Table 2 of the manuscript:
    Delta     = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0]
    eps(Delta)= [0.499, 0.470, 0.393, 0.298, 0.210, 0.140, 0.090, 0.013, 1e-4]
"""

from config import epsilon_from_delta, delta_from_epsilon

DELTA_GRID = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0]


def main():
    print("Analytical epsilon(Delta) at T = 1 (G/T excluded)")
    print("-" * 60)
    print(f"{'Delta':>8}{'eps(Delta)':>14}")
    for d in DELTA_GRID:
        print(f"{d:>8.2f}{epsilon_from_delta(d):>14.4f}")

    threshold = delta_from_epsilon(0.05)
    print()
    print(f"Fidelity threshold: eps = 0.05 at Delta = {threshold:.3f} kT")
    print("Prediction D: 'faithful' regime begins at Delta >~ 3.7 kT.")


if __name__ == "__main__":
    main()