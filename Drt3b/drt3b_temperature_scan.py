"""
Temperature scan at fixed Delta (kT units).

Note: in kT units the Boltzmann factor exp(-E/kT) is invariant under
rescaling T while holding Delta fixed. This script is therefore a
numerical sanity check, not the epsilon-Delta curve of Figure 2.
For Figure 2, use drt3b_delta_scan.py.
"""

from config import epsilon_from_delta, DELTA_E26A, DELTA_DOUBLE

T_GRID = [0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]


def main():
    print("Temperature scan at fixed Delta (kT units)")
    print("-" * 60)
    print(f"{'T':>6}{'eps(WT proxy)':>16}{'eps(E26A)':>14}{'eps(double)':>14}")
    for t in T_GRID:
        # Energies fixed in kT units; the effective Delta/kT ratio is Delta/T.
        eps_wt = epsilon_from_delta(100.0 / t)
        eps_e26a = epsilon_from_delta(DELTA_E26A / t)
        eps_dbl = epsilon_from_delta(DELTA_DOUBLE / t)
        print(f"{t:>6.1f}{eps_wt:>16.4f}{eps_e26a:>14.4f}{eps_dbl:>14.4f}")


if __name__ == "__main__":
    main()
