"""
Error rate as a function of temperature at fixed Delta (kT units).

Note on interpretation: with energies expressed in kT units, varying T
while holding Delta fixed does not change the Boltzmann factor
exp(-E/kT). The temperature scan is therefore provided as a numerical
check that the implementation is well-behaved; it is NOT the same as
Fig. 2 of the manuscript, which plots eps vs. Delta at T = 1.

For the manuscript's Fig. 2, see drt3b_delta_scan.py.
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
