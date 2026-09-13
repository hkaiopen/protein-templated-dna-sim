"""
Ensemble stochastic simulator for the A/C-only Drt3b model.

Produces mean +/- standard deviation of the adjacent-identical error
rate over 50 independent replicates, matching Table 1 of the manuscript.
"""

import random
from statistics import mean, stdev
from config import (
    DELTA_AC_DEFAULTS,
    DEFAULT_SEQ_LENGTH,
    ENSEMBLE_SEEDS,
    epsilon_from_delta,
    q_from_delta,
)


def simulate_once(delta, seed, length=DEFAULT_SEQ_LENGTH):
    rng = random.Random(seed)
    q = q_from_delta(delta)
    state_is_A = True
    seq = []
    for _ in range(length):
        preferred = "A" if state_is_A else "C"
        seq.append(preferred if rng.random() < q else ("C" if state_is_A else "A"))
        state_is_A = not state_is_A
    return seq


def error_rate(seq):
    return sum(seq[i] == seq[i + 1] for i in range(len(seq) - 1)) / (len(seq) - 1)


def ensemble(delta, seeds=ENSEMBLE_SEEDS, length=DEFAULT_SEQ_LENGTH):
    rates = [error_rate(simulate_once(delta, s, length)) for s in seeds]
    return mean(rates), (stdev(rates) if len(rates) > 1 else 0.0)


def main():
    print("Ensemble A/C-only Drt3b simulator (50 replicates, L=400 nt)")
    print("-" * 72)
    print(f"{'Mutant':<12}{'Delta':>8}{'Sim mean':>12}{'Sim std':>10}{'Analytic':>12}")
    for name, delta in DELTA_AC_DEFAULTS.items():
        m, sd = ensemble(delta)
        eps_an = epsilon_from_delta(delta) if delta < 50 else 0.0
        print(f"{name:<12}{delta:>8.2f}{m:>12.4f}{sd:>10.4f}{eps_an:>12.4f}")


if __name__ == "__main__":
    main()