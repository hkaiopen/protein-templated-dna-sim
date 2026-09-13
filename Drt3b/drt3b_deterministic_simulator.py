"""
Exact (deterministic) simulator for the A/C-only Drt3b model.

The state alternates deterministically S_A -> S_C -> S_A. At each step,
the emitted nucleotide is selected to be the state-preferred one when
the drawn uniform is below the Boltzmann-corrected threshold, otherwise
the other A/C nucleotide.
"""

import random
from config import (
    DELTA_AC_DEFAULTS,
    NUCLEOTIDES_AC,
    DEFAULT_SEQ_LENGTH,
    DEMO_SEED,
    epsilon_from_delta,
    q_from_delta,
)


def simulate_sequence(delta, length=DEFAULT_SEQ_LENGTH, seed=DEMO_SEED):
    """Simulate an A/C-only sequence under the deterministic two-state model."""
    rng = random.Random(seed)
    q = q_from_delta(delta)
    state_is_A = True
    seq = []
    for _ in range(length):
        preferred = "A" if state_is_A else "C"
        if rng.random() < q:
            seq.append(preferred)
        else:
            seq.append("C" if state_is_A else "A")
        state_is_A = not state_is_A
    return "".join(seq)


def error_rate_adjacent(seq):
    """Fraction of adjacent identical bases."""
    if len(seq) < 2:
        return 0.0
    return sum(seq[i] == seq[i + 1] for i in range(len(seq) - 1)) / (len(seq) - 1)


def main():
    print("Deterministic A/C-only Drt3b simulator")
    print("-" * 60)
    print(f"{'Mutant':<12}{'Delta':>8}{'Sim eps':>12}{'Analytic':>12}")
    for name, delta in DELTA_AC_DEFAULTS.items():
        seq = simulate_sequence(delta, length=DEFAULT_SEQ_LENGTH, seed=DEMO_SEED)
        eps_sim = error_rate_adjacent(seq)
        eps_an = epsilon_from_delta(delta) if delta < 50 else 0.0
        print(f"{name:<12}{delta:>8.2f}{eps_sim:>12.4f}{eps_an:>12.4f}")


if __name__ == "__main__":
    main()
