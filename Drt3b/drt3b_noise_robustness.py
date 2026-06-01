#!/usr/bin/env python3
"""
Noise Robustness Experiment for Drt3b
======================================
Simulates non‑ideal conditions:
- Substrate concentration fluctuations (dATP/dCTP ratio)
- Temperature changes (modulates the energy barrier discrimination)
- Random mis‑incorporation probability

Measures the error rate (fraction of non‑alternating positions) as a function of noise level.
"""

import random
import numpy as np
from collections import Counter
from drt3b_deterministic_simulator import Drt3bDeterministicSimulator


class NoisyDrt3bSimulator(Drt3bDeterministicSimulator):
    """
    Extends the deterministic simulator with tunable noise parameters.
    """

    def __init__(self, max_length=500, temperature=0.5, concentration_ratio=1.0):
        """
        Args:
            max_length: maximum chain length
            temperature: controls the softness of the Boltzmann selection.
                         Higher T -> more random, lower T -> more deterministic.
            concentration_ratio: relative concentration of the non‑matching
                                 but still valid nucleotide (e.g., when state=A,
                                 concentration_ratio multiplies the effective
                                 "fuel" of dCTP relative to dATP). Here we
                                 implement it as a bias in the energy.
        """
        super().__init__(max_length)
        self.temperature = temperature
        self.concentration_ratio = concentration_ratio

    def _energy_barrier(self, nucleotide: str) -> float:
        """
        Energy barrier with temperature scaling and concentration bias.
        For invalid nucleotides (G/T) still infinity.
        For valid nucleotides:
            - perfect match: base energy 0.0
            - other valid: base energy E_mismatch = 100.0 / concentration_ratio
                           (higher ratio makes mismatch easier)
        Then the effective energy is E / temperature (Boltzmann factor).
        """
        if nucleotide not in self.valid_nucleotides:
            return float('inf')
        if nucleotide == self.state:
            return 0.0 / self.temperature   # still zero, but kept for consistency
        # Mismatch energy: higher concentration_ratio -> lower barrier
        base_mismatch = 100.0 / max(self.concentration_ratio, 0.001)
        return base_mismatch / self.temperature

    def step_stochastic(self) -> bool:
        """
        Stochastic step using Boltzmann probabilities instead of deterministic greedy.
        Probability of selecting a nucleotide is proportional to exp(-E / T)
        (with T already baked into energies).
        """
        nucleotides = ['A', 'C', 'G', 'T']
        energies = [self._energy_barrier(nt) for nt in nucleotides]
        # Softmax (Boltzmann)
        valid = [np.exp(-e) if e < 1e6 else 0.0 for e in energies]
        if sum(valid) == 0:
            return False
        probs = np.array(valid) / sum(valid)
        chosen = np.random.choice(nucleotides, p=probs)
        if chosen not in self.valid_nucleotides:
            return False
        self.chain.append(chosen)
        self.state = self.transition[self.state]
        return True

    def run_stochastic(self):
        """Run synthesis using stochastic steps."""
        for _ in range(self.max_length):
            if not self.step_stochastic():
                break
        return ''.join(self.chain)


def compute_error_rate(sequence: str) -> float:
    """Calculate fraction of positions where a nucleotide repeats (AA or CC)."""
    if len(sequence) < 2:
        return 0.0
    errors = 0
    for i in range(len(sequence)-1):
        if sequence[i] == sequence[i+1]:
            errors += 1
    return errors / (len(sequence)-1)


def run_noise_experiment():
    """Test robustness under various noise levels."""
    print("=== Drt3b Noise Robustness Experiment ===")
    print("Measuring error rate (non‑alternating adjacent pairs) as a function of noise.\n")

    # Parameters
    max_len = 200
    n_trials = 20

    # Vary temperature (higher = more randomness)
    temps = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    # Vary concentration ratio ( >1 makes mismatched nucleotide more available)
    conc_ratios = [0.5, 1.0, 2.0, 5.0, 10.0]

    print("--- Effect of Temperature (T) ---")
    for T in temps:
        error_rates = []
        for _ in range(n_trials):
            sim = NoisyDrt3bSimulator(max_length=max_len, temperature=T, concentration_ratio=1.0)
            dna = sim.run_stochastic()
            err = compute_error_rate(dna)
            error_rates.append(err)
        mean_err = np.mean(error_rates)
        std_err = np.std(error_rates)
        print(f"T = {T:4.1f} : mean error rate = {mean_err:.4f} ± {std_err:.4f}")

    print("\n--- Effect of Concentration Ratio (bias for mismatched nucleotide) ---")
    for cr in conc_ratios:
        error_rates = []
        for _ in range(n_trials):
            sim = NoisyDrt3bSimulator(max_length=max_len, temperature=0.5, concentration_ratio=cr)
            dna = sim.run_stochastic()
            err = compute_error_rate(dna)
            error_rates.append(err)
        mean_err = np.mean(error_rates)
        std_err = np.std(error_rates)
        print(f"Conc ratio = {cr:4.1f} : mean error rate = {mean_err:.4f} ± {std_err:.4f}")

    # Additional test: combined extreme noise
    print("\n--- Extreme noise (T=10, conc_ratio=10) ---")
    sim = NoisyDrt3bSimulator(max_length=max_len, temperature=10.0, concentration_ratio=10.0)
    dna = sim.run_stochastic()
    err = compute_error_rate(dna)
    print(f"Error rate = {err:.4f}")
    print(f"Sequence preview: {dna[:80]}...")


if __name__ == "__main__":
    run_noise_experiment()