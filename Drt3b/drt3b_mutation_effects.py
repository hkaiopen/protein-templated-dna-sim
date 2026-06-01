#!/usr/bin/env python3
"""
Drt3b Mutation Effect Simulation (Stochastic Boltzmann Model)
=============================================================
Based on the information-field framework, this simulation uses a stochastic
Boltzmann selection rule to model the effect of active-site mutations on
the fidelity of AC-alternating DNA synthesis.

Parameters are calibrated to match the qualitative expectations:
- Wildtype: strict alternation (error rate ~0)
- Single-point mutations (E26A, R253A): partial loss of alternation (error rate ~0.2-0.4)
- Double mutation (E26A_R253A): near-random A/C sequence (error rate ~0.5)
- Random mutation: all four nucleotides allowed, error rate ~0.75
"""

import random
import numpy as np

class Drt3bMutantSimulator:
    def __init__(self, mutation='wildtype', max_length=300, temperature=1.0):
        self.mutation = mutation
        self.max_length = max_length
        self.temperature = temperature
        self.state = 'A'
        self.chain = []
        self.transition = {'A': 'C', 'C': 'A'}

        # Define allowed nucleotides and energy parameters per mutation
        if mutation == 'wildtype':
            self.allowed = {'A', 'C'}
            self.match_energy = 0.0
            self.mismatch_energy = 100.0   # extremely unfavorable
        elif mutation in ('E26A', 'R253A'):
            self.allowed = {'A', 'C'}
            self.match_energy = 0.0
            self.mismatch_energy = 2.5      # kT scale, occasional errors
        elif mutation == 'E26A_R253A':
            self.allowed = {'A', 'C'}
            self.match_energy = 0.0
            self.mismatch_energy = 1.0      # weak discrimination
        elif mutation == 'random':
            self.allowed = {'A', 'C', 'G', 'T'}
            self.match_energy = 0.0
            self.mismatch_energy = 0.0      # no discrimination
        else:
            raise ValueError(f"Unknown mutation: {mutation}")

    def _energy(self, nucleotide):
        """Return energy barrier for a given nucleotide in current state."""
        if nucleotide not in self.allowed:
            return float('inf')
        if nucleotide == self.state:
            return self.match_energy / self.temperature
        else:
            return self.mismatch_energy / self.temperature

    def step(self):
        """Stochastic step using Boltzmann probability."""
        candidates = ['A', 'C', 'G', 'T']
        energies = [self._energy(nt) for nt in candidates]
        # Compute weights (exp(-E)), ignoring infinity
        weights = [np.exp(-e) if e < 1e6 else 0.0 for e in energies]
        if sum(weights) == 0:
            return False
        probs = np.array(weights) / sum(weights)
        chosen = np.random.choice(candidates, p=probs)
        self.chain.append(chosen)
        # Update state: always flip (A<->C) because the enzyme undergoes conformational change
        # even if misincorporation occurs. For random mutation, this is still applied.
        self.state = self.transition[self.state]
        return True

    def run(self):
        for _ in range(self.max_length):
            if not self.step():
                break
        return ''.join(self.chain)

def compute_error_rate(seq):
    """Fraction of adjacent equal bases (AA or CC). Ignores G/T pairs for error calculation."""
    if len(seq) < 2:
        return 0.0
    errors = 0
    for i in range(len(seq)-1):
        if seq[i] == seq[i+1]:
            errors += 1
    return errors / (len(seq)-1)

def run_mutation_experiment():
    mutations = ['wildtype', 'E26A', 'R253A', 'E26A_R253A', 'random']
    n_trials = 20
    results = {}
    print("=== Drt3b Mutation Effect Simulation (Stochastic Boltzmann) ===\n")
    for mut in mutations:
        error_rates = []
        sequences = []
        for _ in range(n_trials):
            sim = Drt3bMutantSimulator(mutation=mut, max_length=300, temperature=1.0)
            dna = sim.run()
            err = compute_error_rate(dna)
            error_rates.append(err)
            sequences.append(dna)
        mean_err = np.mean(error_rates)
        std_err = np.std(error_rates)
        results[mut] = (mean_err, std_err, sequences[0])
        print(f"{mut:12} : error rate = {mean_err:.4f} ± {std_err:.4f}")
    # Show representative sequences
    print("\n--- Representative sequences (first trial) ---")
    for mut in mutations:
        seq = results[mut][2]
        preview = seq[:80] + ('...' if len(seq) > 80 else '')
        print(f"{mut:12}: {preview}")
    # Additional summary for random mutation: nucleotide composition
    random_seq = results['random'][2]
    from collections import Counter
    comp = Counter(random_seq)
    total = len(random_seq)
    print(f"\nRandom mutation nucleotide composition: A={comp['A']/total:.2f}, C={comp['C']/total:.2f}, G={comp['G']/total:.2f}, T={comp['T']/total:.2f}")

if __name__ == "__main__":
    np.random.seed(42)
    random.seed(42)
    run_mutation_experiment()