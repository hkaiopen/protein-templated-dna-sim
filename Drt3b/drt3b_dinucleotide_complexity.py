#!/usr/bin/env python3
"""
Drt3b Mutation Effect: Dinucleotide Frequency & Sequence Complexity
====================================================================
Extends the stochastic Boltzmann model to compute:
1. Dinucleotide frequency matrix (e.g., P(AA), P(AC), ...)
2. Shannon entropy (per symbol) and Lempel-Ziv complexity (normalized)

These metrics quantify the degree of alternation loss and residual order.
"""

import random
import numpy as np
from collections import Counter
from math import log2

# ---------- Core Simulator (same as before) ----------
class Drt3bMutantSimulator:
    def __init__(self, mutation='wildtype', max_length=300, temperature=1.0):
        self.mutation = mutation
        self.max_length = max_length
        self.temperature = temperature
        self.state = 'A'
        self.chain = []
        self.transition = {'A': 'C', 'C': 'A'}

        if mutation == 'wildtype':
            self.allowed = {'A', 'C'}
            self.match_energy = 0.0
            self.mismatch_energy = 100.0
        elif mutation in ('E26A', 'R253A'):
            self.allowed = {'A', 'C'}
            self.match_energy = 0.0
            self.mismatch_energy = 2.5
        elif mutation == 'E26A_R253A':
            self.allowed = {'A', 'C'}
            self.match_energy = 0.0
            self.mismatch_energy = 1.0
        elif mutation == 'random':
            self.allowed = {'A', 'C', 'G', 'T'}
            self.match_energy = 0.0
            self.mismatch_energy = 0.0
        else:
            raise ValueError(f"Unknown mutation: {mutation}")

    def _energy(self, nucleotide):
        if nucleotide not in self.allowed:
            return float('inf')
        if nucleotide == self.state:
            return self.match_energy / self.temperature
        else:
            return self.mismatch_energy / self.temperature

    def step(self):
        candidates = ['A', 'C', 'G', 'T']
        energies = [self._energy(nt) for nt in candidates]
        weights = [np.exp(-e) if e < 1e6 else 0.0 for e in energies]
        if sum(weights) == 0:
            return False
        probs = np.array(weights) / sum(weights)
        chosen = np.random.choice(candidates, p=probs)
        self.chain.append(chosen)
        self.state = self.transition[self.state]
        return True

    def run(self):
        for _ in range(self.max_length):
            if not self.step():
                break
        return ''.join(self.chain)


# ---------- Analysis Functions ----------
def dinucleotide_frequencies(seq):
    """Return a dictionary of normalized dinucleotide frequencies (2x2 for A/C; full 4x4 if G/T present)."""
    counts = Counter()
    for i in range(len(seq)-1):
        counts[seq[i:i+2]] += 1
    total = sum(counts.values())
    if total == 0:
        return {}
    return {k: v/total for k, v in counts.items()}

def shannon_entropy(seq, base_probs=None):
    """Shannon entropy of the sequence (bits per symbol)."""
    if not seq:
        return 0.0
    counts = Counter(seq)
    total = len(seq)
    if base_probs is None:
        probs = [c/total for c in counts.values()]
    else:
        probs = [base_probs.get(b, 0.0) for b in sorted(base_probs.keys())]
    return -sum(p * log2(p) for p in probs if p > 0)

def lempel_ziv_complexity(seq):
    """Lempel-Ziv complexity (number of distinct substrings)."""
    n = len(seq)
    if n == 0:
        return 0
    substrings = set()
    i = 0
    while i < n:
        j = i + 1
        while j <= n and seq[i:j] in substrings:
            j += 1
        if j <= n:
            substrings.add(seq[i:j])
        i = j
    # Normalize by n / log2(n) (classic LZ complexity per symbol)
    if n == 1:
        return 1.0
    return len(substrings) / (n / log2(n))

def average_metrics(sequences):
    """Compute average dinucleotide frequencies, entropy, and LZ complexity over a list of sequences."""
    dinuc_accum = Counter()
    entropies = []
    lz_complexities = []
    for seq in sequences:
        dinuc = dinucleotide_frequencies(seq)
        dinuc_accum.update(dinuc)
        entropies.append(shannon_entropy(seq))
        lz_complexities.append(lempel_ziv_complexity(seq))
    avg_entropy = np.mean(entropies)
    avg_lz = np.mean(lz_complexities)
    # Average dinucleotide frequencies (normalized over all sequences)
    total_dinuc = sum(dinuc_accum.values())
    if total_dinuc > 0:
        avg_dinuc = {k: v/total_dinuc for k, v in dinuc_accum.items()}
    else:
        avg_dinuc = {}
    return avg_dinuc, avg_entropy, avg_lz


# ---------- Main Experiment ----------
def run_complexity_analysis():
    mutations = ['wildtype', 'E26A', 'R253A', 'E26A_R253A', 'random']
    n_trials = 50          # more trials for stable statistics
    seq_length = 300
    temperature = 1.0

    results = {}
    for mut in mutations:
        sequences = []
        for _ in range(n_trials):
            sim = Drt3bMutantSimulator(mutation=mut, max_length=seq_length, temperature=temperature)
            seq = sim.run()
            sequences.append(seq)
        avg_dinuc, avg_entropy, avg_lz = average_metrics(sequences)
        results[mut] = (avg_dinuc, avg_entropy, avg_lz)

    # Print summary tables
    print("=== Drt3b Mutation: Dinucleotide Frequencies (average over {} trials) ===\n".format(n_trials))
    for mut, (dinuc, _, _) in results.items():
        print(f"{mut:12}: ", end='')
        if not dinuc:
            print("no data")
            continue
        # Format only top 4 most frequent dinucleotides or all if small
        items = sorted(dinuc.items(), key=lambda x: -x[1])[:6]
        print(', '.join(f"{k}={v:.3f}" for k, v in items))

    print("\n=== Sequence Complexity (average over {} trials) ===\n".format(n_trials))
    print(f"{'Mutation':<12} | {'Shannon Entropy (bits)':<22} | {'LZ Complexity (norm)':<20}")
    print("-" * 60)
    for mut, (_, entropy, lz) in results.items():
        print(f"{mut:<12} | {entropy:<22.4f} | {lz:<20.4f}")

    # Additional: for random mutation, expected entropy for uniform 4 symbols = 2 bits
    # For A/C only uniform = 1 bit. Our values should approach these.
    print("\n--- Interpretation ---")
    print("Wildtype: entropy ~ 1.0 bit (only A/C, perfect alternation gives some structure, "
          "so entropy may be slightly less than 1 due to constraints?) Actually perfect alternating AC"
          "has entropy 1.0 bit per symbol (equal probability of A and C, no correlation). "
          "LZ complexity < 1 indicates redundancy.")
    print("E26A/R253A: entropy still near 1.0 but LZ may increase due to occasional repeats.")
    print("E26A_R253A: entropy near 1.0, LZ higher, approaching random A/C.")
    print("Random: entropy near 2.0 bits (four symbols), LZ near 1.0 (fully random).")

if __name__ == "__main__":
    np.random.seed(42)
    random.seed(42)
    run_complexity_analysis()