#!/usr/bin/env python3
"""
Deterministic Greedy Drt3b Simulator
=====================================

Implements the information-field algorithm for Drt3b protein-templated
DNA synthesis, in the large-barrier limit (Δ → ∞).

- Constraint space: 2-state automaton (S_A / S_C) with strict
  alternation S_A → S_C → S_A. The state→nucleotide map is
  m(S_A) = A, m(S_C) = C (Reviewer 1, P2).
- Real space: cellular dNTP pool (all four nucleotides present, with
  state-dependent matching energies set by the active-site geometry).
- Coupling: greedy selection of the nucleotide with the lowest energy
  barrier. No randomness, no free parameters.

This model reproduces the experimental observation that Drt3b
synthesizes strictly alternating ACACAC... DNA chains
(Science 2026, Deng et al.).

Note on Δ = 100 kT (Reviewer 1, P4)
-----------------------------------
Δ = 100 kT is a numerical proxy for the hard-exclusion limit
(Δ → ∞), NOT a literal single-bond energy. The wild-type conclusion
(ε ≈ 0) is unchanged for any Δ ≳ 10: the closed form
ε(Δ) = 2e^(−Δ)/(1+e^(−Δ))² falls below 10⁻³ for Δ ≳ 10.
"""

DELTA_WT_PROXY = 100.0   # numerical proxy for the hard-exclusion limit


class Drt3bDeterministicSimulator:
    def __init__(self, max_length: int = 500):
        self.state = 'A'
        self.chain = []
        self.max_length = max_length

        # Constraint space: absolute transition rules
        self.transition = {'A': 'C', 'C': 'A'}
        # State→nucleotide map (Reviewer 1, P2)
        self.state_to_nt = {'A': 'A', 'C': 'C'}
        # Valid nucleotides (A/C-only layer; G, T have E = ∞)
        self.valid_nucleotides = {'A', 'C'}

    def _energy_barrier(self, nucleotide: str) -> float:
        """
        0.0 for the state-preferred nucleotide m(s);
        DELTA_WT_PROXY for the other valid A/C base;
        +∞ for G, T.
        """
        if nucleotide not in self.valid_nucleotides:
            return float('inf')
        if nucleotide == self.state_to_nt[self.state]:
            return 0.0
        return DELTA_WT_PROXY

    def step(self) -> bool:
        energies = {nt: self._energy_barrier(nt) for nt in ['A', 'C', 'G', 'T']}
        best_nt = min(energies, key=energies.get)
        if best_nt not in self.valid_nucleotides or energies[best_nt] >= float('inf'):
            return False
        self.chain.append(best_nt)
        self.state = self.transition[self.state]
        return True

    def run(self) -> str:
        for _ in range(self.max_length):
            if not self.step():
                break
        return ''.join(self.chain)

    def statistics(self) -> dict:
        seq = ''.join(self.chain)
        length = len(seq)
        if length == 0:
            return {"length": 0, "is_alternating": True, "A_fraction": 0.0}
        a_count = seq.count('A')
        c_count = seq.count('C')
        is_alternating = all(seq[i] != seq[i + 1] for i in range(length - 1))
        return {
            "length": length,
            "is_strictly_alternating": is_alternating,
            "A_fraction": a_count / length,
            "C_fraction": c_count / length,
            "preview": seq[:80] + ("..." if length > 80 else ""),
        }


if __name__ == "__main__":
    sim = Drt3bDeterministicSimulator(max_length=500)
    dna = sim.run()
    stats = sim.statistics()
    print("=== Deterministic Greedy Drt3b Simulator ===")
    print(f"Synthesized DNA chain length: {stats['length']} bp")
    print(f"Strictly alternating AC? {stats['is_strictly_alternating']}")
    print(f"A fraction: {stats['A_fraction']:.2f}, C fraction: {stats['C_fraction']:.2f}")
    print(f"Preview: {stats['preview']}")
    if stats['is_strictly_alternating'] and stats['length'] == 500:
        print("\nOK: strict AC alternation reproduced.")
    else:
        print("\nWARNING: output deviates from the ideal alternation.")
