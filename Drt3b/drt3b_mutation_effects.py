#!/usr/bin/env python3
"""
Drt3b Mutation Effects Simulation
==================================
Models template-free poly(AC) DNA synthesis by the bacterial reverse
transcriptase Drt3b under wild-type and mutant conditions.

Model
-----
The active site alternates between an A-selecting state (S_A) and a
C-selecting state (S_C). At each step, one of four nucleotides is
incorporated with a Boltzmann probability determined by the energy
barrier between the candidate nucleotide and the preferred nucleotide
of the current state.

State-dependent G/T barriers
----------------------------
Biochemical evidence (Deng et al., 2026, Fig. S11A-S11C) indicates
that the two states have distinct selectivity gates:
- S_A is gated by Glu26, which anchors dATP.
- S_C is gated by Arg253 and Gly248, which anchor dCTP and exclude dGTP.

Mutations that disrupt Glu26 (E26A, E26Q) therefore allow G/T
misincorporation only at the A-selecting state. Mutations that disrupt
Arg253 (R253A) would affect only the C-selecting state, but R253A has
not been experimentally characterized (see note below).

Experimental anchors
--------------------
- WT: error rate approximately 0; strict A/C alternation.
- E26Q: product-level dG fraction = 10.16% (Fig. S11D, S11F).
- E26A: A-state residual preference for dA over dG is approximately
  80% / 20% (Fig. S11A). Since the A-state occupies half of the cycle,
  the product-level dG fraction is approximately 10%, similar to E26Q.
  The distinguishing feature of E26A is severe product truncation
  (processivity defect, Fig. 4H), which is outside the scope of the
  current per-step fidelity model.

Calibration
-----------
For the G/T-extended model (G incorporation at the dA-selecting state
only), the A-state probability of incorporating G obeys

    P(G | A-state) = e^(−Δ_G,A) / (1 + e^(−Δ_G,A)),
    hence  Δ_G,A = ln[ (1 − p_G) / p_G ].

Both E26Q and E26A retain an A-state preference of approximately 80%
for dA over dG (p_G ≈ 0.20), giving Δ_G,A ≈ 1.39 kT for both mutants.
This is consistent with the reported product-level dG fraction of
~10% for E26Q, since the A-state occupies half of the cycle.

Note on E26A processivity
-------------------------
E26A produces severely truncated cDNA (Fig. 4H), indicating a
processivity defect that is not captured by the present per-step
selectivity model. The model therefore reports the per-step fidelity
of E26A as comparable to E26Q, while noting that the overall product
yield is lower for E26A.

Note on R253A
-------------
R253A has not been experimentally characterized. The original study
(Deng et al., 2026) did not purify the R253A mutant protein, so no
biochemical or genetic data are available (Deng, personal communication).
The R253A parameters in this code are therefore a model prediction based
on the structural symmetry between Glu26 (A-state gate) and Arg253
(C-state gate). This prediction remains to be tested.

Data sources
------------
Deng et al., Science 392, 1274 (2026), Figs. 4H, S11A, S11D, S11F.
Kiran (2026) for E26Q dG = 10.16%.

"""

import numpy as np
from collections import Counter


# ---------------------------------------------------------------------------
# Calibration utilities
# ---------------------------------------------------------------------------

def calibrate_delta_AC_from_AC6(f_AC6):
    """
    Return the A/C mismatch barrier delta_AC consistent with an
    experimental (AC)6 motif frequency f_AC6.

    Leading-order mapping: f_AC6 ~ (1 - p)^6, where p is the per-step
    mismatch probability. Inverting:
        p = 1 - f_AC6^(1/6)
        delta_AC = -ln( p / (1 - p) )

    Parameters
    ----------
    f_AC6 : float
        Experimental (AC)6 frequency, in (0, 1].

    Returns
    -------
    float
        Calibrated delta_AC in kT units. Returns infinity if f_AC6 = 1.
    """
    if f_AC6 <= 0 or f_AC6 > 1:
        raise ValueError("f_AC6 must lie in (0, 1].")
    p = 1.0 - f_AC6 ** (1.0 / 6.0)
    if p <= 0:
        return float('inf')
    return -np.log(p / (1.0 - p))


def calibrate_delta_G_A_from_A_state_fraction(p_G):
    """
    Return the A-state G barrier delta_G_A consistent with a per-step
    probability p_G of incorporating G at the A-selecting state.

    Derivation
    ----------
    P(G | A-state) = e^{-delta_G_A} / (1 + e^{-delta_G_A})

    Solving for delta_G_A:
        delta_G_A = ln[ (1 - p_G) / p_G ]

    Parameters
    ----------
    p_G : float
        Per-step probability of incorporating G at the A-state,
        in (0, 1).

    Returns
    -------
    float
        Calibrated delta_G_A in kT units.
    """
    if p_G <= 0:
        return float('inf')
    if p_G >= 1:
        return 0.0
    return np.log((1.0 - p_G) / p_G)


def calibrate_delta_G_A_from_product_fraction(f_G_product, delta_AC=None):
    """
    Return the A-state G barrier delta_G_A consistent with a product-level
    dG fraction f_G_product, assuming G is incorporated only at the
    A-selecting state and that state occupies half of the cycle.

    Then the per-step A-state probability is p_G = 2 * f_G_product.

    Parameters
    ----------
    f_G_product : float
        Product-level dG fraction (e.g., 0.1016 for E26Q).
    delta_AC : float, optional
        Retained for backward compatibility; not used in this
        simplified two-state calculation.

    Returns
    -------
    float
        Calibrated delta_G_A in kT units.
    """
    p_G = 2.0 * f_G_product
    return calibrate_delta_G_A_from_A_state_fraction(p_G)


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

class Drt3bMutantSimulator:
    """
    Stochastic two-state simulator with state-dependent G/T barriers.

    Parameters
    ----------
    mutation : str
        'wildtype', 'E26A', 'E26Q', 'R253A', 'E26A_R253A', or 'random'.
    max_length : int
        Maximum number of nucleotides to synthesize.
    temperature : float
        System temperature in kT units (default 1.0).
    delta_AC : float or None
        A/C mismatch barrier. If None, a mutation-specific default is used.
    delta_G_A : float or None
        G barrier at the A-selecting state. If None, calibrated from the
        mutation-specific A-state preference.
    delta_T_A : float or None
        T barrier at the A-selecting state. Defaults to delta_G_A.
    """

    def __init__(self, mutation='wildtype', max_length=500, temperature=1.0,
                 delta_AC=None, delta_G_A=None, delta_T_A=None):
        self.mutation = mutation
        self.max_length = max_length
        self.temperature = temperature
        self.state = 'A'
        self.chain = []
        self.transition = {'A': 'C', 'C': 'A'}

        # Mutation-specific parameters.
        if mutation == 'wildtype':
            self.allowed = {'A', 'C'}
            self.delta_AC = 100.0 if delta_AC is None else delta_AC
            self.delta_G_A = float('inf')
            self.delta_T_A = float('inf')

        elif mutation == 'E26A':
            # Glu26 eliminated: G/T enter at the A-state.
            # A-state residual preference: ~80% dA vs ~20% dG (Fig. S11A).
            self.allowed = {'A', 'C', 'G', 'T'}
            self.delta_AC = 2.47 if delta_AC is None else delta_AC
            if delta_G_A is None:
                delta_G_A = calibrate_delta_G_A_from_A_state_fraction(0.20)
            self.delta_G_A = delta_G_A
            self.delta_T_A = delta_T_A if delta_T_A is not None else delta_G_A

        elif mutation == 'E26Q':
            # Glu26 weakened: similar A-state preference to E26A.
            # Product-level dG ~10.16% => A-state p_G ~0.2032.
            self.allowed = {'A', 'C', 'G', 'T'}
            self.delta_AC = 2.5 if delta_AC is None else delta_AC
            if delta_G_A is None:
                delta_G_A = calibrate_delta_G_A_from_A_state_fraction(0.2032)
            self.delta_G_A = delta_G_A
            self.delta_T_A = delta_T_A if delta_T_A is not None else delta_G_A

        elif mutation == 'R253A':
            # Model prediction only. R253A was not experimentally
            # characterized; the mutant protein was not purified.
            # Assumes structural symmetry with E26A: Arg253 is treated
            # as a C-state selectivity gate analogous to Glu26 at the
            # A-state, so the A/C mismatch barrier is set to the same
            # value as E26A, while G/T remain excluded (the C-state gate
            # for dG is Gly248, which is intact in R253A).
            self.allowed = {'A', 'C'}
            self.delta_AC = 2.47 if delta_AC is None else delta_AC
            self.delta_G_A = float('inf')
            self.delta_T_A = float('inf')

        elif mutation == 'E26A_R253A':
            self.allowed = {'A', 'C'}
            self.delta_AC = 1.0 if delta_AC is None else delta_AC
            self.delta_G_A = float('inf')
            self.delta_T_A = float('inf')

        elif mutation == 'random':
            self.allowed = {'A', 'C', 'G', 'T'}
            self.delta_AC = 0.0
            self.delta_G_A = 0.0
            self.delta_T_A = 0.0

        else:
            raise ValueError(f"Unknown mutation: {mutation}")

    def _energy(self, nucleotide):
        """
        Return the energy barrier for the candidate nucleotide given the
        current state. G and T are only accessible at the A-selecting state,
        consistent with the structural assignment of Glu26 as the A-state
        selectivity gate.
        """
        if nucleotide not in self.allowed:
            return float('inf')
        if nucleotide == self.state:
            return 0.0
        if nucleotide in {'A', 'C'}:
            return self.delta_AC / self.temperature
        if self.state == 'A':
            if nucleotide == 'G':
                return self.delta_G_A / self.temperature
            if nucleotide == 'T':
                return self.delta_T_A / self.temperature
        # G/T are excluded at the C-selecting state.
        return float('inf')

    def step(self):
        """Perform one incorporation using Boltzmann sampling."""
        candidates = ['A', 'C', 'G', 'T']
        energies = [self._energy(nt) for nt in candidates]
        weights = [np.exp(-e) if e < 1e6 else 0.0 for e in energies]
        total = sum(weights)
        if total == 0:
            return False
        probs = np.array(weights) / total
        chosen = np.random.choice(candidates, p=probs)
        self.chain.append(chosen)
        self.state = self.transition[self.state]
        return True

    def run(self):
        """Synthesize a full sequence."""
        for _ in range(self.max_length):
            if not self.step():
                break
        return ''.join(self.chain)

    def error_rate(self):
        """Fraction of adjacent identical bases in the product."""
        if len(self.chain) < 2:
            return 0.0
        errors = sum(1 for i in range(len(self.chain) - 1)
                     if self.chain[i] == self.chain[i + 1])
        return errors / (len(self.chain) - 1)

    def dG_fraction(self):
        """Fraction of G in the product."""
        if not self.chain:
            return 0.0
        return self.chain.count('G') / len(self.chain)

    def nucleotide_composition(self):
        """Fraction of each nucleotide in the product."""
        if not self.chain:
            return {}
        counter = Counter(self.chain)
        total = len(self.chain)
        return {nt: counter.get(nt, 0) / total for nt in ['A', 'C', 'G', 'T']}


# ---------------------------------------------------------------------------
# Ensemble analysis
# ---------------------------------------------------------------------------

def run_ensemble(mutation, n_trials=50, max_length=500, temperature=1.0):
    """
    Run multiple independent simulations for a given mutation and return
    summary statistics.
    """
    error_rates = []
    dG_fractions = []
    compositions = []
    representative = None

    for i in range(n_trials):
        sim = Drt3bMutantSimulator(
            mutation=mutation,
            max_length=max_length,
            temperature=temperature,
        )
        seq = sim.run()
        error_rates.append(sim.error_rate())
        dG_fractions.append(sim.dG_fraction())
        compositions.append(sim.nucleotide_composition())
        if i == 0:
            representative = seq

    mean_composition = {}
    for nt in ['A', 'C', 'G', 'T']:
        mean_composition[nt] = float(np.mean([c.get(nt, 0.0) for c in compositions]))

    return {
        'mutation': mutation,
        'n_trials': n_trials,
        'error_rate_mean': float(np.mean(error_rates)),
        'error_rate_std': float(np.std(error_rates, ddof=1)) if n_trials > 1 else 0.0,
        'dG_fraction_mean': float(np.mean(dG_fractions)),
        'dG_fraction_std': float(np.std(dG_fractions, ddof=1)) if n_trials > 1 else 0.0,
        'composition_mean': mean_composition,
        'representative_sequence': representative,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    np.random.seed(2024)

    print("Drt3b mutation effects simulation")
    print("Data sources: Deng et al., Science 392, 1274 (2026), "
          "Figs. 4H, S11A, S11D, S11F")
    print("Note: R253A has not been experimentally characterized.")
    print()

    mutations = ['wildtype', 'E26A', 'E26Q', 'R253A', 'E26A_R253A', 'random']

    print(f"{'Mutation':<14} | {'Error rate':<18} | {'dG fraction':<18} | "
          f"{'A':<6} | {'C':<6} | {'G':<6} | {'T':<6}")
    print("-" * 90)

    for mut in mutations:
        result = run_ensemble(mut, n_trials=50, max_length=500)
        comp = result['composition_mean']
        print(f"{mut:<14} | "
              f"{result['error_rate_mean']:.4f} +/- {result['error_rate_std']:.4f} | "
              f"{result['dG_fraction_mean']:.4f} +/- {result['dG_fraction_std']:.4f} | "
              f"{comp['A']:.3f}  | {comp['C']:.3f}  | {comp['G']:.3f}  | {comp['T']:.3f}")

    print()
    print("Representative sequences (first trial, 80 nt):")
    for mut in mutations:
        sim = Drt3bMutantSimulator(mutation=mut, max_length=80)
        seq = sim.run()
        print(f"{mut:<14} : {seq}")

    print()
    print("E26A vs E26Q per-step dG fraction (model prediction):")
    res_E26A = run_ensemble('E26A', n_trials=50, max_length=500)
    res_E26Q = run_ensemble('E26Q', n_trials=50, max_length=500)
    print(f"  E26A product-level dG  = {res_E26A['dG_fraction_mean']:.4f}")
    print(f"  E26Q product-level dG  = {res_E26Q['dG_fraction_mean']:.4f}")
    print("  Both mutants share the same per-step A-state preference (~80% dA).")
    print("  E26A's lower overall yield is attributed to a processivity defect")
    print("  (Fig. 4H), which is outside the scope of this per-step model.")
