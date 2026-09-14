"""
Central configuration for the Drt3b information-dynamics simulations.

The model has a **deterministic state transition** (S_A -> S_C -> S_A)
and **stochastic nucleotide emission** governed by a Boltzmann rule:
at each step the state-preferred nucleotide is emitted with probability
q = 1/(1+exp(-Delta)), and the alternative A/C nucleotide otherwise.
The state chain is therefore periodic with period 2 and admits a unique
stationary occupation distribution pi_A = pi_C = 1/2 (not a static
fixed point).

All energy values are in kT units (kT = 1).
"""

# ---------------------------------------------------------------
# A/C-only model: energy barriers (kT units)
# ---------------------------------------------------------------
# Wild-type is the hard-exclusion limit (Delta -> infinity). The
# numerical proxy 100 is a safe surrogate; results are insensitive
# to the exact value for any Delta > 10.
DELTA_WT_PROXY = 100.0

# E26A barrier is calibrated from the experimental error rate
# epsilon = 0.144 via q = 1/(1+exp(-Delta)) and Delta = ln(q/(1-q)).
DELTA_E26A = 2.47

# R253A is predicted to be symmetric to E26A.
DELTA_R253A = DELTA_E26A

# E26A_R253A is a model prediction contingent on residual-barrier
# assumption. Sensitivity to this value is reported in the paper.
DELTA_DOUBLE = 1.0

# Random mutant: no discrimination.
DELTA_RANDOM = 0.0

DELTA_AC_DEFAULTS = {
    "WT": DELTA_WT_PROXY,
    "E26A": DELTA_E26A,
    "R253A": DELTA_R253A,
    "E26A_R253A": DELTA_DOUBLE,
    "Random": DELTA_RANDOM,
}

# ---------------------------------------------------------------
# G/T-extended model (E26Q)
# ---------------------------------------------------------------
# At the A-selecting state, only A and G compete under the two-state
# approximation stated in the manuscript (Sec. 3.2). The A-state G
# incorporation probability is 0.20, from which gamma_G = ln(4).
P_G_AT_A_STATE_E26Q = 0.20
GAMMA_G_E26Q = -__import__("math").log(P_G_AT_A_STATE_E26Q
                                    / (1.0 - P_G_AT_A_STATE_E26Q))  # ln 4

# For E26A, the Gly248 steric gate still excludes G and T.
P_G_AT_A_STATE_E26A = 0.0

# ---------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------
DEFAULT_SEQ_LENGTH = 400          # nt; paper uses 300-500
DEFAULT_N_REPLICATES = 50
DEMO_SEED = 42
ENSEMBLE_SEEDS = list(range(1000, 1000 + DEFAULT_N_REPLICATES))

NUCLEOTIDES_AC = ("A", "C")
NUCLEOTIDES_ACGT = ("A", "C", "G", "T")

# ---------------------------------------------------------------
# Analytical helpers
# ---------------------------------------------------------------
import math


def q_from_delta(delta):
    """Per-step correctness probability q = 1/(1+exp(-Delta))."""
    return 1.0 / (1.0 + math.exp(-delta))


def delta_from_q(q):
    """Inverse of q_from_delta: Delta = ln(q/(1-q))."""
    return math.log(q / (1.0 - q))


def epsilon_from_delta(delta):
    """Closed-form error rate epsilon(Delta) = 2*exp(-Delta)/(1+exp(-Delta))^2."""
    e = math.exp(-delta)
    return 2.0 * e / (1.0 + e) ** 2


def delta_from_epsilon(eps):
    """Invert epsilon(Delta) = 0.05 to obtain the fidelity threshold."""
    # Solve 2*x/(1+x)^2 = eps with x = exp(-Delta), take smaller root.
    disc = 1.0 - 2.0 * eps
    if disc < 0:
        raise ValueError("epsilon must be in [0, 0.5]")
    x = (1.0 - math.sqrt(disc)) / (1.0 + math.sqrt(disc)) if disc > 0 else 1.0
    return -math.log(x)


def h2(p):
    """Binary entropy in bits."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log2(p) - (1.0 - p) * math.log2(1.0 - p)


def H2_dinucleotide(epsilon):
    """Dinucleotide entropy for the A/C-only model (merged phases).

    P(AA)=P(CC)=eps/2, P(AC)=P(CA)=(1-eps)/2.
    """
    if epsilon <= 0.0:
        return 1.0
    if epsilon >= 1.0:
        return 1.0
    return 1.0 + h2(epsilon)


def conditional_entropy(epsilon):
    """H(X2|X1) = h2(epsilon) for the A/C-only model."""
    return h2(epsilon)
