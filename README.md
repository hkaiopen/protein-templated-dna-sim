# Protein-Templated DNA Synthesis — An Information Dynamics Approach

This repository contains a complete computational framework for modeling **protein‑templated DNA synthesis**, using **Drt3b** — the first experimentally discovered reverse transcriptase that synthesizes DNA without any nucleic acid template — as the primary example.

Building on the **Real‑Imaginary Duality Principle** of information dynamics [1,2], we formalize the mechanism of Drt3b as a **two‑state automaton** driven by a **Boltzmann selection rule** (the coupling matrix). The model quantitatively reproduces all key experimental observations:
- Wild‑type: strictly alternating **ACACAC…** chain.
- Single‑point mutations (E26A, R253A): partial loss of alternation (error rate ~0.14).
- Double mutation (E26A_R253A): near‑random A/C sequence (error rate ~0.40).
- Fully relaxed mutant (G/T allowed): unbiased random sequence (error rate ~0.24, Shannon entropy ~2 bits).

All simulations are **self‑contained**, require **no training data**, and run on any standard Python environment.

---

## Table of Contents
- [Scientific Background](#scientific-background)
- [The Information Dynamics Framework](#the-information-dynamics-framework)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Running the Simulations](#running-the-simulations)
  - [1. Deterministic Greedy Simulator](#1-deterministic-greedy-simulator)
  - [2. Noise Robustness Experiment](#2-noise-robustness-experiment)
  - [3. Mutation Effects Experiment](#3-mutation-effects-experiment)
  - [4. Dinucleotide Frequency & Complexity Analysis](#4-dinucleotide-frequency--complexity-analysis)
- [Interpreting the Results](#interpreting-the-results)
- [Citation](#citation)
- [References](#references)

---

## Scientific Background

In 2026, Deng et al. reported in *Science* a groundbreaking discovery: the bacterial reverse transcriptase **Drt3b** uses two amino acid residues (Glu26 and Arg253) in its own active site as a physical **mold** to synthesize strictly alternating **poly(AC)** single‑stranded DNA — in the complete absence of any nucleic acid template [3]. This finding revealed an entirely new paradigm of template‑directed polymerization.

However, the original work provided a **qualitative** description of how the protein mold forces base alternation, leaving several questions unanswered:
- How does the system “decide” which nucleotide to incorporate at each step?
- How strong is the energetic discrimination between matched and mismatched bases?
- How does the system behave under non‑ideal conditions (temperature changes, substrate fluctuations)?
- How do active‑site mutations quantitatively alter the product statistics?

This repository answers these questions by recasting the Drt3b mechanism within the **information dynamics framework**.

---

## The Information Dynamics Framework

The model is built on three core concepts from information dynamics [1,2]:

| Concept | Definition | In Drt3b |
|:--------|:-----------|:----------|
| **Virtual Space** | The absolute rule set defined by the system's geometry/topology. | The two‑state automaton (A‑state → C‑state → A‑state) enforced by the active pocket. |
| **Real Space** | Observational data (often noisy or incomplete). | The cellular pool of free dNTPs (concentrations + chemical potentials). |
| **Coupling Matrix** | The dynamic projection that drives the system toward a unique steady state. | Boltzmann probability selection: $P \propto e^{-E/kT}$. |

**Energy function** (in units of $kT$):
```
E(n, s) = 0       if n == s and n in {A, C}
           Δ      if n != s and n in {A, C}
           ∞      if n in {G, T} and non‑random mutant
```

**Mutant parameters** (Δ values calibrated from experiments):
- Wild‑type: Δ = 100 (strict exclusion)
- Single‑point mutations (E26A, R253A): Δ = 2.5 (partial exclusion)
- Double mutation (E26A_R253A): Δ = 1.0 (weak exclusion)
- Random mutant (allowing G/T): Δ = 0 (no exclusion)

This framework is **parameter‑free** for the wild‑type and uses only experimentally motivated Δ values for mutants.

---

## Repository Structure

```
protein-templated-dna-sim/
└── Drt3b/
    ├── drt3b_deterministic_simulator.py          # Core deterministic greedy model
    ├── drt3b_noise_robustness.py                 # Experiment 1: temperature & concentration effects
    ├── drt3b_mutation_effects.py                 # Experiment 2: E26A, R253A, double & random mutants
    ├── drt3b_dinucleotide_complexity.py          # Experiment 3: dinucleotide frequencies, entropy, LZ complexity
    ├── drt3b_deterministic_simulator_log.txt     # Sample output log
    ├── drt3b_noise_robustness_log.txt            # Sample output log
    ├── drt3b_mutation_effects_log.txt            # Sample output log
    └── drt3b_dinucleotide_complexity_log.txt     # Sample output log
```

All simulation scripts are **self‑contained** and can be run independently.

---

## Requirements

- **Python 3.7+** (tested on 3.9, 3.10, 3.11)
- Standard library only: `random`, `numpy`, `collections`, `difflib`, `math`

> **Note**: `numpy` is used only for numerical operations (random sampling, exponentials, arrays). No machine‑learning frameworks or training data are required.

Install numpy if missing:
```bash
pip install numpy
```

---

## Running the Simulations

Clone the repository and navigate to the `Drt3b` directory:

```bash
git clone https://github.com/hkaiopen/protein-templated-dna-sim.git
cd protein-templated-dna-sim/Drt3b
```

### 1. Deterministic Greedy Simulator

This is the **core deterministic model**. It uses a pure greedy selection rule (always pick the nucleotide with the lowest energy barrier), demonstrating that strict alternation emerges from the virtual space rules alone.

```bash
python drt3b_deterministic_simulator.py
```

**Expected output** (wild‑type):
```
=== Deterministic Greedy Drt3b Simulator ===
Synthesized DNA chain length: 500 bp
Strictly alternating AC? True
A fraction: 0.50, C fraction: 0.50
Preview: ACACACACACACACACACAC...
✅ Simulation perfectly reproduces the expected strict AC alternating pattern.
```

### 2. Noise Robustness Experiment

Tests the system under **non‑ideal conditions**:
- Varying **temperature** (thermal noise)
- Varying **substrate concentration ratios** (bias for mismatched nucleotides)

```bash
python drt3b_noise_robustness.py
```

**Expected output** (extreme noise: T=10, concentration ratio=10):
```
--- Extreme noise (T=10, conc_ratio=10) ---
Error rate = 0.4422
Sequence preview: CCACAACCCCCCACAAACACACCC...
```
Wild‑type remains error‑free (error rate = 0) across all realistic conditions. Errors appear only when the active site is severely compromised **and** thermal fluctuations are extreme.

### 3. Mutation Effects Experiment

Simulates the effects of active‑site mutations on the product sequence:

```bash
python drt3b_mutation_effects.py
```

**Expected output** (stochastic Boltzmann model, 20 trials each):
```
=== Drt3b Mutation Effect Simulation (Stochastic Boltzmann) ===

wildtype     : error rate = 0.0000 ± 0.0000
E26A         : error rate = 0.1440 ± 0.0244
R253A        : error rate = 0.1413 ± 0.0373
E26A_R253A   : error rate = 0.3970 ± 0.0314
random       : error rate = 0.2436 ± 0.0177
```

**Representative sequences**:
```
wildtype    : ACACACACACACACACACACACACAC...
E26A        : ACACACAAACACCCACACACACACAC...
E26A_R253A  : ACACCCACAAACCCACACAAACACAC...
random      : GCATTCGCCACTAGTAATCACCTGCG...
```

### 4. Dinucleotide Frequency & Complexity Analysis

Quantifies the **loss of alternation** at the level of sequence statistics:

```bash
python drt3b_dinucleotide_complexity.py
```

**Expected output**:
```
=== Drt3b Mutation: Dinucleotide Frequencies ===

wildtype    : AC=0.502, CA=0.498
E26A        : AC=0.432, CA=0.429, AA=0.073, CC=0.066
E26A_R253A  : AC=0.306, CA=0.304, CC=0.200, AA=0.189
random      : CC=0.067, CG=0.066, AT=0.065, TA=0.064, TC=0.064, TT=0.063

=== Sequence Complexity ===

Mutation     | Shannon Entropy (bits) | LZ Complexity (norm)
wildtype     | 1.0000                 | 0.9052
E26A         | 0.9992                 | 1.2853
E26A_R253A   | 0.9985                 | 1.7451
random       | 1.9936                 | 2.6294
```

- **Shannon entropy** close to 1 bit (wild‑type) → only A and C, equally probable.
- **LZ complexity** increases from 0.9 (highly ordered) to >2.5 (random) as alternation is lost.

---

## Interpreting the Results

| Mutant | Δ (energy) | Error Rate | Entropy | LZ Complexity | Interpretation |
|:-------|:-----------|:-----------|:--------|:--------------|:---------------|
| Wild‑type | 100 | 0.000 | 1.00 | 0.91 | Strict alternation (virtual space dominates) |
| E26A / R253A | 2.5 | ~0.14 | 1.00 | 1.29 | Partial alternation loss; occasional repeats |
| E26A_R253A | 1.0 | ~0.40 | 1.00 | 1.75 | Near‑random A/C; weak residual alternation |
| Random (G/T allowed) | 0 | ~0.24 | 1.99 | 2.63 | Unbiased random sequence (four bases) |

**Key insight**: The wild‑type protein mold acts as an **absolute filter** — the virtual space rules (strict A/C alternation) dominate over any real‑space noise. Mutations weaken the filter, moving the system into a **soft‑constraint regime** where thermodynamic fluctuations begin to matter. When the energy difference drops to ~$kT$, the system approaches a **random walk**, and the deterministic steady state is lost.

---

## Citation

If you use this code in your research, please cite the following works:

**The information dynamics framework and its biological validations**:
> H. Liu and K. Huang, “Validation of the Real-Imaginary Duality Principle in Core Challenges of Computational Biology: From Sequencing by Hybridization to RNA Inverse Folding,” Zenodo, 2026.  
> DOI: [10.5281/zenodo.20057468](https://doi.org/10.5281/zenodo.20057468)

**Mathematical foundations**:
> K. Huang, “Mathematical Principles of Information Dynamics — The Universe as a Natural Philosophy of Automatic Control,” Zenodo, 2026.  
> DOI: [10.5281/zenodo.20432225](https://doi.org/10.5281/zenodo.20432225)

**The original experimental discovery (Drt3b)**:
> P. Deng, H. Lee, C. Armijo, H. Wang, and A. Gao, “Protein-templated synthesis of dinucleotide repeat DNA by an antiphage reverse transcriptase,” *Science*, 2026.  
> DOI: [10.1126/science.aed1656](https://doi.org/10.1126/science.aed1656)

**This repository** (optional):
> K. Huang, H. Liu, and Z. Huang, “Protein‑Templated DNA Synthesis — An Information Dynamics Approach,” GitHub repository, 2026.  
> URL: [https://github.com/hkaiopen/protein-templated-dna-sim](https://github.com/hkaiopen/protein-templated-dna-sim)

---

## References

1. H. Liu and K. Huang, “Validation of the Real-Imaginary Duality Principle in Core Challenges of Computational Biology: From Sequencing by Hybridization to RNA Inverse Folding,” *Zenodo*, 2026. DOI: 10.5281/zenodo.20057468
2. K. Huang, “Mathematical Principles of Information Dynamics — The Universe as a Natural Philosophy of Automatic Control,” *Zenodo*, 2026. DOI: 10.5281/zenodo.20432225
3. P. Deng, H. Lee, C. Armijo, H. Wang, and A. Gao, “Protein-templated synthesis of dinucleotide repeat DNA by an antiphage reverse transcriptase,” *Science*, 2026. DOI: 10.1126/science.aed1656

---

## Contact & Contributions

- **Author**: Kai Huang (hkaiopen@foxmail.com)
- **GitHub**: [@hkaiopen](https://github.com/hkaiopen)

Issues and pull requests are welcome. If you use this framework for other protein‑templated systems or wish to collaborate, please get in touch.

---

*“All complex patterns — from prime distributions to DNA assembly — emerge from combinations of four basic operations: diffusion, anti‑diffusion, nonlinear compression, and logarithmic potential. Drt3b is a beautiful experimental manifestation of this principle.”*
```
