# MicrobiomeDrug

**Pitch**: Systematic analysis of drug-microbiome interactions — predict how gut microbial enzymes metabolize drugs, using Tanimoto similarity to known metabolizing enzymes and gene family abundance profiles.

**Why it matters**: The gut microbiome modulates drug efficacy and toxicity. No existing claw4s tool predicts drug-metabolizing enzyme profiles from metagenomic data, despite this being critical for precision medicine.

**Status**: 60% skeleton — subagent completes implementation.

**Dependencies**: Python 3.9+, NumPy, SciPy, pandas, scikit-learn, RDKit, matplotlib
**Runtime**: ~5-10 min per sample (CPU-only)
