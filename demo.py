#!/usr/bin/env python3
"""
Demo for MicrobiomeDrug.

Generates synthetic Pfam abundance matrix to demonstrate the workflow.
Expected output: enzyme_scores.csv + drug_interactions.csv + heatmap.html
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path

def generate_demo_data():
    """Generate synthetic metagenomic Pfam abundance data."""
    
    np.random.seed(42)
    
    n_samples = 20
    n_pfams = 100
    
    # Pfam families (including some drug-metabolizing ones)
    pfam_ids = [f"PF{(i+1):05d}" for i in range(n_pfams)]
    
    # Some relevant Pfams
    relevant_pfams = {
        "PF00067": "CYP450",    # Cytochrome P450
        "PF02798": "GST",       # GST N-terminal
        "PF00043": "GST",        # GST C-terminal  
        "PF00685": "SULT",       # Sulfotransferase
        "PF00201": "UGT",        # UDPGT
        "PF00881": "bacterial_reductases",  # Nitroreductase
    }
    
    # Inject relevant Pfams into the matrix
    for pfam_id in relevant_pfams:
        if pfam_id not in pfam_ids:
            pfam_ids.append(pfam_id)
    
    pfam_ids = sorted(pfam_ids)
    
    # Generate abundance profiles (samples have varying drug metabolism potential)
    samples = [f"SAMPLE_{i:03d}" for i in range(n_samples)]
    
    data = np.zeros((n_samples, len(pfam_ids)))
    
    for i, sample in enumerate(samples):
        # Random baseline
        data[i] = np.random.exponential(0.1, len(pfam_ids))
        
        # Some samples have elevated drug-metabolizing enzymes
        if i < 5:  # samples 0-4: high CYP450
            pfam_idx = pfam_ids.index("PF00067")
            data[i, pfam_idx] = np.random.uniform(0.5, 1.0)
        
        if i >= 5 and i < 10:  # samples 5-9: high bacterial reductases
            pfam_idx = pfam_ids.index("PF00881")
            data[i, pfam_idx] = np.random.uniform(0.5, 1.0)
        
        if i >= 10 and i < 15:  # samples 10-14: high GST
            for pf in ["PF02798", "PF00043"]:
                if pf in pfam_ids:
                    pfam_idx = pfam_ids.index(pf)
                    data[i, pfam_idx] = np.random.uniform(0.3, 0.7)
    
    pfam_df = pd.DataFrame(data, index=samples, columns=pfam_ids)
    pfam_df.to_csv("/tmp/MicrobiomeDrug_demo_pfam.csv")
    
    print("[Demo] Generated:")
    print(f"  - /tmp/MicrobiomeDrug_demo_pfam.csv ({n_samples} samples x {len(pfam_ids)} Pfams)")
    print(f"  - Relevant Pfams: {list(relevant_pfams.keys())}")
    
    return "/tmp/MicrobiomeDrug_demo_pfam.csv"


if __name__ == "__main__":
    pfam_matrix = generate_demo_data()
    
    print("\n[Demo] Expected output structure:")
    print("  microbiome_drug_results/")
    print("  ├── enzyme_scores.csv           # samples x enzyme families")
    print("  ├── drug_interactions.csv       # samples x drugs")
    print("  ├── drug_similarity.csv         # drug x drug Tanimoto similarity")
    print("  └── drug_enzyme_heatmap.html    # interactive heatmap")
    
    print("\n[Demo] Run with:")
    print("  python SKELETON.py \\")
    print("    --pfam-matrix /tmp/MicrobiomeDrug_demo_pfam.csv \\")
    print("    --drugs acetaminophen warfarin metronidazole \\")
    print("    --output-dir microbiome_drug_results")
