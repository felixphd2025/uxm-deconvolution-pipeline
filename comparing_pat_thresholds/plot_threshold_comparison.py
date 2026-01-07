#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read detailed comparison
df = pd.read_csv('threshold_comparison_detailed.csv')

# Get unique samples
samples = df['sample'].unique()

# For each sample, plot top 15 cell types side by side
for sample in samples:
    sample_df = df[df['sample'] == sample].copy()
    
    # Get top 15 by average of both thresholds
    sample_df['avg'] = (sample_df.iloc[:, 1] + sample_df.iloc[:, 2]) / 2
    top15 = sample_df.nlargest(15, 'avg')
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    x = np.arange(len(top15))
    width = 0.35
    
    ax.bar(x - width/2, top15.iloc[:, 1], width, label='np_thresh=0.8', alpha=0.8)
    ax.bar(x + width/2, top15.iloc[:, 2], width, label='np_thresh=0.67 (default)', alpha=0.8)
    
    ax.set_xlabel('Cell Type', fontsize=12)
    ax.set_ylabel('Fraction', fontsize=12)
    ax.set_title(f'Threshold Comparison: {sample}', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(top15['CellType'], rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'threshold_comparison_{sample}.pdf')
    print(f"Saved: threshold_comparison_{sample}.pdf")
    plt.close()

print("\nAll comparison plots created!")
