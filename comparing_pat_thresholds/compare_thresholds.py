#!/usr/bin/env python3
import pandas as pd
import glob
import os

print("="*80)
print("UXM DECONVOLUTION COMPARISON: np_thresh 0.8 vs 0.67 (default)")
print("="*80)

# Directories
dir_080 = "uxm_deconvolution_results"
dir_067 = "uxm_deconvolution_results_default_thresh"

# Find all samples
samples = []
for d in glob.glob(f"{dir_080}/*/deconvolution"):
    sample = d.split('/')[1]
    samples.append(sample)

print(f"\nFound {len(samples)} samples to compare:")
for s in samples:
    print(f"  - {s}")

# Compare each sample
all_comparisons = []

for sample in samples:
    print(f"\n{'='*80}")
    print(f"SAMPLE: {sample}")
    print(f"{'='*80}")
    
    # Read U25 results (we'll focus on U25 for simplicity)
    file_080 = f"{dir_080}/{sample}/deconvolution/{sample}_deconv_U25.csv"
    file_067 = f"{dir_067}/{sample}/deconvolution/{sample}_deconv_U25.csv"
    
    if not os.path.exists(file_080) or not os.path.exists(file_067):
        print(f"  WARNING: Missing files for {sample}, skipping...")
        continue
    
    df_080 = pd.read_csv(file_080)
    df_067 = pd.read_csv(file_067)
    
    # Merge on cell type
    comparison = pd.merge(df_080, df_067, on='CellType', suffixes=('_np080', '_np067'))
    
    # Calculate differences
    comparison['difference'] = comparison.iloc[:, 2] - comparison.iloc[:, 1]
    comparison['abs_diff'] = comparison['difference'].abs()
    comparison['percent_change'] = (comparison['difference'] / comparison.iloc[:, 1].replace(0, 1)) * 100
    
    # Add sample name
    comparison['sample'] = sample
    
    # Summary statistics
    print(f"\n  Mean absolute difference: {comparison['abs_diff'].mean():.4f} ({comparison['abs_diff'].mean()*100:.2f}%)")
    print(f"  Max absolute difference: {comparison['abs_diff'].max():.4f} ({comparison['abs_diff'].max()*100:.2f}%)")
    print(f"  Cell types with >1% change: {sum(comparison['abs_diff'] > 0.01)}")
    print(f"  Cell types with >2% change: {sum(comparison['abs_diff'] > 0.02)}")
    
    # Top differences
    print(f"\n  Top 10 cell types with biggest changes:")
    top_changes = comparison.nlargest(10, 'abs_diff')[['CellType', 'difference', 'percent_change']].copy()
    top_changes['difference'] = top_changes['difference'].apply(lambda x: f"{x:+.4f}")
    top_changes['percent_change'] = top_changes['percent_change'].apply(lambda x: f"{x:+.2f}%")
    print(top_changes.to_string(index=False))
    
    # Key brain cell types
    brain_cells = ['Neuron', 'Oligodend', 'Blood-Mono+Macro', 'Endothel', 'Blood-NK', 'Blood-T', 'Blood-B']
    print(f"\n  Key cell types for brain tumors:")
    for cell in brain_cells:
        if cell in comparison['CellType'].values:
            row = comparison[comparison['CellType'] == cell].iloc[0]
            print(f"    {cell:20s}: 0.8={row.iloc[1]:.4f}  →  0.67={row.iloc[2]:.4f}  (Δ={row['difference']:+.4f})")
    
    all_comparisons.append(comparison)

# Combine all samples
if all_comparisons:
    all_df = pd.concat(all_comparisons, ignore_index=True)
    
    print(f"\n{'='*80}")
    print("OVERALL SUMMARY ACROSS ALL SAMPLES")
    print(f"{'='*80}")
    
    print(f"\nOverall mean absolute difference: {all_df['abs_diff'].mean():.4f} ({all_df['abs_diff'].mean()*100:.2f}%)")
    print(f"Overall max absolute difference: {all_df['abs_diff'].max():.4f} ({all_df['abs_diff'].max()*100:.2f}%)")
    
    print(f"\nCell types most affected across all samples:")
    avg_by_cell = all_df.groupby('CellType')['abs_diff'].mean().sort_values(ascending=False)
    print(avg_by_cell.head(15).to_string())
    
    # Save detailed comparison
    all_df.to_csv('threshold_comparison_detailed.csv', index=False)
    print(f"\nDetailed comparison saved to: threshold_comparison_detailed.csv")
    
    # Create summary by sample
    summary = all_df.groupby('sample').agg({
        'abs_diff': ['mean', 'max'],
        'CellType': 'count'
    }).round(4)
    summary.columns = ['Mean_Abs_Diff', 'Max_Abs_Diff', 'Num_CellTypes']
    summary.to_csv('threshold_comparison_summary.csv')
    print(f"Summary by sample saved to: threshold_comparison_summary.csv")
    
    print(f"\n{summary}")

print(f"\n{'='*80}")
print("CONCLUSION")
print(f"{'='*80}")
print("\nThe np_thresh parameter affects the number of CpGs included in the analysis.")
print("- np_thresh=0.8: More conservative, only high-confidence CpG calls")
print("- np_thresh=0.67 (default): More inclusive, balanced confidence/coverage")
print("\nSmall differences (<2-3%) suggest both thresholds produce robust results.")
print("Large differences (>5%) would suggest sensitivity to threshold choice.")
