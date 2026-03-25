# fresco_p_adjust.py

import os
import pandas as pd
import sys
from statsmodels.stats.multitest import multipletests

def adjust_fresco_pvalues(input_file, output_dir):
    # Load FRESCo output, skipping the first line if it's metadata (e.g., window size)
    df = pd.read_csv(input_file, sep='\t', skiprows=1)

    if df.shape[1] < 11:
        raise ValueError("Input file must have at least 11 columns including p_value.")

    # Count total and valid p-values
    pval_col = df.columns[10]
    df[pval_col] = pd.to_numeric(df[pval_col], errors='coerce')
    total_rows = len(df)
    valid_pvals = df[pval_col].notnull()
    valid_count = valid_pvals.sum()
    print(f"Total windows in file: {total_rows}")
    print(f"Number of p-values used for Bonferroni correction: {valid_count}")

    # Bonferroni adjustment on valid p-values
    raw_pvals = df.loc[valid_pvals, pval_col].values
    p_adj_bonf = multipletests(raw_pvals, method='bonferroni')[1]
    df['p_adj_bonf'] = p_adj_bonf

    # Apply SAE/SCE classification logic
    syn_rate_alt = df.iloc[:, 1]  # syn_rate_alternate_model is column 2 (index 1)

    def classify_region(p_adj, syn_rate):
        if p_adj < 0.05 and syn_rate < 1:
            return 'SCE'
        elif p_adj < 0.05 and syn_rate > 1:
            return 'SAE'
        else:
            return 'normal'

    df['classification'] = [classify_region(p, s) for p, s in zip(df['p_adj_bonf'], syn_rate_alt)]

    # Prepare output filename
    os.makedirs(output_dir, exist_ok=True)
    basename = os.path.basename(input_file)
    output_file = os.path.join(output_dir, basename)

    # Save output
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Saved adjusted output to: {output_file}")

if __name__ == '__main__':
    # Example usage: python fresco_p_adjust.py A3GALT2_CCDS60080.1_fresco.txt chr1_p_adj
    if len(sys.argv) != 3:
        print("Usage: python fresco_p_adjust.py <input_fresco_file> <output_directory>")
        sys.exit(1)

    input_fresco = sys.argv[1]
    output_dir = sys.argv[2]
    adjust_fresco_pvalues(input_fresco, output_dir)
