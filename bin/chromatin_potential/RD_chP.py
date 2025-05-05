#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats.multitest import multipletests
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Calculate chromatin potential and generate plots')
    parser.add_argument('--input_path', type=str, required=True,
                      help='Input directory path')
    parser.add_argument('--output_path', type=str, required=True,
                      help='Output directory path')
    parser.add_argument('--input_path_RNAseq', type=str, required=True,
                      help='Input RNA-seq directory path')
    parser.add_argument('--counts_RNAseq', type=str, required=True,
                      help='RNA-seq counts filename')
    parser.add_argument('--counts_contacts', type=str, required=True,
                      help='Contacts counts filename')
    parser.add_argument('--N_contacts_min', type=int, default=100,
                      help='Minimum number of contacts (default: 100)')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                      help='FDR threshold (default: 0.05)')
    parser.add_argument('--type', type=str, required=True,
                      help='Type of contacts (UU_all or UU_dist or UM_all or UM_dist)')
    return parser.parse_args()

def calculate_stats(row, Ne, Nc):
    #P_rna_RNAseq = row['N_counts_RNAseq'] / Ne
    #P_rna_contacts = row['N_contacts'] / Nc
    #P = (row['N_counts_RNAseq'] + row['N_contacts']) / (Ne + Nc)
    #Z = (P_rna_contacts - P_rna_RNAseq) / np.sqrt(P * (1 - P) * (1/Nc + 1/Ne))
    
    # Proportion test
    Z, pval = proportions_ztest([row['N_contacts'], row['N_counts_RNAseq']], [Nc, Ne])
    
    return pd.Series({'chP': Z,'pval': pval})

def main():
    args = parse_args()
    
    # Read input files
    counts_RNAseq = pd.read_csv(os.path.join(args.input_path_RNAseq, args.counts_RNAseq), sep='\t')
    counts_RNAseq = counts_RNAseq[['gene_name', 'gene_type', 'N_counts']].rename(columns={'N_counts': 'N_counts_RNAseq'})

    counts_contacts = pd.read_csv(os.path.join(args.input_path, args.counts_contacts), sep='\t')
    counts_contacts = counts_contacts[['gene_name', 'gene_type', 'N_counts']].rename(columns={'N_counts': 'N_contacts'})

    # Merge and filter data
    counts_merge = counts_RNAseq.merge(counts_contacts, how='left')
    counts_merge = counts_merge[counts_merge['N_contacts'] > args.N_contacts_min]
    counts_merge[['N_counts_RNAseq', 'N_contacts']] += 1

    # Calculate statistics
    Nc = counts_merge['N_contacts'].sum()
    Ne = counts_merge['N_counts_RNAseq'].sum()
    
    chP = counts_merge.join(counts_merge.apply(lambda x: calculate_stats(x, Ne, Nc), axis=1))
    chP['fdr_bh'] = multipletests(chP['pval'], method='fdr_bh')[1]

    # Save results
    chP.to_csv(os.path.join(args.output_path, 'chP_{type}.tab'.format(type=args.type)), sep='\t', index=False)

    # Create plots
    unique_gene_types = chP[chP['fdr_bh'] < args.fdr_threshold].sort_values(by="N_contacts")['gene_type'].unique()
    num_plots = len(unique_gene_types)

    # Calculate number of rows based on number of plots
    if num_plots <= 3:
        nrows = 1
    elif num_plots <= 6:
        nrows = 2
    elif num_plots <= 9:
        nrows = 3
    elif num_plots <= 12:
        nrows = 4
    elif num_plots <= 15:
        nrows = 5
    else:
        nrows = 6

    # Calculate number of columns needed
    ncols = (num_plots + nrows - 1) // nrows

    if chP[chP['fdr_bh'] < args.fdr_threshold].shape[0] != 0:

        fig, axes = plt.subplots(nrows=nrows, figsize=(15, nrows*3), ncols=ncols, sharey=True)
        axes = axes.flatten()

        # Set X-axis limits based on filtered data
        xlim = (chP[chP['fdr_bh'] < args.fdr_threshold].chP.min(), chP[chP['fdr_bh'] < args.fdr_threshold].chP.max())

        # For each unique value of gene_type we construct a scatterplot
        for ax, gene in zip(axes, unique_gene_types):
            subset = chP[(chP['gene_type'] == gene) & (chP['fdr_bh'] < args.fdr_threshold)]
            ax.scatter(subset['chP'], subset['N_contacts'], label=gene, alpha=0.9, edgecolor='white')
            ax.set_xlabel('Chromatin potential')
            ax.set_ylabel('N contacts')
            ax.set_xlim(xlim)
            ax.set_yscale('log')
            ax.axvline(x=0, color='red', linestyle='--')
            ax.grid()
            ax.legend()

        # Hide unused subplots
        for i in range(len(unique_gene_types), len(axes)):
            fig.delaxes(axes[i])

        plt.suptitle('Distribution of chromatin potential across all RNA biotypes (FDR < {fdr_threshold}, N contacts > {N_contacts_min}, {type})'.format(fdr_threshold=args.fdr_threshold, N_contacts_min=args.N_contacts_min, type=args.type), fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        plt.savefig(os.path.join(args.output_path, 'chP_{type}.png'.format(type=args.type)), dpi=360, bbox_inches="tight")
        plt.close()

if __name__ == '__main__':
    main()