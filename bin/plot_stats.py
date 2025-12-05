#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

sns.set_style("whitegrid")

def plot_df(df, ax, title):
    """Plot every numeric column as a line (for per-replicate) or bar (for merged)"""
    # Drop the sample column, keep everything numeric
    numeric_cols = df.select_dtypes(include='number').columns
    data = df.set_index('sample')[numeric_cols]

    if len(data) <= 8:               # few samples → lines (replicates)
        data.T.plot(ax=ax, marker='o', linewidth=2, markersize=6)
    else:                            # many samples → bars (usually merged stats)
        data.plot(kind='bar', ax=ax, width=0.8, alpha=0.85)

    ax.set_title(title, fontsize=14, pad=15)
    ax.set_ylabel("Count")
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    ax.grid(True, axis='y', linestyle='--', alpha=0.7)

# ------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: plot_stats.py <replica_stats.tsv> <merged_stats.tsv>")
        sys.exit(1)

    replica_file = sys.argv[1]
    merged_file   = sys.argv[2]

    # Read files (ignore blank lines / comments)
    df_rep = pd.read_csv(replica_file, sep=r'\s+', comment='#', engine='python')
    df_mer = pd.read_csv(merged_file,   sep=r'\s+', comment='#', engine='python', on_bad_lines='skip')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    plot_df(df_rep, ax1, "Per-replicate statistics")
    ax1.tick_params(axis='x', rotation=50)

    if not df_mer.empty and len(df_mer.columns) > 1:
        plot_df(df_mer, ax2, "After merging replicates")
    else:
        ax2.text(0.5, 0.5, "No merged\nstatistics", ha='center', va='center',
                 transform=ax2.transAxes, fontsize=16, color='gray')
        ax2.set_xticks([])
        ax2.set_yticks([])

    ax2.tick_params(axis='x', rotation=50)
    plt.tight_layout()
    plt.savefig("combined_stats.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("→ combined_stats.png created")