import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import numpy as np

def create_combined_plots(replica_file, merged_file):

    df_replica = pd.read_csv(replica_file, sep='\s+')
    df_merged = pd.read_csv(merged_file, sep='\s+')
   
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))
    
    # Set style
    sns.set_style("whitegrid")
    
    steps = ['Raw', 'SmartSeqFilter', 'Dedup', 'Trimming', 'OverlapMerged', 
             'Debridged', 'RestrSites', 'UniqueRawContacts', 'FilteredUniqueRawContacts']
    
    for idx, row in df_replica.iterrows():
        sample_name = row['sample']  
        values = [row[step] for step in steps]
        ax1.plot(range(len(steps)), values, marker='o', label=sample_name, 
                linewidth=2, markersize=8)
    
    ax1.set_title('Read/Contact Extinction Plot', fontsize=16, pad=20)
    ax1.set_xlabel('Processing Step', fontsize=12)
    ax1.set_ylabel('Number of Reads/Contacts', fontsize=12)
    ax1.set_xticks(range(len(steps)))
    ax1.set_xticklabels(steps, rotation=45, ha='right')
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=8)  
    
    merged_steps = ['MergedReplicas', 'Voted', 'Singletons']
    
    x = np.arange(len(merged_steps))
    width = 0.8 / len(df_merged)
    
    for idx, row in df_merged.iterrows():
        sample_name = row['sample']  
        values = [row[step] for step in merged_steps]
        ax2.bar(x + idx * width, values, width, label=sample_name, alpha=0.8)
    
    ax2.set_title('Merged Statistics Plot', fontsize=16, pad=20)
    ax2.set_xlabel('Processing Step', fontsize=12)
    ax2.set_ylabel('Number of Contacts', fontsize=12)
    ax2.set_xticks(x + (width * (len(df_merged) - 1)) / 2)
    ax2.set_xticklabels(merged_steps, rotation=45, ha='right')
    ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=8)  

    ax2.yaxis.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    
    plt.savefig('combined_plots.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <replica_stats.tsv> <merged_stats.tsv>")
        sys.exit(1)
    
    create_combined_plots(sys.argv[1], sys.argv[2])