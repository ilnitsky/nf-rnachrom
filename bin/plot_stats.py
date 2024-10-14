import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('stats.tsv', sep='\t', index_col="sample               ")

df_transposed = df.T

plt.figure(figsize=(10, 6))

for column in df_transposed.columns:
    initial_value = df_transposed[column].iloc[0]  
    print(initial_value)
    percentages = df_transposed[column] / initial_value * 100  # Calculate percentages
    
    plt.plot(df_transposed.index, df_transposed[column], marker='o', label=column)
    
    # Annotate each point with its percentage
    for x, y, pct in zip(df_transposed.index, df_transposed[column], percentages):
        plt.text(x, y, f'{pct:.0f}%', ha='center', va='bottom')  

plt.title('Read Counts at Each Stage')
plt.ylabel('Read Counts')
plt.xlabel('Stage')
plt.xticks(rotation=45)  
plt.legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left') 
plt.grid(True)  

plt.tight_layout() 
plt.savefig('read_counts_stage.png', dpi=300) 




# import pandas as pd
# import matplotlib.pyplot as plt

# # Manually transcribing the data into a DataFrame for simplicity
# data = {
#     'Step': ['Dedup', 'Trimming', 'OverlapMerged', 'OverlapUNmerged', 'Debridged', 'RestrSites', 'FilteredAlignedRNA', 'FilteredAlignedDNA', 'RawContacts', 'CigarFiltered'],
#     'arabi_wt_SRR13304067': [980692, 972388, 939634, 32754, 633382, 238492, 41145, 143055, 27881, 27880],
#     'arabi_wt_SRR13304068': [981153, 970798, 936013, 34785, 656867, 326694, 87871, 205032, 61190, 61189]
# }

# df = pd.DataFrame(data)

# df_percentage = df.copy()
# for col in df.columns[1:]:
#     df_percentage[col] = (df[col] / df[col].max()) * 100


# plt.figure(figsize=(10, 6))
# for col in df.columns[1:]:
#     plt.fill_between(df['Step'], df[col], label=col, step='mid', alpha=0.7)
# plt.xticks(rotation=45, ha="right")
# plt.ylabel('Absolute Values')
# plt.title('Process Steps - Absolute Values')
# plt.legend()
# plt.tight_layout()
# plt.savefig("stages_absolute_values.png")


# plt.figure(figsize=(10, 6))
# for col in df_percentage.columns[1:]:
#     plt.fill_between(df_percentage['Step'], df_percentage[col], label=col, step='mid', alpha=0.7)
# plt.xticks(rotation=45, ha="right")
# plt.ylabel('Percentage')
# plt.title('Process Steps - Percentage')
# plt.legend()
# plt.tight_layout()
# plt.savefig("stages_percentage.png")


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Plot bridge summaries from provided files")
#     parser.add_argument("se_file", help="The filename for single-end bridge summary")
#     parser.add_argument("pe_file", help="The filename for paired-end bridge summary")
#     parser.add_argument("save_plot", help="")

#     args = parser.parse_args()

#     main(args.se_file, args.pe_file, args.save_plot)