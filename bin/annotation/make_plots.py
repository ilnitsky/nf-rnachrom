import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def make_plots(input_file, output_dir):
    contacts = pd.read_csv(input_file, sep="\t")
    contacts_grouped = contacts[['gene_type', 'N_counts']].groupby(['gene_type']).sum().reset_index().sort_values(by="N_counts", ascending=False)
    contacts_grouped['%'] = 100 * contacts_grouped['N_counts'] / contacts_grouped['N_counts'].sum()
    contacts_grouped = pd.concat([contacts_grouped[contacts_grouped['%'] >= 5], pd.DataFrame.from_dict({'gene_type': ['others biotypes'], 'N_counts': [contacts_grouped[contacts_grouped['%'] < 5]['N_counts'].sum()], '%': [0]})])
    contacts_grouped['%'] = 100 * contacts_grouped['N_counts'] / contacts_grouped['N_counts'].sum()

    plt.figure(figsize=(12, 6))
    sns.barplot(data=contacts_grouped, x="gene_type", y="%")
    plt.xticks(rotation=90)
    plt.ylabel('%')
    plt.xlabel('RNA biotypes')
    plt.title('Percentage of contacts that are accounted for by RNA biotypes')
    # Save the plot
    plt.tight_layout()
    plt.savefig(f"{output_dir}/representation_of_RNA_biotypes.png", dpi=360, bbox_inches="tight")
    plt.close()



    contacts_not_zero = contacts[contacts['N_counts'] != 0]
    N_counts_median = int(contacts_not_zero['N_counts'].median())

    sns.ecdfplot(data=contacts_not_zero, x="N_counts")
    plt.plot([1,contacts_not_zero['N_counts'].max()], [0.5,0.5], 'r-')
    plt.plot([N_counts_median,N_counts_median], [0,1], 'r-')
    plt.xscale("log")
    plt.ylabel('Proportion of RNAs with the corresponding\nnumber of contacts or less')
    plt.xlabel('Number of contacts in i-th RNA')
    plt.title('RNAs representation in data. Half of RNAs has less than {} contacts'.format(N_counts_median))
    # Save the plot
    plt.tight_layout()
    plt.savefig(f"{output_dir}/RNAs_representation_in_data.png", dpi=360, bbox_inches="tight")
    plt.close()


    table_df = contacts.groupby(['gene_type']).size().reset_index().merge(contacts_not_zero.groupby(['gene_type']).size().reset_index(), on='gene_type', how='left').fillna(0)
    table_df = table_df.rename({'gene_type': 'RNA biotype', '0_x': 'The amount of RNAs in gene annotation', '0_y': 'The amount of RNAs that have contacts with chromatin'}, axis='columns')
    table_df['The amount of RNAs that have contacts with chromatin'] = table_df['The amount of RNAs that have contacts with chromatin'].apply(lambda x: int(x))
    table_df = table_df.sort_values(by="The amount of RNAs that have contacts with chromatin", ascending=False)

    fig, axs = plt.subplots()
    axs.axis('tight')
    axs.axis('off')
    the_table = axs.table(cellText=table_df[['RNA biotype', 'The amount of RNAs that have contacts with chromatin', 'The amount of RNAs in gene annotation']].values, colLabels=['RNA biotype', 'The amount of RNAs that have contacts with chromatin', 'The amount of RNAs in gene annotation'], loc='center', cellLoc='left')
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(7)
    the_table.auto_set_column_width(col=list(range(len(table_df.columns))))
    # Save the plot
    plt.savefig(f"{output_dir}/RNAs_representation_in_data_2.png", dpi=360, bbox_inches="tight")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Plot RNA biotypes representation')
    parser.add_argument('-i', '--input', required=True, help='Input contacts file')
    parser.add_argument('-o', '--output-dir', required=True, help='Output directory')
    
    args = parser.parse_args()
    make_plots(args.input, args.output_dir)

if __name__ == "__main__":
    main() 