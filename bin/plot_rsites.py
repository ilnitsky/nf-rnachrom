import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys

# Set matplotlib to use the 'Agg' backend for non-interactive plots
plt.switch_backend('Agg')

def plot_oligos(ax, data, title, max_count):
    bars = ax.bar(data["Sequence"].astype(str), data["Count"].tolist())
    ax.set_title(title)
    ax.set_xticks(range(len(data["Sequence"])))
    ax.set_xticklabels(data["Sequence"], rotation=45, ha="right")
    ax.set_ylabel('Count')
    ax.set_ylim(0, max_count)

def main(id, prefix):
    file_path = f"{prefix}"
    df = pd.read_csv(file_path, sep="\t", header=None, names=["Sequence", "Type", "Count"], dtype={"Sequence":str, "Type":str, "Count":int})

    dna_5_dinucleotides = df[df["Type"].str.contains("DNA_5'end_last_2")]
    dna_3_dinucleotides = df[df["Type"].str.contains("DNA_3'end_last_2")]
    dna_5_trinucleotides = df[df["Type"].str.contains("DNA_5'end_last_3")]
    dna_3_trinucleotides = df[df["Type"].str.contains("DNA_3'end_last_3")]
    rna_5_dinucleotides = df[df["Type"].str.contains("RNA_5'end_last_2")]
    rna_3_dinucleotides = df[df["Type"].str.contains("RNA_3'end_last_2")]
    rna_5_trinucleotides = df[df["Type"].str.contains("RNA_5'end_last_3")]
    rna_3_trinucleotides = df[df["Type"].str.contains("RNA_3'end_last_3")]

    dna_5_dinucleotides = dna_5_dinucleotides.sort_values(by="Count", ascending=False)
    dna_3_dinucleotides = dna_3_dinucleotides.sort_values(by="Count", ascending=False)
    dna_5_trinucleotides = dna_5_trinucleotides.sort_values(by="Count", ascending=False).head(35)
    dna_3_trinucleotides = dna_3_trinucleotides.sort_values(by="Count", ascending=False).head(35)
    rna_5_dinucleotides = rna_5_dinucleotides.sort_values(by="Count", ascending=False)
    rna_3_dinucleotides = rna_3_dinucleotides.sort_values(by="Count", ascending=False)
    rna_5_trinucleotides = rna_5_trinucleotides.sort_values(by="Count", ascending=False).head(35)
    rna_3_trinucleotides = rna_3_trinucleotides.sort_values(by="Count", ascending=False).head(35)
    max_count = max(df["Count"])

    fig, axs = plt.subplots(4, 2, figsize=(15, 20))

    plot_oligos(axs[0, 0], dna_5_dinucleotides, "5' DNA Dinucleotides", max_count)
    plot_oligos(axs[0, 1], dna_3_dinucleotides, "3' DNA Dinucleotides", max_count)
    plot_oligos(axs[1, 0], dna_5_trinucleotides, "5' DNA Trinucleotides", max_count)
    plot_oligos(axs[1, 1], dna_3_trinucleotides, "3' DNA Trinucleotides", max_count)
    plot_oligos(axs[2, 0], rna_5_dinucleotides, "5' RNA Dinucleotides", max_count)
    plot_oligos(axs[2, 1], rna_3_dinucleotides, "3' RNA Dinucleotides", max_count)
    plot_oligos(axs[3, 0], rna_5_trinucleotides, "5' RNA Trinucleotides", max_count)
    plot_oligos(axs[3, 1], rna_3_trinucleotides, "3' RNA Trinucleotides", max_count)

    plt.suptitle(f'{id} {prefix} Nucleotide Distribution Plots', fontsize=22, y=0.99)

    plt.tight_layout()
    plt.savefig(f"{prefix}.Nucl_distr.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot restr sites on read ends")
    parser.add_argument("id", help="Id")
    parser.add_argument("prefix", help="Prefix")
    args = parser.parse_args()

    main(args.id, args.prefix)