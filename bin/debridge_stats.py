import numpy as np
import matplotlib.pyplot as plt
import re
import argparse

def parse_count_matrix(file_content):
    lines = file_content.split('\n')
    index = lines.index('>>Count matrix') + 1
    lengths_index = lines.index('>>Count matrix lengths filter') + 1
    mtx = lines[index].strip("[]")
    mtx_lengths = lines[lengths_index].strip("[]")
    mtx1 = re.split('; | ', mtx)
    mtx1_lengths = re.split('; | ', mtx_lengths)
    return [int(item) for item in mtx1], [int(item) for item in mtx1_lengths]

def plot_data(ax, data, title, labels, ylim=None):
    ax.bar(range(len(data)), data, tick_label=labels)
    ax.set_title(title)
    ax.set_xlabel('Category')
    ax.set_ylabel('Count')
    if ylim is not None:
        ax.set_ylim(0, ylim)  # Set the y-axis limit

def main(se_filename, pe_filename, plot_filename):
    with open(se_filename, 'r') as file:
        se_file_content = file.read()

    with open(pe_filename, 'r') as file:
        pe_file_content = file.read()

    pe_matrix, pe_matrix_lengths = parse_count_matrix(pe_file_content)
    se_matrix, se_matrix_lengths = parse_count_matrix(se_file_content)
    max_val_pe = max(max(pe_matrix), max(pe_matrix_lengths))
    max_val_se = max(max(se_matrix), max(se_matrix_lengths))
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    labels_pe = ['00', '0F', '0R', '0M', 'F0', 'FF', 'FR', 'FM', 'R0', 'RF', 'RR', 'RM', 'M0', 'MF', 'MR', 'MM']
    labels_se = ['0', 'F', 'R', 'M']

    description = 'F - Forward, R - Reverse, 0 - Not Found, M - Double Bridge'
    plt.figtext(0.5, -0.05, description, ha='center', va='center', fontsize=12)

    plot_data(axs[0, 0], se_matrix, 'SE Matrix Counts', labels_se, max_val_se)
    plot_data(axs[1, 0], pe_matrix, 'PE Matrix Counts', labels_pe, max_val_pe)
    plot_data(axs[0, 1], se_matrix_lengths, 'SE Matrix Lengths Filter', labels_se, max_val_se)
    plot_data(axs[1, 1], pe_matrix_lengths, 'PE Matrix Lengths Filter', labels_pe, max_val_pe)

    plt.tight_layout()
    plt.savefig(plot_filename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot bridge summaries from provided files")
    parser.add_argument("se_file", help="The filename for single-end bridge summary")
    parser.add_argument("pe_file", help="The filename for paired-end bridge summary")
    parser.add_argument("save_plot", help="")

    args = parser.parse_args()

    main(args.se_file, args.pe_file, args.save_plot)