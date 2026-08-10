#!/usr/bin/env python3
"""Summary PDF for one or more ucarna_assembly.sh .ucaRNAs.tab outputs.

Usage:
    plot_ucarna_stats.py -o report.pdf --tab LABEL path/to/uca.ucaRNAs.tab [--tab LABEL2 path2 ...]

Multiple --tab entries sharing the same LABEL (e.g. two biological replicates
of the same experiment) are pooled into one group for the summary/faceted plots.
"""
import argparse
from collections import OrderedDict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import seaborn as sns

FIGSIZE = (11.7, 8.27)
PVAL_SIG = 0.05


def set_style() -> None:
    sns.set_style('white')
    sns.set_palette('husl')
    plt.rc('font', size=10)
    plt.rc('axes', labelsize=13, labelweight='bold', titlesize=15, titleweight='bold')
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9


def load_groups(entries):
    groups = OrderedDict()
    for label, path in entries:
        df = pd.read_csv(path, sep='\t')
        # A header-only (zero-ucaRNA) .tab has no rows to infer dtype from,
        # so pandas leaves numeric columns as object - coerce explicitly so
        # facet plots don't choke on log10/hist of an object-dtype column.
        numeric_cols = [c for c in df.columns if c in ('length', 'n_reads', 'p_value') or c.endswith('_TPM')]
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        groups.setdefault(label, []).append(df)
    return OrderedDict((label, pd.concat(dfs, ignore_index=True)) for label, dfs in groups.items())


def facet_grid(groups, plot_fn, title, pdf):
    n = len(groups)
    ncols = min(4, n) or 1
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=FIGSIZE, squeeze=False)
    for ax, (label, df) in zip(axes.flat, groups.items()):
        plot_fn(ax, label, df)
    for ax in axes.flat[n:]:
        ax.axis('off')
    fig.suptitle(title, fontsize=16, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig)
    plt.close(fig)


def plot_summary_counts(groups, pdf):
    labels = list(groups.keys())
    totals = [len(df) for df in groups.values()]
    sig = [int((df['p_value'] < PVAL_SIG).sum()) for df in groups.values()]

    fig, ax = plt.subplots(figsize=FIGSIZE)
    x = np.arange(len(labels))
    width = 0.35
    ax.bar(x - width / 2, totals, width, label='total ucaRNAs')
    ax.bar(x + width / 2, sig, width, label=f'p < {PVAL_SIG}')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=60, ha='right')
    ax.set_ylabel('# ucaRNAs')
    ax.set_title('Assembled ucaRNAs per experiment')
    ax.legend()
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_pvalue_facets(groups, pdf):
    def _plot(ax, label, df):
        ax.hist(df['p_value'], bins=30, range=(0, 1))
        ax.set_yscale('log')
        ax.set_title(label, fontsize=11)
        ax.set_xlabel('p-value')
        ax.set_ylabel('count')
    facet_grid(groups, _plot, 'p-value distribution', pdf)


def plot_length_facets(groups, pdf):
    def _plot(ax, label, df):
        ax.hist(np.log10(df['length']), bins=30)
        ax.set_title(label, fontsize=11)
        ax.set_xlabel('log10(length, bp)')
        ax.set_ylabel('count')
    facet_grid(groups, _plot, 'ucaRNA length distribution', pdf)


def plot_tpm_boxplot(groups, pdf):
    rows = []
    for label, df in groups.items():
        tpm_cols = [c for c in df.columns if c.endswith('_TPM')]
        if not tpm_cols:
            continue
        mean_tpm = df[tpm_cols].mean(axis=1)
        for v in mean_tpm:
            rows.append({'experiment': label, 'mean_TPM': v})
    if not rows:
        return
    tpm_df = pd.DataFrame(rows)
    tpm_df['log10_mean_TPM'] = np.log10(tpm_df['mean_TPM'] + 1e-3)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    order = list(groups.keys())
    sns.boxplot(data=tpm_df, x='experiment', y='log10_mean_TPM', order=order, ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=60, ha='right')
    ax.set_ylabel('log10(mean TPM across replicates + 1e-3)')
    ax.set_title('ucaRNA expression (StringTie TPM)')
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--out', required=True, help='output PDF path')
    ap.add_argument('--tab', nargs=2, action='append', metavar=('LABEL', 'PATH'),
                     required=True, help='label + path to a ucaRNAs.tab file, repeatable')
    args = ap.parse_args()

    set_style()
    groups = load_groups(args.tab)

    with PdfPages(args.out) as pdf:
        plot_summary_counts(groups, pdf)
        plot_pvalue_facets(groups, pdf)
        plot_length_facets(groups, pdf)
        plot_tpm_boxplot(groups, pdf)


if __name__ == '__main__':
    main()
