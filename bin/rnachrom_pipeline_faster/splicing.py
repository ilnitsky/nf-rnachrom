import pandas as pd
import numpy as np
import re
import argparse
import warnings

from natsort import natsort_keygen
from io_.fileops import load_rdc, save_file_custom
from io_.schemas import rdc_BED, rdc_dtypes, CIGAR_suffix, CIGAR_suffix_stat


warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Parses CIGAR field and processes spliced contacts')
parser.add_argument('rdc', type=str, help='Path to the primary rna-dna contacts file')
parser.add_argument('--rna_parts', action='store_true', help='Whether file contains dna-parts')
parser.add_argument('--no_header', action='store_true', help='Headerless flag')
parser.add_argument('--outdir', type=str, nargs='?', default = './results', help='A folder to save the results')



def run_cigar_processing(args):
    """
    Parse CIGAR field and processes spliced contacts: 
    (1) remove complex splicing cases ('I', 'D', > 1 'N' in CIGAR string)
    (2) retain the longer, ungapped part of the read when read contains one spanning splicing even (1 'N' in CIGAR field) 
    
    """
    if not args.rna_parts:
        coln_rdc = rdc_BED.copy()
    else:
        coln_rdc = rdc_BED.copy()[:6]

    if not args.no_header:
        rdc = load_rdc(args.rdc, header = 0, names = coln_rdc, dtype = rdc_dtypes)
    else:
        rdc = load_rdc(args.rdc, header = None, names = coln_rdc, dtype = rdc_dtypes)

    rdc = splicing_cases(rdc)
    calculate_stats(rdc, args.outdir)
    
    if not (rdc['N_cnt'] == 0).all() & (rdc['ID_cnt'] == 0).all():
        rdc = process_simple_splicing(rdc)
    save_file_custom(rdc[coln_rdc], args.outdir, args.rdc, CIGAR_suffix)
    
    return


def splicing_cases(df: pd.DataFrame) -> pd.DataFrame:
    """
    Count the N, I, D occurences in the CIGAR string.
    
    """
    df['N_cnt'] = df.rna_cigar.str.count('N')
    df['ID_cnt'] = df.rna_cigar.str.count(r'[ID]')
    return df


def process_simple_splicing(df: pd.DataFrame) -> pd.DataFrame: 

    """
    Keep longer part of spliced contacts and re-calculates start/end coordinates.
    
    """

    df_ss = df[(df.N_cnt == 1) & (df.ID_cnt == 0)]

    df_ss[['N1','N2']] = df_ss['rna_cigar'].str.split(r'[NM]', expand=True).iloc[:,[0,2]]
    df_ss = df_ss.astype({'N1':'int','N2':'int'})
    
    df_ss['rna_end'] = np.select([df_ss.N1 >= df_ss.N2], [df_ss.rna_bgn + df_ss.N1], default = df_ss.rna_end)
    df_ss['rna_bgn'] = np.select([df_ss.N1 < df_ss.N2], [df_ss.rna_end - df_ss.N2], default = df_ss.rna_bgn)
    df_ss.drop(['N1', 'N2'], axis = 1, inplace = True)
    
    return pd.concat([df[(df.N_cnt == 0) & (df.ID_cnt == 0)], df_ss])


def calculate_stats(df, outdir): 

    """
    Save file with statistics on splicing.
    
    """
    stats = {}
    
    for chrom in df['rna_chr'].unique():
        stats[chrom] = {}
        df_chr = df[df['rna_chr'] == chrom]
        stats[chrom]['raw'] = df_chr.shape[0]
        stats[chrom]['no_splicing'] = df_chr[(df_chr.N_cnt == 0) & (df_chr.ID_cnt == 0)].shape[0]
        stats[chrom]['splicing_corrected'] = df_chr[(df_chr.N_cnt == 1) & (df_chr.ID_cnt == 0)].shape[0]
        stats[chrom]['splicing_removed'] = stats[chrom]['raw'] - stats[chrom]['no_splicing'] - stats[chrom]['splicing_corrected']
        stats[chrom]['splicing_out'] = stats[chrom]['no_splicing'] + stats[chrom]['splicing_corrected']
        del df_chr
        
    save_file_custom(pd.DataFrame.from_dict(stats, orient='index').reset_index(names = "rna_chr").sort_values(by="rna_chr",      key=natsort_keygen()), outdir, args.rdc, CIGAR_suffix_stat)
    return 

args = parser.parse_args()
run_cigar_processing(args)
