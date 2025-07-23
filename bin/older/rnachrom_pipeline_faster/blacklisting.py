import pandas as pd
import argparse
import pyranges as pr

from io_.fileops import *
from io_.schemas import *
from natsort import natsort_keygen

parser = argparse.ArgumentParser(description='Removes contacts whose DNA parts overlap with blacklisted regions')
parser.add_argument('rdc', type=str, help='Path to the RNA-DNA contact file')
parser.add_argument('blacklist', help='Path to the BED file')
parser.add_argument('--outdir', type=str, nargs='?', default = './results', help='')
args = parser.parse_args()


def remove_blacklisted_regions(args):

    """
    Removes contacts whose DNA parts overlap with blacklisted regions
    
    """
    cnts_orig = load_rdc(args.rdc, header = 0, names = rdc_BED, dtype = rdc_dtypes)
    cnts_dna_first = swap_rna_dna_order(cnts_orig) #put dna contacts first
    del cnts_orig
    cnts_py, old_names = make_pyranges_format(cnts_dna_first)
    del cnts_dna_first
    
    blist = make_pyranges_format(load_BED(args.blacklist, names = BED3), save_old_names = False)
    blck_removed = cnts_py.overlap(blist, invert = True).as_df()
    blck_removed.columns = old_names #restore original column names

    calculate_stats(cnts_py.as_df(), blck_removed)

    save_file_custom(blck_removed[rdc_BED], args.outdir, args.rdc, suff = blacklist_suffix)
    
    return 


def calculate_stats(before_blckl: pd.DataFrame, after_blckl: pd.DataFrame):
    """
    Returns file with statistics on blacklisted regions.
    
    """
    stats = {}

    for chrom in before_blckl['rna_chr'].unique():
        stats[chrom] = {}
        stats[chrom]['raw'] = before_blckl.query('rna_chr == @chrom').shape[0]
        stats[chrom]['after_blacklisting'] = after_blckl.query('rna_chr == @chrom').shape[0]
    
    save_file_custom(pd.DataFrame.from_dict(stats, orient='index').reset_index(names = "rna_chr").sort_values(by="rna_chr", key=natsort_keygen()), 
                     args.outdir, args.rdc, suff = blacklist_suffix_stat)
    return 
    

remove_blacklisted_regions(args)
