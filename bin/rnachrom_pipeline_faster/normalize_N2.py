import argparse
import pandas as pd
import warnings

from natsort import natsort_keygen, natsorted
from io_.fileops import *
from io_.schemas import *

warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Calculate N2 background normalization')
parser.add_argument('rdc', type=str, help='Path to the RNA-DNA contact file')
parser.add_argument('chrsizes', type=str, help='Path to the file with chromosome sizes')
parser.add_argument('--by_chr', type=str, nargs='?', default = "", help='Path to raw N2 statistics if files are split by DNA chromosome')
parser.add_argument('--outdir', type=str, nargs='?', default = './results')
args = parser.parse_args()


def re_calculate(args):

    rdc = load_rdc(args.rdc, header = 0, names = normalized_BED_raw, dtype = rdc_dtypes)
    chrsizes = load_rdc(args.chrsizes, names = ['chr', 'len']) 

    if args.by_chr:
        stats_chr = load_rdc(args.by_chr, header = 0)
        library_size = stats_chr['n_contacts'].sum() 
        sum_of_weights = stats_chr['N2_raw_sum'].sum()

    else:
        library_size = rdc.shape[0]
        sum_of_weights = rdc["N2_raw"].sum()

    rdc["N2"] = round(rdc["N2_raw"] * (library_size / sum_of_weights), 3)
    
    save_file_custom(rdc[normalized_BED], args.outdir, args.rdc, N2_suffix)
    save_file_custom(rdc.groupby(['rna_chr','dna_chr']).agg(n_contacts = ('gene_name', 'count'),
                                                       N2_sum=('N2', 'sum')).sort_values(['dna_chr'], ascending = True, key=natsort_keygen()),
                                                       args.outdir, args.rdc, N2_final_stats, idx = True)
        
    return

re_calculate(args)
