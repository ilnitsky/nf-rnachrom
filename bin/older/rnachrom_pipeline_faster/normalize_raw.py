import argparse
import pandas as pd
import numpy as np
import warnings

from natsort import natsort_keygen

from io_.fileops import *
from io_.schemas import *

warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Calculate raw N2 background normalization')
parser.add_argument('rdc', type=str, help='Path to the RNA-DNA contact file')
parser.add_argument('bg_sm', type=str, help='Path to smoothed background')
parser.add_argument('--cpus', type=int, nargs='?', default = 1)
parser.add_argument('--outdir', type=str, nargs='?', default = './results', help='')
args = parser.parse_args()

def run_N2raw(args, how = 'bgn'):

    rdc = load_rdc(args.rdc, header = 0, names = voted_BED, dtype = rdc_dtypes)


    sm_bg = load_BED(args.bg_sm, names = BED_Stereogene)
    sm_bg = make_pyranges_format(sm_bg, save_old_names = False)

    rdc = backannot(sm_bg, rdc, cpus = args.cpus)
    calculate_N2_normalization(rdc)

    return


def backannot(bg_smoothed: pr.PyRanges,
              rdc_: pd.DataFrame, cpus = 1, how = "bgn"):
    """
    Prescribe background weights to contacts.
    
    """
    if how == 'mean':
        rdc_ = rdc_.assign(dna_bgn_bg = lambda x: (x.dna_bgn + x.dna_end)/2,
                  dna_end_bg = lambda x: x.dna_bgn_bg + 1 )
    else:
        rdc_ = rdc_.assign(dna_bgn_bg = lambda x: x.dna_bgn,
                  dna_end_bg = lambda x: x.dna_bgn_bg + 1 )

    dna_first, names = make_pyranges_format(swap_rna_dna_order(rdc_, suffix_dna = ['chr', 'bgn_bg', 'end_bg']))
    dna_first = dna_first.join(bg_smoothed, suffix = '_bg', how = 'left', nb_cpu = cpus)
    dna_first = dna_first.drop(like="bg$").as_df()
    names.append('bg_sm') 
    dna_first.columns = names
    return dna_first[normalized_BED[:-1]]


def calculate_N2_normalization(rdc_back: pd.DataFrame):
    
    """
    Calculate N2 metrics and renormalize weights to preserve the total number of RNA-DNA contacts.
    
    """

    rdc_back['N2_raw'] = np.where(rdc_back["bg_sm"] != -1, 1 / (rdc_back["bg_sm"] + 0.5), 0)

    save_file_custom(rdc_back.groupby(['rna_chr','dna_chr']).agg(n_contacts = ('gene_name', 'count'), 
                                                       N2_raw_sum=('N2_raw', 'sum')).sort_values(['dna_chr'], ascending = True, key=natsort_keygen()),
                                                       args.outdir, args.rdc, N2_raw_stats, idx = True)
    
    save_file_custom(rdc_back[normalized_BED_raw], args.outdir, args.rdc, N2_raw_suffix)
    return


run_N2raw(args)
