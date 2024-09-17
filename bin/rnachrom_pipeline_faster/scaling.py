import argparse
import pandas as pd
import numpy as np
from natsort import natsort_keygen
from pandas.api.types import union_categoricals

pd.options.mode.chained_assignment = None

from io_.schemas import *
from io_.fileops import load_rdc, load_BED, save_file_custom


parser = argparse.ArgumentParser(description='Calculate inner- and inter-chromosomal scaling based on mRNA contacts density')
parser.add_argument('rdc', type=str, help='Path to the RNA-DNA contact file, background_normalized')
parser.add_argument('chrsizes', type=str, help='Path to the file with chromosome sizes')
parser.add_argument('--factor', type=float, nargs='?', default = 1.2, help='Scale factor')
parser.add_argument('--outdir', type=str, nargs='?', default = './results')
args = parser.parse_args()


def run_all_scaling(args, b0 = 0, b_n = 100, factor = 1.2, threads = 1):

    cnts_all = load_rdc(args.rdc, header = 0, names = normalized_BED, dtype = rdc_dtypes)
    
    union = union_categoricals([cnts_all.rna_chr, cnts_all.dna_chr]).categories
    cnts_all['rna_chr'] = cnts_all.rna_chr.cat.set_categories(union)
    cnts_all['dna_chr'] = cnts_all.dna_chr.cat.set_categories(union)

    chrsizes = load_rdc(args.chrsizes, names = chrsizes_format)
    bins_ = make_geom_bins(length = max(chrsizes.len), factor = args.factor)
    
    #inner scaling
    scaled_cis = run_inner_scaling(bins_, cnts_all)

    if not scaled_cis.empty:
        print(scaled_cis[['rna_chr', 'bin', 'n_PC_bin', 'dens_PC_bin']])
        save_file_custom(scaled_cis[['rna_chr', 'bin', 'n_PC_bin', 'dens_PC_bin']].sort_values(['rna_chr', 'bin'], ascending = True, key=natsort_keygen()).drop_duplicates(),
       	             args.outdir, args.rdc, PC_inner_suffix)
        scaled_cis = scaled_cis[normalized_BED_scaling]

    save_file_custom(scaled_cis, args.outdir, args.rdc, scaling_inner_suffix)
    del scaled_cis
    
    scaled_trans = run_inter_scaling(chrsizes, cnts_all)

    if not scaled_trans.empty:   
        save_file_custom(scaled_trans[['rna_chr', 'dna_chr', 'n_PC_trans', 'dens_PC']].sort_values(['rna_chr', 'dna_chr'], ascending = True, key=natsort_keygen()).drop_duplicates(),
                     args.outdir, args.rdc, PC_inter_suffix)
        scaled_trans = scaled_trans[normalized_BED_scaling]
    save_file_custom(scaled_trans[normalized_BED_scaling], args.outdir, args.rdc, scaling_inter_suffix)
    
    return 



def make_geom_bins(length: int, factor: float, b0 = 0, b_n = 100) -> list:
    """
    Make geometrically increasing bins

    Parameters
    ----------
    length : int
        Chromosome length
        
    b0 :
        the 1st term of a geometric sequence
    
    b_n :
        the first step of a geometric sequence
    
    factor:
         the common ratio
    
    """
    
    bins = [0]
    
    while b_n < length:
        bins.append(b_n)
        b_n = b0 + factor * b_n    
    bins.append(length)
    
    return bins


def rna_dna_distance_binning(cnts: pd.DataFrame, bins: list) -> pd.DataFrame:
    
    """
    Assign bin according to distance between RNA and DNA parts of a contact.

    """
    cnts = cnts.query("dna_chr == rna_chr")
    
    cnts['dist'] = (cnts['rna_bgn'] - cnts['dna_bgn']).abs()
    
    cnts['bin'] = pd.cut(cnts['dist'], bins = bins, include_lowest = True, 
                   precision = 2, duplicates = 'drop')

    cnts['bin_size'] = pd.to_numeric(cnts['bin'].apply(lambda x: (x.right - x.left)))
    cnts['bin_mean'] = pd.to_numeric(cnts['bin'].apply(lambda x: x.mid))
    
    return cnts


def inner_PC_density(cnts_b: pd.DataFrame) -> pd.DataFrame:
    
    """
    Calculate protein-coding RNA contacts density in each bin, based on RNA-DNA distance.

    """
    cnts_b = cnts_b.query('gene_type == "protein_coding"')
    libS_pc = cnts_b['N2'].sum()
        
    cnts_b['n_PC_bin'] = cnts_b.groupby(by = ['bin'])['N2'].transform('sum')
    cnts_b['dens_PC_bin'] = (cnts_b['n_PC_bin'] / cnts_b['bin_size']) / libS_pc
        
    cnts_b['log10Dist_PC'] = np.log10(cnts_b['bin_mean'])
    cnts_b['log10Dens_PC'] = np.log10(cnts_b['dens_PC_bin'])
    bins_pc_stat = cnts_b.sort_values('bin').drop_duplicates('bin', keep='first')
        
    return bins_pc_stat[['bin', 'n_PC_bin', 'dens_PC_bin','log10Dist_PC', 'log10Dens_PC']]


def inter_PC_density(cnts: pd.DataFrame, chrsizes: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate protein-coding RNA contacts density on non-parental chromosomes

    """
    
    #select only mRNA trans-contacts
    cnts = cnts.query('dna_chr != rna_chr and gene_type == "protein_coding"')
    
    cnts['n_PC_trans'] = cnts.groupby(by = ['rna_chr', 'dna_chr'])['N2'].transform('sum')
    cnts = cnts.merge(chrsizes, left_on = 'dna_chr', right_on = 'chr')
    cnts['dens_PC'] = (cnts['n_PC_trans'] / cnts['len']) * 1000
    
    cnts.drop(columns=['chr'], inplace = True)
    
    bins_pc_stat = cnts.sort_values('dna_chr').drop_duplicates('dna_chr', keep='first')
    
    return bins_pc_stat[['dna_chr', 'n_PC_trans', 'dens_PC']].reset_index(drop=True)


def inner_scaling(cnts_cis: pd.DataFrame, dens_pc_cis: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate scaling for all RNA-DNA cis-contacts.
        
    """
    cnts_cis = cnts_cis.query('rna_chr == dna_chr')
    
    if not cnts_cis.empty:
        libS_cis = cnts_cis['N2'].sum()
    
        cnts_cis = cnts_cis.merge(dens_pc_cis, on = 'bin')
    
        cnts_cis['SCnorm_raw'] = (cnts_cis['N2'] / cnts_cis['dens_PC_bin']) 
        cnts_cis['SCnorm'] = cnts_cis['SCnorm_raw'] * (libS_cis / cnts_cis['SCnorm_raw'].sum())
        cnts_cis.drop(columns=['SCnorm_raw'], inplace = True)
    
    return cnts_cis


def inter_scaling(cnts_trans: pd.DataFrame, dens_pc_trans: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate scaling for all RNA-DNA trans-contacts.
        
    """
    cnts_trans = cnts_trans.query('rna_chr != dna_chr')
    cnts_trans = cnts_trans.merge(dens_pc_trans, on = 'dna_chr')
    
    cnts_trans['SCnorm_raw'] = (cnts_trans['N2'] / cnts_trans['dens_PC'])
    cnts_trans['libS_trans'] = cnts_trans.groupby('rna_chr')['N2'].apply(lambda x: x.sum())
    cnts_trans['SCnorm'] = cnts_trans['SCnorm_raw'] * (cnts_trans['libS_trans'] / cnts_trans['SCnorm_raw'].sum())
    
    cnts_trans.drop(columns=['SCnorm_raw', 'libS_trans'], inplace = True)
    
    return cnts_trans


def run_inner_scaling(bins: list, cnts_all: pd.DataFrame) -> pd.DataFrame:
    """
    Run all commands to calculate inner scaling.
        
    """
    cnts_inner_all = rna_dna_distance_binning(cnts_all, bins)
    #dataset with cis-contacts only, RNA-DNA distances are binned
    
    pc_values = inner_PC_density(cnts_inner_all)    
    cnts_inner_scaled = inner_scaling(cnts_inner_all, pc_values)
    
    return cnts_inner_scaled


def run_inter_scaling(chrsizes: pd.DataFrame, cnts_all: pd.DataFrame) -> pd.DataFrame:
    """
    Run all commands to calculate inter-chromosomal scaling.
        
    """
    pc_values = inter_PC_density(cnts_all, chrsizes).sort_values(['dna_chr'], ascending = True, key=natsort_keygen())
    cnts_inter_scaled = inter_scaling(cnts_all, pc_values)
    
    return cnts_inter_scaled


run_all_scaling(args)