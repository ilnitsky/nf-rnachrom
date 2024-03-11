import pandas as pd
import numpy as np
import argparse

pd.options.mode.chained_assignment = None

from io_.fileops import *
from io_.schemas import *
from BaRDIC.bardic.api.schemas import BedValidator

parser = argparse.ArgumentParser(description='Validates RNA annotation file')
parser.add_argument('annot', type=str, help='Path to RNA annotation file')
parser.add_argument('chrsizes', type=str, help='Path to file with chromosome sizes')
parser.add_argument('--annot_format', nargs='?', default='GTF', type=str, choices = ['GTF', 'BED'], help='File format: GTF or BED (default: GTF)')
parser.add_argument('--outdir', nargs='?', type=str, default = './results', help='')
args = parser.parse_args()

types = {
         'gencode': 1,
         'GB_snomirna': 0,
         'GB_repM': 0,
         'vlinc': 0,
         'GB_trna': 0,
         'Xrna_human': 1,
        }

def run_validation(args):

    if args.annot_format == 'GTF':
        gene_ann = pr.read_gtf(args.annot)

    elif args.annot_format == 'BED':
        gene_ann = load_BED(args.annot)
    else:
        raise Exception('Gene annotation format is not supported')
        
    chrsizes = load_rdc(args.chrsizes,
                   names = chrsizes_format).set_index(chrsizes_format[0]).to_dict()[chrsizes_format[1]]

    gene_ann = zero_based_to_one(gene_ann.drop_duplicates(), types)
    
    # save initial names in separate column
    gene_ann['name_unmod'] = gene_ann.loc[:, 'name']
    gene_ann = modify_dupl_names(gene_ann, 'name')
    save_file_custom(gene_ann, args.outdir, args.annot, '.test_modify.tab')
    
    good_annot, bad_annot = BedValidator(chrsizes).validate_df(gene_ann.assign(score = '.'))
    
    # Save BED6 for BaRDIC
    save_file_custom(good_annot[BED6], args.outdir, args.annot, annot_bed6_suffix)
    # Save in pipeline format
    save_file_custom(good_annot[annot_BED], args.outdir, args.annot, annot_suffix)
    # Save unmodified names
    save_file_custom(good_annot, args.outdir, args.annot, annot_old_suffix)
    # Save erroneous strings
    save_file_custom(bad_annot[annot_BED], args.outdir, args.annot, bad_annot_suffix)
    
    return
    
def zero_based_to_one(df, dict_coord_system):
    """
    Convert zero-based entries to one-based.
    """
    # when not specified, source is assumed to be one-based
    df['coord_system'] = df['source'].map(dict_coord_system).fillna(1).astype(int)
    df.loc[df.coord_system == 0,["start"]] += 1
    return df

def modify_dupl_names(df: pd.DataFrame, col: str) -> pd.DataFrame:
    """
    Adds suffix to gene name to avoid duplicated entries.
    """
    df[col] = df[col] + '__' + df.groupby(col).cumcount().astype(str)
    df[col] = df[col].map(lambda x: x.replace('__0', ''))
    df[col] = df[col].str.replace('/','_')
    return df


run_validation(args)
    
    
    
    