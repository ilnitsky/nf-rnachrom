import argparse
import pandas as pd
import numpy as np
import pyranges as pr
import warnings

from io_.fileops import *
from io_.schemas import *

#from BaRDIC.bardic.api.schemas import BedValidator, ContactsValidator
#from bardic.cli import bardic_parser, bardic_subparsers

warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None

parser_input = argparse.ArgumentParser(
                    prog='input',
                    description='Prepares background file')

parser_input.add_argument('--bgtype', type=str, nargs='?', default = 'mrnas', choices = ('input_contacts', 'rnas_list', 'mrnas'), help='Type of backround. If "mrnas", then will calculate endogeneous background from mRNAs trans contacts extracted from RNA-DNA file (bgdata parameter is not needed). If "input_contacts", will create bedgraph track from bed3 provided in "bgdata". If "rnas", then will calculate background from trans contacts of RNAs supplied in file as "bgdata".  (default: "mrnas")')
parser_input.add_argument('--bgdata', type = str, nargs='?', help='A file with data on background. If --bgtype="rnas_list", this is a file with names of RNAs with one RNA name per line. '
'If --bgtype="custom", this is a bed3 file with input contacts.')

parser_input_group = parser_input.add_argument_group('Input')
parser_input_group.add_argument('rdc_file', type=str, help='Path to the RNA-DNA contact file')
parser_input_group.add_argument('annot', type=str, help='Path to the file with gene annotation')
parser_input_group.add_argument('chrsizes', type=str, help='Path to the file with chromosome sizes')
parser_output_group = parser_input.add_argument_group('Output')
parser_output_group.add_argument('--outdir', type=str, nargs='?', default = './results', help='A folder in which to store the results')


def prepare_bg(rdc_file, annot, chrsizes, bgdata, bgtype, outdir):
    """
    Runs validation steps and prepare backgroud track.

    """
    
    chrsz = load_chrsizes_as_pyranges(chrsizes)

    ann = load_BED(annot, names = annot_BED)
    
    rdc = load_rdc(rdc_file, header = 0, names = voted_BED,
               use_cols = ['dna_chr', 'dna_bgn','dna_end', 'gene_name', 'rna_chr'])

    genome_b = genome_binning(chrsz)


    dna_bed3 = ['dna_chr', 'dna_bgn', 'dna_end']
    if bgtype == 'mrnas':
        pc = ann.query("gene_type == 'protein_coding'")['name'].unique().tolist()
        bg_bedgraph = background_track(rdc.query('gene_name.isin(@pc) & dna_chr != rna_chr')[dna_bed3], genome_b)
    
    elif bgtype == 'rnas_list':
        rnas = pd.read_csv(bgdata, names = ['name'])['name']
        bg_bedgraph = background_track(rdc.query('gene_name.isin(@rnas) & dna_chr != rna_chr')[dna_bed3], genome_b)

    elif bgtype == 'input_contacts':
        bg_cnts = load_BED(bgdata, names = dna_bed3)
        bg_bedgraph = background_track(bg_cnts, genome_b)

    save_file_custom(bg_bedgraph.as_df(), outdir, rdc_file, bardic_bg, hdr = False)
    
    rdc_inp = prepare_input(rdc)
    save_file_custom(rdc_inp, outdir, rdc_file, bardic_input_rdc, hdr = False)

    #good_annot, bad_annot, good_rdc, bad_rdc = validate_input(ann_inp, rdc_inp, chrsz)
    #save_file_custom(good_annot, outdir, annot, bardic_input_annot, hdr = False)
    #save_file_custom(bad_annot, outdir, annot, bardic_bad_annot)
    #save_file_custom(good_rdc, outdir, rdc_file, bardic_input_rdc, hdr = False)
    #save_file_custom(bad_rdc, outdir, rdc_file, bardic_bad_rdc)
    return

parser_input.set_defaults(func=prepare_bg)   

def genome_binning(chrsizes_p: pr.PyRanges, binsize = 1000) -> pr.PyRanges:
    """
    Partitions all chromosomes into fixed-size bins.

    Parameters
    ----------
    chrsizes_p : pr.PyRanges
        Chromosome lengths in pr.PyRanges format
    
    binsize : int
        Bin width (default: 1 Kb)

    """   
    genome_binned = pr.genomicfeatures.tile_genome(chrsizes_p, binsize)
    
    return genome_binned


def background_track(bg_cnts: pd.DataFrame, binned_genome: pr.PyRanges) -> pr.PyRanges:
    """
    Count background contacts in bins.
    
    """
    return binned_genome.count_overlaps(make_pyranges_format(bg_cnts, save_old_names = False), overlap_col="Count")

def prepare_input(rdc: pd.DataFrame, annot: pd.DataFrame):
    """
    Convert RDC and RNA annotation files to BED6 format.
    
    """
    rdc_input = rdc.drop(columns = ['rna_chr'], inplace = False).assign(score='.', strand = '.')
    rdc_input.columns = BED6

    #annot_input = annot[annot_BED[:-2]].assign(score='.')
    
    return rdc_input #, annot_input[BED6]
    


# def validate_input(rna_annot: pd.DataFrame, rdc: pd.DataFrame, chrsizes: pr.PyRanges):
#     """
#     Checkes bad lines in RDC file and RNA gene annotation file.

#     """
#     chrsizes_dict = chrsizes.as_df().set_index('Chromosome').to_dict()['End']
#     g_annot, b_annot = BedValidator(chrsizes_dict).validate_df(rna_annot)
#     print(rna_annot[rna_annot['name'] == 'RNU5E-1'].head())
#     print(rdc.head())
#     print(dict(zip(g_annot['name'], [None]*len(g_annot['name']))))
#     g_rdc, b_rdc = ContactsValidator(chrsizes_dict, dict(zip(g_annot['name'], [None]*len(g_annot['name'])))).validate_df(rdc)
#     return g_annot, b_annot, g_rdc, b_rdc

def main():
    args = parser_input.parse_args()
    #print(args)
    func = args.func
    del args.func
    kwargs = vars(args)
    func(**kwargs)


if __name__ == "__main__":
    main()
