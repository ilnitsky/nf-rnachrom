import pyranges as pr
import pandas as pd
import numpy as np
import argparse
import warnings

from stats_voting import *
from io_.schemas import *
from io_.fileops import load_rdc, load_BED, make_pyranges_format, save_file_custom

warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Annotate RNA-parts of contacts')
parser.add_argument('rdc', type=str, help='Path to the RNA-DNA contact file')
parser.add_argument('annot', type=str, help='Path to the gene annotation')
parser.add_argument('--no_stat', action="store_true", help='Flag to skip statistics calculation')
parser.add_argument('--rna_parts', action='store_true', help='Whether file contains rna-parts')
parser.add_argument('--no_ribo', action="store_true", help='Flag to exclude contacts of rRNA')
parser.add_argument('--cpus', nargs='?', type=int, default = 1)
parser.add_argument('--outdir', nargs='?', type=str, default = './results', help='A folder to save results')
args = parser.parse_args()


def run_annotation_and_voting(args):
    """
    Annotate RNA-parts of contacts, apply a voting procedure and remove
    contacts derived from ribosomal DNA if specified.
    
    """

    d = {}
    save_names = {'singletons':       singleton_suffix,
                  'selected_annot':   selected_ann_suffix,
                  'voted':            voted_suffix, 
                  'complement_annot': complement_ann_suffix,
                  'voted_noribo':     voted_nr_suffix
                  }
    empty = []
    
    gene_ann = load_BED(args.annot).rename(columns={'name': 'gene_name', 'source': 'Source'})
    gene_ann = make_pyranges_format(gene_ann, strand = True, save_old_names = False)


    gene_ann.gene_length = gene_ann.lengths()
    gene_ann = gene_ann[['gene_name', 'gene_type', 'Source', 'gene_length']]

    
    if args.rna_parts:
        coln_rdc = rdc_BED_sample[:5] + [rdc_BED_sample[-1]]
        coln_voted = voted_BED[:5] + voted_BED[9:13]
        cnts = load_rdc(args.rdc, header = 0, names = coln_rdc, dtype = rdc_dtypes, use_cols = [0,1,2,3,4,6])

    else:
        coln_rdc = rdc_BED_sample.copy()
        coln_voted = voted_BED
        cnts = load_rdc(args.rdc, header = 0, names = coln_rdc, dtype = rdc_dtypes)

    cnts.drop(cnts.filter(regex='_cigar').columns, axis=1, inplace=True)
    cnts, old_cols = make_pyranges_format(cnts, strand = True)

    d['selected_annot'], d['complement_annot'] = annotate_rdc(cnts, gene_ann, cpus = args.cpus)
    del gene_ann

    d = dict(map(lambda item: (item[0], item[1].drop(like="annot$").as_df()), d.items()))
    save_file_custom(d['selected_annot'].iloc[:,:len(coln_voted)], args.outdir, args.rdc, save_names['selected_annot'], hdr=coln_voted)
    d['voted'] = vote(d['selected_annot'])
    del d['selected_annot']

    cnts = cnts.as_df()
    cnts.colnames = old_cols

    d['singletons'] = cnts[~cnts['id'].isin(d['voted']['id'])]
    
    if args.no_ribo:
        d['voted_noribo'] = remove_ribo(d['voted'])

    for name, dfr in d.items():
        if not dfr.empty:
            if name == 'singletons':
                save_file_custom(dfr, args.outdir, args.rdc, save_names[name], hdr=coln_voted[:-3])
            else:
                save_file_custom(dfr.iloc[:,:len(coln_voted)], args.outdir, args.rdc, save_names[name], hdr=coln_voted)
                    
        else:
            save_file_custom(dfr, args.outdir, args.rdc, save_names[name], hdr=coln_voted)
            empty.append(name)

    for empty_df in empty:
        d.pop(empty_df, None)

    if not args.no_stat:
        calculate_stats(d, args.rdc, args.outdir)

    return 



def annotate_rdc(contacts: pr.PyRanges, annot: pr.PyRanges, cpus: int) -> pr.PyRanges:
    """
    Create annotated genomic intervals from a RNA-DNA contacts file.

    Parameters
    ----------
    contacts: PyRanges
        RDC-like file converted to PyRanges
        
    annot: PyRanges
        A GTF- or BED-like file converted to PyRanges

    Returns
    -------
    PyRanges: 
         1. Intervals annotated by genes on selected strand (selected_annot)
         2. Intervals annotated by genes on the complementing strand (complement_annot)

    """

    
    all_annot = contacts.join(annot, how = 'left', strandedness = False, suffix = '_annot', nb_cpu = cpus)
    print(contacts.head())

    strand_match = all_annot.Strand == all_annot.Strand_annot
    gene_found = all_annot.gene_name != '-1'
    
    return all_annot[(strand_match & gene_found)], all_annot[(~strand_match & gene_found)]


def vote(an_contacts: pd.DataFrame) -> pd.DataFrame:
    """
    Apply voting procedure to RNA parts of contacts. In the case of unambiguous annotation 
    the preference is given to the gene with more dense coverage of RNA-parts.

    """
    an_contacts['count'] = an_contacts['gene_name'].map(an_contacts['gene_name'].value_counts())
    an_contacts['density'] = (an_contacts['count'] / 
                                 an_contacts['gene_length']) * 100000

    if any('gencode' in s for s in an_contacts['Source'].unique()):
        an_contacts['gencode'] = np.where(an_contacts['Source'].str.contains('gencode'), 'gencode', 'other')
        return an_contacts.sort_values(['density', 'gencode'], ascending=(False, True)).drop_duplicates('id', keep='first').drop('gencode', axis=1)
    
    else:
        return an_contacts.sort_values('density', ascending=False).drop_duplicates('id', keep='first')



def remove_ribo(contacts: pd.DataFrame, ribo_list = ['rRNA', 'rRNA_pseudogene', 'rRNA_RepM']) -> pd.DataFrame:
    """
    Removes contacts mapping to ribosomal DNA 
    
    """
    return contacts[~contacts.gene_type.isin(ribo_list)]

run_annotation_and_voting(args)
