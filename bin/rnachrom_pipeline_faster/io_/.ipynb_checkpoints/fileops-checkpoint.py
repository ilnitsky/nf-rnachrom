import pyranges as pr
import pandas as pd
import os
import pathlib
from pathlib import Path

from .schemas import annot_BED, BED3, BED6


def load_BED(path: str, names = annot_BED, dtype = None) -> pd.DataFrame:
    """
    Read intervals and metadata from a bed and bed-like files.

    """
    if names == BED3:
        return pd.read_csv(path, header = None, index_col=None, sep = '\t', usecols = [0,1,2], names = names, dtype = dtype)
    elif names == BED6:
        return pd.read_csv(path, header = None, index_col=None, sep = '\t', usecols = [0,1,2,3,5,6], names = names, dtype = dtype)
    else:
        return pd.read_csv(path, header = 0, index_col=None, sep='\t', names = names, dtype = dtype)

    

def load_rdc(path: str, header = None, names = [], use_cols = False, dtype = None) -> pd.DataFrame:
    """
    Read intervals and metadata from a RNA-DNA contacts file.

    """
    if header == 0 and len(names) == 0:
        
        rdc = pd.read_csv(path, header=header, index_col=False, sep='\t', usecols = [0,1,2,3,5,6], dtype = dtype)
    elif not use_cols:
        rdc = pd.read_csv(path, header=header, index_col=False, sep='\t', names=names, dtype = dtype)
    else:
        rdc = pd.read_csv(path, header=header, index_col=False, sep='\t', usecols = [0,1,2,3,5,6], names=names, dtype = dtype)
    return rdc



def load_chrsizes_as_pyranges(chrsizes_path: str) -> pr.PyRanges:
    """
    Convert file with chromosome sizes to PyRanges object.

    """
    
    chrsizes = pd.read_csv(chrsizes_path, header = None, sep = '\t')
    chrsizes.insert(loc = 1, column = 'Start', value = 0)
    chrsizes.columns = BED3
    
    return make_pyranges_format(chrsizes, save_old_names = False)


def swap_rna_dna_order(rdc: pd.DataFrame, suffix_rna = ['chr', 'bgn', 'end'],
                                                          suffix_dna = ['chr', 'bgn', 'end']) -> pd.DataFrame:
    """
    Swap the order of DNA and RNA-parts of contacts.

    Parameters
    ----------
    rdc : DataFrame
    
    first : str
        Which contacts go first, column prefix.
        
    after: str
        Which contacts go after, column prefix.
        
    suffix: list
        Which columns to swap.
    

    Returns
    -------
    DataFrame

    """
    
    def contacts_first(rdc: pd.DataFrame, column_ord: list, suf_rna: list, suf_dna: list) -> list:
        
        n1 = "_".join(['rna', suf_rna])
        n2 = "_".join(['dna', suf_dna])

        x, y = column_ord.index(n1), column_ord.index(n2)
        column_ord[y], column_ord[x] = column_ord[x], column_ord[y]
        
        return column_ord
    
    order = list(rdc.columns)
    for s_rna, s_dna in zip(suffix_rna, suffix_dna):
        order = contacts_first(rdc, order, s_rna, s_dna)
        
    return rdc[order]


def make_pyranges_format(rdc_: pd.DataFrame, name = False, strand = False, save_old_names = True)  -> pr.PyRanges:
    """
    Make RDC column names suitable for PyRanges format.

    """
    coln = list(rdc_.columns)
    if ('chr' in coln[0]) & ('bgn' in coln[1] or 'start' in coln[1]) & ('end' in coln[2]):
        rdc_.rename(columns = {coln[0]: 'Chromosome', coln[1]: 'Start', coln[2]: 'End'}, inplace = True)
    else:
        return 'Conversion to PyRanges may be incorrect'

    if name:
        rdc_.rename(columns = {coln[3]: 'Name'}, inplace = True)
    if strand:
        rdc_.rename(columns = {coln[4]: 'Strand'}, inplace = True)
        rdc_ = rdc_[rdc_['Strand'].isin(["+", "-", "."])]
    if save_old_names:
        return pr.PyRanges(rdc_), list(coln)
    else:
        return pr.PyRanges(rdc_)



def save_file_custom(df_: pd.DataFrame, outdir_: str, prev_path: str, suff: str, hdr = True, idx = False, save = True) -> pr.PyRanges:
    """
    Create custom file extension

    """
    pathlib.Path(outdir_).mkdir(parents=True, exist_ok=True)
    save_name = os.path.join(outdir_, Path(prev_path).stem.rsplit('.', maxsplit=1)[0] + suff)
    if save:
        df_.to_csv(save_name, sep = '\t', header=hdr, index=idx)
        return
    else:
        return save_name


