import pandas as pd
import pathlib
import os

from natsort import natsort_keygen
from pathlib import Path

from io_.fileops import save_file_custom
from io_.schemas import *


def calculate_stats(frames_d: dict, rdc_path: str, outdir_path = './results'):
    
    frames_d['voted'] = frames_d['voted'].sort_values('gene_name', ascending=True)
    
    save_file_custom(stats_by_gene_type(frames_d['voted']), outdir_path, rdc_path, voted_stat_cnts_gene_type, idx = True)
    
    save_file_custom(density_by_gene(frames_d['voted']), outdir_path, rdc_path, voted_stat_dens)
    
    save_file_custom(counts_by_gene(frames_d['voted']), outdir_path, rdc_path, voted_stat_cnts_by_genes, idx = True)
    
    save_file_custom(counts_by_source(frames_d['voted']), outdir_path, rdc_path, voted_stat_cnts_by_source, idx = True)

    save_file_custom(counts_by_chr(frames_d), outdir_path, rdc_path, voted_stat_cnts_by_chr, idx = True)

    return
    


def stats_by_gene_type(vote_sorted):
    
    return vote_sorted.groupby(['Chromosome', 'gene_type']).agg(n_contacts=('gene_name', 'count'), 
                                           n_genes=('gene_name', 'nunique'))
                                                                                                                                                            
def density_by_gene(vote_sorted):
    
    return vote_sorted[['Chromosome', 'gene_name', 'gene_type', 'density']].sort_values(['Chromosome', 'gene_name'], ascending = True, key=natsort_keygen()).drop_duplicates('gene_name')   
    


def counts_by_gene(vote_sorted):
    
    return vote_sorted.groupby(['Chromosome', 'gene_name']).agg(count=('gene_name', 'count'))     


def counts_by_source(vote_sorted):

    return vote_sorted.groupby(['Chromosome', 'Source', 'gene_type']).agg(n_contacts=('gene_name', 'count'),
                                           n_genes=('gene_name', 'nunique')).query('n_contacts > 0')


def counts_by_chr(d):
    all_frames_stat = dict(map(lambda item: (item[0], item[1].groupby('Chromosome')['id'].nunique()), d.items()))
    
    return pd.DataFrame.from_dict(all_frames_stat).rename(columns={'complement_annot': 'Nopposite',
                                                                   'voted': 'Nvoted',
                                                                   'singletons': 'Nsingeltons',
                                                                   'voted_noribo': 'NnoRibo'}).sort_values(by="Chromosome", key=natsort_keygen())
