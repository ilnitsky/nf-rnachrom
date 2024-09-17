import argparse
import pandas as pd
import subprocess
import os 
import warnings
from natsort import natsort_keygen
from pandas.api.types import union_categoricals

from io_.fileops import *
from io_.schemas import *

warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Construct smoothed background profile')
parser.add_argument('rdc', type=str, help='Path to the RNA-DNA contacts file')
parser.add_argument('Stereogene', type=str, help='Path to the Stereogene Smoother')
parser.add_argument('config', type=str, help='Path to the config template file')
parser.add_argument('chrsizes', type=str, help='Path to the file with chromosome sizes')
parser.add_argument('--filter_top_n', type=int, nargs='?', default = 50, help='Drop the most contacting RNAs (rank, default: 50)')
parser.add_argument('--filter_tail_n', type=int, nargs='?', default = 1000, help='Drop the least contacting RNAs (rank, default: 1000)')
parser.add_argument('--outdir', type=str, nargs='?', default = './results', help='')
args = parser.parse_args()


def run_background(args):

    rdc = load_rdc(args.rdc, header = 0, names = voted_BED, use_cols = ['rna_chr', 'dna_chr', 'dna_bgn', 'dna_end', 'gene_name', 'gene_type', 'id'], dtype = rdc_dtypes)
    # print(rdc.head)
    union = union_categoricals([rdc.rna_chr, rdc.dna_chr]).categories
    df = pd.DataFrame()
    df['rna_chr'] = rdc.rna_chr.cat.set_categories(union)
    df['dna_chr'] = rdc.dna_chr.cat.set_categories(union)
    
    all_stat = stats_by_chr(rdc, colname = 'all_contacts')

    back_pc_rnas = rdc.query('gene_type == "protein_coding"')['gene_name'].value_counts()[args.filter_top_n:-args.filter_tail_n].index
    rdc = rdc.query("gene_name in @back_pc_rnas & rna_chr != dna_chr")
    back_stat = stats_by_chr(rdc, colname = 'background_contacts')

    save_file_custom(rdc[['dna_chr', 'dna_bgn', 'dna_end', 'id']], args.outdir, args.rdc, Stereogene_bed, hdr =False)
    save_file_custom(stats_agg(all_stat, back_stat), args.outdir, args.rdc, bg_suffix_stat, idx = True, hdr = True)

    # mod_cfg = 'cfg_mod.cfg'
    mod_cfg = 'smoother.cfg'

    prepeare_cfg(cfg_path = args.config, 
                 chrsizes_path = args.chrsizes,
                 mod_cfg_path = mod_cfg)
    
    cmd = prepare_cmd(path_smoother = args.Stereogene,
                  path_cfg = mod_cfg,
                  path_back = args.rdc,
                  outdir = args.outdir)
    print(cmd)
    run_smoother(cmd)

    return
    

def prepeare_cfg(cfg_path: str, chrsizes_path: str, mod_cfg_path: str):
    """
    Modify Stereogene config file
    
    """
    
    with open(cfg_path, 'r') as raw_cfg:
        content = raw_cfg.read()
        
        for name, value in [('CHRSIZE', chrsizes_path), ('REPLACE', args.outdir)]:
            content = content.replace(name, value)

    with open(mod_cfg_path, 'w') as mod_cfg:
        mod_cfg.write(content)
    
    return 


def prepare_cmd(path_smoother: str, path_cfg: str, path_back: str, outdir: str) -> list:
    
    """
    Prepare bash command to run Smoother.

    Parameters
    ----------
    path_smoother : str
        Path to Smoother
        
    path_cfg: str
        Path to the modified config
        
    path_back: str
        Path to the background contacts file
         
    Returns
    -------
    List with a command
    
    """
    args_for_sg = [path_smoother,
        os.path.abspath(path_cfg),
        save_file_custom('', '', path_back, Stereogene_bed, save = False)]
    
    cmd = "{0} cfg={1} {2}".format(args_for_sg[0], args_for_sg[1], args_for_sg[2])
    print(cmd)
    return [cmd]


def run_smoother(cmd: list):
    
    sb = subprocess.Popen(cmd, shell=True, executable='/bin/bash', text=True)
    
    sb.wait()
    
    if sb.returncode != 0:
        print('Non-zero exit status')
        return sb.communicate()[0]
    
    return 


def stats_by_chr(cnts: pd.DataFrame, colname = 'contacts') -> pd.DataFrame:
    
    return cnts.groupby(['rna_chr']).agg(contacts=('gene_name', 'count')).sort_values(['rna_chr'], ascending = True, key=natsort_keygen()).rename(columns = {'contacts': colname})
    

def stats_agg(all_stat: pd.DataFrame, back_stat: pd.DataFrame) -> pd.DataFrame:
    
    stats = all_stat.join(back_stat)
    stats['background_contacts, %'] = round(back_stat.iloc[:, 0] / all_stat.iloc[:, 0] * 100, 3)
    
    return pd.DataFrame(stats)

run_background(args)
