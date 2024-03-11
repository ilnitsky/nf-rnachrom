import argparse
import pandas as pd
import numpy as np
import pyranges as pr

from io_.fileops import *
from io_.schemas import *

from BaRDIC.bardic.api.schemas import BedValidator, ContactsValidator

parser_input = argparse.ArgumentParser(
                    prog='input',
                    description='Validates contacts file')

parser_input_group = parser_input.add_argument_group('Input')
parser_input_group.add_argument('rdc_file', type=str, help='Path to the RNA-DNA contact file')
parser_input_group.add_argument('annot', type=str, help='Path to the file with gene annotation')
parser_input_group.add_argument('chrsizes', type=str, help='Path to the file with chromosome sizes')
parser_output_group = parser_input.add_argument_group('Output')
parser_output_group.add_argument('--outdir', type=str, nargs='?', default = './results', help='A folder in which to store the results')

args = parser_input.parse_args()

chrsz = load_chrsizes_as_pyranges(args.chrsizes)
chrsizes_dict = chrsz.as_df().set_index('Chromosome').to_dict()['End']

ann = load_BED(args.annot, names = BED6)
print(ann.head())

rdc = load_BED(args.rdc_file, names = BED6)

g_rdc, b_rdc = ContactsValidator(chrsizes_dict, dict(zip(ann['name'], [None]*len(ann['name'])))).validate_df(rdc)

save_file_custom(g_rdc, args.outdir, args.rdc_file, bardic_input_rdc, hdr = False)
save_file_custom(b_rdc, args.outdir, args.rdc_file, bardic_bad_rdc)
