import argparse
import pandas as pd
import numpy as np
import pyranges as pr

from io_.fileops import *
from io_.schemas import *

from BaRDIC.bardic.api.schemas import BedValidator, ContactsValidator

parser_input = argparse.ArgumentParser(
                    prog='input',
                    description='Validates annotation file')

parser_input_group = parser_input.add_argument_group('Input')
parser_input_group.add_argument('annot', type=str, help='Path to the file with gene annotation')
parser_input_group.add_argument('chrsizes', type=str, help='Path to the file with chromosome sizes')
parser_output_group = parser_input.add_argument_group('Output')
parser_output_group.add_argument('--outdir', type=str, nargs='?', default = './results', help='A folder in which to store the results')
args = parser_input.parse_args()

chrsz = load_chrsizes_as_pyranges(args.chrsizes)
chrsizes_dict = chrsz.as_df().set_index('Chromosome').to_dict()['End']

ann = load_BED(args.annot, names = BED6)


good_annot, bad_annot = BedValidator(chrsizes_dict).validate_df(ann)
save_file_custom(good_annot, args.outdir, args.annot, bardic_input_annot, hdr = False)
save_file_custom(bad_annot, args.outdir, args.annot, bardic_bad_annot)
