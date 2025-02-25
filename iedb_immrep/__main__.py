from iedb_immrep.total_process import IEDBDataset
from iedb_immrep.utils import *
from argparse import ArgumentParser
import os
import sys
import warnings
warnings.filterwarnings("ignore")

parser = ArgumentParser()
parser.add_argument('--export_file', dest='export_file', help='Parquet export file from IEDB.')
parser.add_argument('--output_dir', dest='output_dir', help='directory to put output')
parser.add_argument('--chains', dest='chains', default='paired', choices=['paired', 'beta'])
args = parser.parse_args()

if os.path.exists(args.output_dir) is False:
    os.mkdir(args.output_dir)
if os.path.exists(args.export_file) is False:
    sys.exit()

chains = ['Chain 1', 'Chain 2'] if args.chains == 'paired' else ['Chain 2']

print(f'Reading export file...')
iedb = IEDBDataset(args.export_file, file_type='parquet', chains=chains)
print(f'Filtering on MHC allele names...')
iedb.filter_mhc()
print(f'Prioritizing Calculated over Curated genes...')
iedb.prioritize_entries(prioritize_func=prioritize_vj)
print(f'Gene name standardization...')
iedb.standardize_genes(correct_j=True)
print(f'CDR assignment...')
iedb.assign_cdrs()
iedb.format()
print(f'Running Thimble to get full length sequences...')
iedb.run_stitchr()
iedb.format_for_competition()
print(f'Writing.')
iedb.format_for_competition(output_path=os.path.join(args.output_dir, 'immrep_IEDB.csv'))
iedb.format_airr(output_path=os.path.join(args.output_dir, 'immrep_IEDB_airr.tsv'))