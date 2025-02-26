from itertools import product
import pandas as pd
from Bio import SeqIO
import os
import tidytcells as tt
import subprocess
import requests
import numpy as np
import iedb_immrep
import secrets
from Stitchr import stitchrfunctions as fxn
from Bio.Seq import translate




default_allele = '*01'
def get_leader_constant_seqs():
    leaders = {}
    constant = {}
    for chain in ['TRA', 'TRB']:
        imgt_dat, tcr_functionality, partial = fxn.get_imgt_data(chain, ['LEADER', 'CONSTANT', 'VARIABLE', 'JOINING'],
                                                                 'HUMAN')
        leaders[chain] = set()
        constant[chain] = set()
        for gene in imgt_dat['LEADER']:
            for allele in imgt_dat['LEADER'][gene]:
                nt = imgt_dat['LEADER'][gene][allele]
                offset = len(nt) % 3
                leaders[chain].add(translate(nt[offset:]).split('*')[0])
        for gene in imgt_dat['CONSTANT']:
            for allele in imgt_dat['CONSTANT'][gene]:
                nt = imgt_dat['CONSTANT'][gene][allele]
                for offset in range(0, 3):
                    nt_seq = nt[offset:]
                    ns = len(nt_seq) % 3
                    seq = translate(nt_seq+'N'*ns).split('*')[0]
                    if len(seq) > 50:
                        constant[chain].add(seq)
    return leaders, constant


leaders, constants = get_leader_constant_seqs()
data_path = '/'.join(iedb_immrep.__file__.split('/')[:-1])

def get_cdrs(path: str = os.path.join(data_path, 'dat/trv_cdrs.tsv')):
    output = pd.read_csv(path, sep='\t')
    cdr1s = dict(output[['sequence_id', 'cdr1_aa']].values)
    cdr2s = dict(output[['sequence_id', 'cdr2_aa']].values)
    return cdr1s, cdr2s


def studies_to_flag(xlsx_path: str = 'immrep_exclusion.xlsx', method='IEDB'):
    exclude = pd.read_excel(xlsx_path)
    if method == 'IEDB':
        return set(exclude[(exclude['Exclude'] == True)&(exclude['Reference - IEDB IRI'].notna())]["Reference - IEDB IRI"])
    else:
        return set(exclude[(exclude['Exclude'] == True)&(exclude['PMID'].notna())]['PMID']).union(set(exclude[(exclude['Exclude'] == True)&(exclude['VDJdb_alias'].notna())]['VDJdb_alias']))

def summarize_counts(df, receptor_key: str = 'Receptor - IEDB Receptor ID', epitope_key: str = 'Epitope - Name'):
    return {'total': df.shape[0], 'receptors': df[receptor_key].nunique(), 'epitopes': df[epitope_key].nunique()}

def prioritize_vj(df, chain):
    v_call = df[f'{chain} - Calculated V Gene'].where(df[f'{chain} - Calculated V Gene'].notna(),
                                                      df[f'{chain} - Curated V Gene'])
    j_call = df[f'{chain} - Calculated J Gene'].where(df[f'{chain} - Calculated J Gene'].notna(),
                                                      df[f'{chain} - Curated J Gene'])
    cdr3 = df[f'{chain} - CDR3 Calculated'].where(df[f'{chain} - CDR3 Calculated'].notna(),
                                                  df[f'{chain} - CDR3 Curated'])
    return pd.DataFrame({'v_call': v_call, 'junction': cdr3, 'j_call': j_call})


def assert_correct_chain_order(df):
    if 'TRB' in str(df['Chain 1 - v_call']):
        df['Chain 1 - v_call'], df['Chain 2 - v_call'] = df['Chain 2 - v_call'], df['Chain 1 - v_call']
        df['Chain 1 - j_call'], df['Chain 2 - j_call'] = df['Chain 2 - j_call'], df['Chain 1 - j_call']
        df['Chain 1 - junction'], df['Chain 2 - junction'] = df['Chain 2 - junction'], df['Chain 1 - junction']
    if 'TRA' in str(df['Chain 2 - v_call']):
        df['Chain 1 - v_call'], df['Chain 2 - v_call'] = df['Chain 2 - v_call'], df['Chain 1 - v_call']
        df['Chain 1 - junction'], df['Chain 2 - junction'] = df['Chain 2 - junction'], df['Chain 1 - junction']
        df['Chain 1 - j_call'], df['Chain 2 - j_call'] = df['Chain 2 - j_call'], df['Chain 1 - j_call']
    return df

def run_thimble(df, chains=['A', 'B'], missing_strategy='remove'):
    outputs = {}
    for chain in chains:
        format_for_thimble(df, chain=chain, output=f'thimble{chain}.tsv', missing_strategy=missing_strategy)
        full_length_trs, new_cdr3s = run_thimble_single_chain(f'thimble{chain}.tsv', chain=chain)
        full_length_trs = {p: remove_leader_constant(full_length_trs[p], 'TR' + chain) for p in full_length_trs}
        os.remove(f'thimble{chain}.tsv',)
        outputs[chain] = full_length_trs, new_cdr3s
    return outputs


def run_thimble_single_chain(thimble_input_path, chain='A', retry=True):
    print(f'Running Thimble for chain {chain} with terminal Fs.')
    output_path = f'output_thimble.tsv'
    out_subprocess = subprocess.call(f'thimble -in {thimble_input_path} -o {output_path} -s human', shell=True)
    if out_subprocess != 0:
        raise FileNotFoundError("Thimble didn't run - did you run stitchrdl?")
    outputs = pd.read_csv(output_path, sep='\t')
    fixed = outputs[outputs[f'TR{chain}_aa'].notna()]
    output = dict(fixed[['TCR_name', f'TR{chain}_aa']].values)
    os.remove(output_path)
    if retry:
        perc_failed = (1 - fixed.shape[0] / outputs.shape[0]) * 100
        print(f'Running Thimble for chain {chain} on failures ({perc_failed}%) with terminal Ws.')
        retry_path = thimble_input_path.replace(".tsv", "_redo.tsv")
        input_file = pd.read_csv(thimble_input_path, sep='\t')
        to_redo = set(input_file['TCR_name'].unique()).difference(set(fixed['TCR_name'].unique()))
        to_redo_df = input_file[input_file['TCR_name'].isin(to_redo)].copy()
        to_redo_df[f'TR{chain}_CDR3'] = to_redo_df[f'TR{chain}_CDR3'].map(lambda x: x[:-1] + 'W')
        to_redo_df.to_csv(retry_path, sep='\t', index=False)
        subprocess.call(f'thimble -in {retry_path} -o {output_path} -s human', shell=True)
        outputs = pd.read_csv(output_path, sep='\t')
        output.update(dict(outputs[outputs[f'TR{chain}_aa'].notna()][['TCR_name', f'TR{chain}_aa']].values))
        new_cdr3bs = dict(outputs[outputs[f'TR{chain}_aa'].notna()][['TCR_name', f'TR{chain}_CDR3']].values)
        os.remove(output_path)
        os.remove(retry_path)
    else:
        new_cdr3bs = {}
        os.remove(thimble_input_path)
    return output, new_cdr3bs


def format_for_thimble(df_input, chain='A', output='thimble.tsv', missing_strategy='remove'):
    df = df_input.copy()
    thimble_columns = ['TCR_name', 'TRAV', 'TRAJ', 'TRA_CDR3', 'TRBV', 'TRBJ', 'TRB_CDR3', 'TRAC', 'TRBC', 'TRA_leader',
                       'TRB_leader',
                       'Linker', 'Link_order', 'TRA_5_prime_seq', 'TRA_3_prime_seq', 'TRB_5_prime_seq',
                       'TRB_3_prime_seq']
    df['TCR_name'] = range(0, df.shape[0])
    df[f'TR{chain}_CDR3'] = 'C' + df[f'CDR{chain}3'] + 'F'
    for column in thimble_columns:
        if column not in df.columns:
            df[column] = np.nan
    for column in ['TRAV', 'TRAJ', 'TRA_CDR3', 'TRBV', 'TRBJ', 'TRB_CDR3']:
        if column[2] != chain:
            df[column] = np.nan
        elif missing_strategy != 'remove':
            df[column] = df[column].map(lambda x: x.replace('unknown', '%'))
        else:
            df = df[df[column] != 'unknown']
    df[thimble_columns].to_csv(output, sep='\t', index=False)

def filter_for_information(df, chains=['Chain 1', 'Chain 2'], fields=['v_call', 'j_call', 'junction']):
    filter_chains = [f'{chain} - {field}' for chain, field in product(chains, fields)]
    return df[df[filter_chains].notna().all(axis=1)]

def fix_junction(cdr3):
    return tt.junction.standardize(cdr3)

def fix_junctions(df, chains=['Chain 1', 'Chain 2']):
    for chain in chains:
        df[f'{chain} - junction'] = df[f'{chain} - junction'].map(fix_junction)
    return df

def write_fasta(df, fname, seq='nucleotides', chains=['Chain 1', 'Chain 2']):
    if seq == 'nucleotides':
        name = 'Nucleotide Sequence'
    else:
        name = 'Protein Sequence'
    seqs = set().union(*[set(df[(df[f'{chain} - {name}'].notna())][f'{chain} - {name}'].unique()) for chain in chains])
    with open(fname, 'w') as k:
        for i, seq in enumerate(seqs):
            k.write(f'>seq_{i}\n{seq}\n')
    return fname

def valid_mhcclassI(mro_path=None):
    if mro_path is None:
        mro_path = f'dat/MRO_molecules.tsv'
        with open(mro_path, 'w') as k:
            k.writelines(requests.get(
                f'https://raw.githubusercontent.com/IEDB/MRO/refs/heads/master/ontology/molecule.tsv').content)
    mro = pd.read_csv(mro_path, sep='\t')
    mro = mro[mro['In Taxon'] == 'human']
    allowed_alleles = set(mro[(mro['Parent'] == 'MHC class I protein complex')]['IEDB Label'].unique())
    return allowed_alleles

def remove_leader_constant(sequence, chain):
    leader_removed = False
    constant_removed = False
    for leader_seq in leaders[chain]:
        if sequence[:len(leader_seq)] == leader_seq:
            leader_removed = True
            sequence = sequence[len(leader_seq):]
    for constant_seq in constants[chain]:
        if sequence[-len(constant_seq):] == constant_seq:
            constant_removed = True
            sequence = sequence[:-len(constant_seq)]
    if (constant_removed == True) & (leader_removed == True):
        return sequence
    else:
        raise Exception("Unable to remove constant region")


def generate_hash():
    return secrets.token_hex(4)

def validate_file(df, cdr_keys=['CDRA1', 'CDRA2', 'CDRA3'], seq_key='TRA'):
    for key in cdr_keys:
        df[f'{key}_val'] = df.apply(lambda x: x[key] in x[seq_key] if pd.notna(x[seq_key]) else True, axis=1)
    return df

def cdrs_from_ndm_file(ndm_path='dat/igblast/internal_data/human/human.ndm.imgt',
                       seq_path='dat/igblast/fasta/imgt_human_tr_v.fasta', outf='dat/trv_cdrs.tsv'):
    input = pd.DataFrame([p.strip('\n').split() for p in open(ndm_path, 'r').readlines() if p.startswith('TR')])
    input.columns = ['allele', 'fwr1_start', 'fwr1_end', 'cdr1_start', 'cdr1_end', 'fwr2_start', 'fwr2_end',
                     'cdr2_start', 'cdr2_end', 'fwr3_start', 'fwr3_end', 'locus', 'coding_start']
    for p in input.columns[1:]:
        if p != 'locus':
            input[p] = input[p].map(int)
    for p in input.columns[1:-2]:
        input[p] = input[p] + input['coding_start']
    seqs = dict((p.id, str(p.seq)) for p in SeqIO.parse(seq_path, format='fasta'))
    cdrs = {}
    for k, p in input.iterrows():
        seq = seqs[p.allele]
        cdr1 = seq[p.cdr1_start - 1:p.cdr1_end]
        cdr2 = seq[p.cdr2_start - 1:p.cdr2_end]
        cdr1_aa = translate(cdr1)
        cdr2_aa = translate(cdr2)
        cdrs[p.allele] = {'cdr1_aa': cdr1_aa, 'cdr2_aa': cdr2_aa}
    pd.DataFrame(cdrs).T.reset_index().rename(columns={'index': 'sequence_id'}).to_csv(outf, sep='\t')


