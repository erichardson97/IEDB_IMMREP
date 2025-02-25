import airr as airrc
from .utils import *

class IEDBDataset:
    def __init__(self, file, file_type='xlsx', chains=['Chain 1', 'Chain 2']):
        if file_type == 'xlsx':
            self.df = pd.read_excel(file)
        else:
            self.df = pd.read_parquet(file)
        self.counts = {'Total': summarize_counts(self.df)}
        self.ids = {'Total': set(self.df['Receptor - IEDB Receptor ID'])}
        cdr1s, cdr2s = get_cdrs()
        self.cdr1s = cdr1s
        self.cdr2s = cdr2s
        self.chains = chains
        self.has_chain = None
        self.final_format = None
        self.airr = None

    def filter_mhc(self):
        self.df = self.df[(self.df['Assay - MHC Allele Names'].fillna('').str.contains(':')) & (
            ~self.df['Assay - MHC Allele Names'].str.contains('mutant'))]
        self.counts['MHC'] = summarize_counts(self.df)
        self.ids['MHC'] = set(self.df['Receptor - IEDB Receptor ID'])

    def prioritize_entries(self, prioritize_func=prioritize_vj):
        df = self.df
        df_chain1 = prioritize_func(df, 'Chain 1')
        df_chain2 = prioritize_func(df, 'Chain 2')
        df_chain1.columns = ['Chain 1 - ' + p for p in df_chain1.columns]
        df_chain2.columns = ['Chain 2 - ' + p for p in df_chain2.columns]
        paired = pd.merge(left=df, right=df_chain1, left_index=True, right_index=True)
        paired = pd.merge(left=paired, right=df_chain2, left_index=True, right_index=True).fillna(value=np.nan)
        has_chain = filter_for_information(paired, chains=self.chains, fields=['v_call', 'junction'])
        self.has_chain = has_chain
        self.ids['VJCDR3 for Relevant chains'] = set(self.has_chain['Receptor - IEDB Receptor ID'])
        self.counts['VJCDR3 for Relevant chains'] = summarize_counts(self.has_chain)

    def standardize_genes(self, correct_j=True):
        has_chain = self.has_chain
        for chain in self.chains:
            has_chain[f'{chain} - v_call_tt'] = has_chain[f'{chain} - v_call'].map(tt.tr.standardize)
            has_chain = has_chain[(has_chain[f'{chain} - v_call_tt'].notna())]
            if correct_j:
                has_chain[f'{chain} - j_call_tt'] = has_chain[f'{chain} - j_call'].map(
                    lambda x: tt.tr.standardize(x) if pd.notna(x) else 'unknown')
                has_chain = has_chain[(has_chain[f'{chain} - j_call_tt'].notna())]
        for chain in self.chains:
            has_chain[f'{chain} - v_call_tt'] = has_chain[f'{chain} - v_call_tt'].map(
                lambda x: x + default_allele if '*' not in x else x)
            if correct_j:
                has_chain[f'{chain} - j_call_tt'] = has_chain[f'{chain} - j_call_tt'].map(
                    lambda x: x + default_allele if ('*' not in x) & (x != 'unknown') else x)
        has_chain = assert_correct_chain_order(has_chain)
        has_chain = fix_junctions(has_chain, chains=self.chains)
        for chain in self.chains:
            has_chain = has_chain[(has_chain[f'{chain} - junction'].notna())]
        self.has_chain = has_chain
        self.ids['Standardized'] = set(self.has_chain['Receptor - IEDB Receptor ID'].unique())
        self.counts['VJCDR3 Standardized'] = summarize_counts(self.has_chain)

    def run_stitchr(self):
        print(f'Running Thimble...')
        chain_mapping = {'Chain 1': 'A', 'Chain 2': 'B'}
        thimble_outputs = run_thimble(self.all_cdrs, [chain_mapping[x] for x in self.chains])
        self.all_cdrs['TCR_name'] = range(self.all_cdrs.shape[0])
        if 'Chain 1' in self.chains:
            self.all_cdrs['TRA'] = self.all_cdrs['TCR_name'].map(thimble_outputs['A'][0])
        if 'Chain 2' in self.chains:
            self.all_cdrs['TRB'] = self.all_cdrs['TCR_name'].map(thimble_outputs['B'][0])

    def assign_cdrs(self):
        has_chain = self.has_chain
        mappings = {'Chain 1': 'A', 'Chain 2': 'B'}
        keys = []
        for chain in self.chains:
            has_chain[f'CDR{mappings[chain]}1'] = has_chain[f'{chain} - v_call_tt'].map(self.cdr1s)
            has_chain[f'CDR{mappings[chain]}2'] = has_chain[f'{chain} - v_call_tt'].map(self.cdr2s)
            has_chain[f'CDR{mappings[chain]}3'] = has_chain[f'{chain} - junction'].map(lambda x: x[1:-1])
            keys += [f'CDR{mappings[chain]}{i}' for i in range(1, 4)]
        all_cdrs = has_chain[has_chain[keys].notna().all(axis=1)]
        self.counts['Assigned CDRs'] = summarize_counts(all_cdrs)
        self.ids['Assigned CDRs'] = set(all_cdrs['Receptor - IEDB Receptor ID'].unique())
        self.all_cdrs = all_cdrs

    def format(self, exclusion_path=os.path.join(data_path, 'dat/immrep_exclusion.xlsx')):
        rename_columns = {'Chain 1 - v_call_tt': 'TRAV', 'Chain 2 - v_call_tt': 'TRBV',
                          'Chain 1 - j_call_tt': 'TRAJ', 'Chain 2 - j_call_tt': 'TRBJ'}
        self.all_cdrs = self.all_cdrs.rename(columns=rename_columns)
        cdr_keys_1 = ['CDRA1', 'CDRA2', 'CDRA3'] if 'Chain 1' in self.chains else []
        cdr_keys_2 = ['CDRB1', 'CDRB2', 'CDRB3'] if 'Chain 2' in self.chains else []
        cdr_keys = cdr_keys_1 + cdr_keys_2
        self.all_cdrs['unique_id'] = self.all_cdrs.apply(lambda x: '_'.join([x[p] for p in cdr_keys + ['Epitope IDs']]),
                                                         axis=1)
        self.all_cdrs = self.flag_10X_studies(self.all_cdrs, exclusion_path=exclusion_path)

    def flag_10X_studies(self, df, exclusion_path=os.path.join(data_path, 'dat/immrep_exclusion.xlsx')):
        studies = studies_to_flag(exclusion_path)
        df['10X'] = df['Reference - IEDB IRI'].map(lambda x:x.replace('http','https')).isin(studies)
        return df

    def write_parquet(self, output):
        self.all_cdrs.to_parquet(output)
        pd.DataFrame(self.counts).to_csv(output.replace('.parquet', '_filtering.csv'))

    def validate(self, df):
        val_columns = []
        if 'Chain 1' in self.chains:
            validated_file = validate_file(df)
            val_columns += ['CDRA1_val', 'CDRA2_val', 'CDRA3_val']
        else:
            validated_file = df.copy()
        if 'Chain 2' in self.chains:
            validated_file = validate_file(validated_file, cdr_keys=['CDRB1', 'CDRB2', 'CDRB3'], seq_key='TRB')
            val_columns += ['CDRB1_val', 'CDRB2_val', 'CDRB3_val']
        return validated_file[validated_file[val_columns].all(axis=1)].drop(val_columns, axis=1)

    def format_for_competition(self, output_path='immrep_IEDB.csv'):
        dedup = {}
        mapping = {"Chain 1": 'A', 'Chain 2': 'B'}
        fields = {'A': {'TRAV', 'TRAJ', 'TRA', 'CDRA1', 'CDRA2', 'CDRA3'}, 'B': {'TRBV', 'TRBJ', 'TRB', 'CDRB1',
                                                                                 'CDRB2', 'CDRB3'}}
        sequence_fields = list(set().union(*[fields[mapping[chain]] for chain in self.chains]))
        self.all_cdrs['full_seq_unique'] = self.all_cdrs.apply(
            lambda x: ':'.join([x[f'TR{mapping[chain]}'] if pd.notna(x[f'TR{mapping[chain]}']) else 'NA'
                                for chain in self.chains]), axis=1)
        for k, df in self.all_cdrs.groupby('unique_id'):
            for i, (x, p) in enumerate(df.groupby(f'full_seq_unique')):
                example = p.iloc[0]
                receptor_ids = ','.join([str(p) for p in set(p['Receptor - IEDB Receptor ID'])])
                references = ','.join(set(p['Reference - IEDB IRI']))
                just_10x = set(p['10X'].unique()).difference({True}) == set()
                for phla in example['Epitope IDs'].split(','):
                    unique_key = generate_hash()
                    dedup[unique_key] = {'receptor_id': receptor_ids, 'references': references, 'just_10X': just_10x,
                                         'epitope': phla}
                    dedup[unique_key].update(
                        dict(example[sequence_fields]))
        dedup = pd.DataFrame(dedup).T
        dedup = self.validate(dedup)
        dedup['mhc'] = dedup['epitope'].map(lambda x: x.split('~')[1])
        dedup['epitope'] = dedup['epitope'].map(lambda x: x.split('~')[0])
        rename_columns = {'epitope': 'Peptide', 'mhc': 'HLA', 'TRAV': 'Va', 'TRBV': 'Vb', 'TRAJ': 'Ja', 'TRBJ': 'Jb',
                          'CDRA1': 'CDR1a', 'CDRA2': 'CDR2a', 'CDRA3': 'CDR3a', 'CDRB1': 'CDR1b', 'CDRB2': 'CDR2b',
                          'CDRB3': 'CDR3b', 'TRA': 'TCRa', 'TRB': 'TCRb'}
        self.final_format = dedup
        self.final_format['CDR3b_extended'] = self.final_format.apply(lambda x: 'C' + x.CDRB3 + 'F' if pd.isna(x['TRB']) else x['TRB'][
                                                                                                 x['TRB'].index(
                                                                                                     x['CDRB3']) - 1:x[
                                                                                                                         'TRB'].index(
                                                                                                     x['CDRB3']) + len(
                                                                                                     x['CDRB3']) + 1],
                                        axis=1)
        self.final_format['CDR3a_extended'] = self.final_format.apply(lambda x: 'C' + x.CDRA3 + 'F' if pd.isna(x['TRA']) else x['TRA'][
                                                                                                 x['TRA'].index(
                                                                                                     x['CDRA3']) - 1:x[
                                                                                                                         'TRA'].index(
                                                                                                     x['CDRA3']) + len(
                                                                                                     x['CDRA3']) + 1],
                                        axis=1)
        self.final_format[['epitope','mhc','TRAV','TRAJ','CDRA1','CDRA2','CDRA3','CDR3a_extended','TRA','TRBV','TRBJ','CDRB1','CDRB2','CDRB3','CDR3b_extended','TRB','references','receptor_id','just_10X']].rename(columns=rename_columns).to_csv(output_path, index=False)

    def format_airr(self, output_path='immrep_IEDB_airr.tsv'):
        df = self.final_format.copy().reset_index()
        df = df.rename(columns={'index': 'cell_id'})
        chain1_cols = ['CDRA1', 'CDRA2', 'CDRA3', 'TRAV', 'TRAJ', 'TRA']
        chain2_cols = ['CDRB1', 'CDRB2', 'CDRB3', 'TRBV', 'TRBJ', 'TRB']
        renaming = {p: p[:3].lower() + p[-1] + '_aa' for p in chain1_cols if p != 'TRAV'}
        renaming.update({p: p[:3].lower() + p[-1] + '_aa' for p in chain2_cols if p != 'TRBV'})
        renaming['TRAV'] = 'v_call'
        renaming['TRBV'] = 'v_call'
        renaming['TRAJ'] = 'j_call'
        renaming['TRBJ'] = 'j_call'
        renaming['TRA'] = 'sequence_aa'
        renaming['TRB'] = 'sequence_aa'
        airr = []
        for k, p in df.groupby('cell_id'):
            if 'Chain 1' in self.chains:
                alpha_chains = p[chain1_cols].rename(columns=renaming)
                alpha_chains['locus'] = 'TRA'
                alpha_chains['epitope'] = p.iloc[0]['epitope']
                alpha_chains['mhc'] = p.iloc[0]['mhc']
                airr.append(alpha_chains.assign(cell_id=k))
            if 'Chain 2' in self.chains:
                beta_chains = p[chain2_cols].rename(columns=renaming)
                beta_chains['locus'] = 'TRB'
                beta_chains['epitope'] = p.iloc[0]['epitope']
                beta_chains['mhc'] = p.iloc[0]['mhc']
                airr.append(beta_chains.assign(cell_id=k))
        airr = pd.concat(airr)
        airr['productive'] = True
        self.airr = airr
        airrc.dump_rearrangement(self.airr, output_path)

