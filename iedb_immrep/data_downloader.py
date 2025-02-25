import requests
import os
import pandas as pd
import time

tcr_search_params = {'mhc_classes':'cs.{I}', 'qualitative_measures':'cs.{Positive}',
                     'host_organism_iris':'cs.{NCBITaxon:9606}', 'order':'receptor_group_id',
                     'receptor_type':'eq.alphabeta', 'offset':0}

base_uri = 'https://query-api.iedb.org'
def rename_func(str_name):
    prefix = str_name.split('__')[0].capitalize().replace('_', ' ')
    suffix = str_name.split('__')[1].replace('_', ' ')
    suffix = ' '.join([p.capitalize() for p in suffix.split(' ')]).replace('Cdr', "CDR").replace('Iedb', 'IEDB').replace('Ids', 'IDs').replace('Id','ID')
    suffix = suffix.replace('Mhc', 'MHC').replace('Iri', 'IRI')
    return prefix + ' - ' + suffix
def api_to_export(api_table):
    api_table.columns = [rename_func(p) if '__' in p else p for p in api_table.columns]
    return api_table
def read_tcr_search_table(tcr_search_params):
    results = []
    result = requests.get(os.path.join(base_uri, 'tcr_search'), params=tcr_search_params)
    results.append(pd.json_normalize(result.json()))
    while result.json() != []:
        tcr_search_params['offset'] += 10000
        result = requests.get(os.path.join(base_uri, 'tcr_search'), params=tcr_search_params)
        results.append(pd.json_normalize(result.json()))
    results = pd.concat(results)
    return results

def get_tcr_export():
    tcr_export_params = {'offset': 0, 'order': 'receptor_group_id'}
    results_export = []
    result = requests.get(os.path.join(base_uri, 'tcr_export'), params=tcr_export_params)
    results_export.append(pd.json_normalize(result.json()))
    while result.json() != []:
        tcr_export_params['offset'] += 10000
        result = requests.get(os.path.join(base_uri, 'tcr_export'), params=tcr_export_params)
        results_export.append(pd.json_normalize(result.json()))
    results_export = pd.concat(results_export)
    return results_export

def read_epitope_table(structure_ids):
    epitope_table = []
    for p in range(0, len(structure_ids), 100):
        assay_id_subset = ','.join(structure_ids[p: p + 100])
        assay_search = {'structure_iri': f'in.({assay_id_subset})'}
        epitope = pd.json_normalize(requests.get(os.path.join(base_uri, 'epitope_search'), params=assay_search).json())
        epitope_table.append(epitope)
    epitope_table = pd.concat(epitope_table)
    return epitope_table

def filter_export(results_export, receptor_ids):
    results_export = results_export[results_export['receptor__iedb_receptor_id'].isin(receptor_ids)]
    results_export = results_export[results_export['assay__mhc_allele_names'].fillna('').str.contains(':')]
    return results_export

def read_assay_table(assay_ids):
    assay_table = []
    for p in range(0, len(assay_ids), 100):
        assay_id_subset = ','.join(assay_ids[p: p + 100])
        assay_search = {'tcell_id': f'in.({assay_id_subset})'}
        epitope = pd.json_normalize(requests.get(os.path.join(base_uri, 'tcell_search'), params=assay_search).json())
        assay_table.append(epitope)
    assay_table = pd.concat(assay_table)
    return assay_table


print(f'Downloading data from IEDB...')
results = read_tcr_search_table(tcr_search_params)
export = get_tcr_export()
receptor_ids = set().union(*[set(p) for p in results['receptor_ids'].values])
print(f'Retrieved {len(receptor_ids)} receptor ids...')
output = filter_export(export, receptor_ids)
print(f'Downloading assay data...')
assay_ids = set(output['assay__iedb_ids'].unique())
assay_ids = list(set().union(*[set(p.split(', ')) for p in assay_ids]))
structure_ids = list(set().union(*[set(p.split(', ')) for p in output['epitope__iedb_iri'].unique()]))
structure_ids = [f'IEDB_EPITOPE:{p.split("/")[-1]}' for p in structure_ids]
assay_table = read_assay_table(assay_ids)
print(f'Downloading epitope data...')
epitope_table = read_epitope_table(structure_ids)
results_export = api_to_export(output)
results_export['Epitope - Structure IRI'] = results_export['Epitope - IEDB IRI'].map(lambda x:'IEDB_EPITOPE:'+x.split('/')[-1])
results_export = pd.merge(left = results_export, right=epitope_table[['e_modification', 'structure_iri', 'structure_type']], left_on='Epitope - Structure IRI', right_on='structure_iri')
results_export = results_export[(results_export['structure_type']=='Linear peptide')&
                                (results_export['e_modification'].map(str)=='None')]

assay_infos = dict(assay_table[['tcell_id','assay_names']].values)
results_export['Assay - Names'] = results_export['Assay - IEDB IDs'].map(lambda x:','.join(set([assay_infos[int(p)] for p in x.split(', ')])))
results_export['Epitope ID'] = results_export['Epitope - Name'] + '~' + results_export['Assay - MHC Allele Names'].replace(', ', '+')
results_export['Epitope IDs'] = results_export.apply(lambda x:','.join([x['Epitope - Name']+ '~' + p for p in x['Assay - MHC Allele Names'].split(', ')]),axis=1)
results_export.to_parquet(f'dat/api_export_results.parquet')
