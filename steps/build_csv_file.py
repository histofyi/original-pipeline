from typing import Dict, List, Tuple

import csv
from functions.files import load_facet
from functions.cli import load_config
from functions.pdb import load_pdb_lists


from rich.progress import Progress
from rich.console import Console

config = load_config()
console = Console()


def stringify(array:List) -> str:
    string = ','.join(array)
    return string


def process_alpha_chain_details(row:Dict, facet:Dict) -> Dict:
    row['locus'] = facet['locus']
    row['allele_slug'] = facet['slug']
    row['allele_match'] = facet['match_type']
    return row


def process_species(row:Dict, facet:Dict) -> Dict:
    row['species_slug'] = facet['slug'] 
    return row


def process_chronology(row:Dict, facet:Dict) -> Dict:
    for item in ['deposition_date','release_date','revision_date']:
        row[item] = facet[item]
    row['deposition_year'] = facet['deposition_date'][:4]
    return row


def process_experiment(row:Dict, facet:Dict) -> Dict:
    row['experimental_method'] = facet['experimental_method']
    resolution = facet['resolution']
    if resolution is not None:
        row['resolution'] = f'{resolution:.2f}'
    return row


def process_complex_types(row:Dict, facet:Dict) -> Dict:
    if facet is not None:
        row['complex_type'] = facet['slug']
        row['components'] = stringify(facet['components'])
    return row


def process_peptide_details(row:Dict, facet:Dict) -> Dict:
    i = 0
    if facet is not None:
        for item in facet:
            if i == 0:
                print (facet[item])
                row['peptide_sequence'] = facet[item]['full_sequence']
                row['peptide_length'] = facet[item]['full_length']
                row['peptide_features'] = stringify(facet[item]['features'])
            i += 1
    return row


def process_tcr_details(row:Dict, facet:Dict) -> Dict:
    if facet is not None:
        row['has_tcr'] = True
    return row


def process_publication_details(row:Dict, facet:Dict) -> Dict:
    row['publication_url'] = facet['bibjson']['url']
    if facet['open_access'] == 'y':
        row['open_access'] = True
    return row

def process_assemblies(row:Dict, facet:Dict) -> Dict:
    if facet is not None:
        row['assembly_count'] = len(facet)
    return row


def process_urls(row:Dict, pdb_code:str) -> Dict:
    row['rcsb_url'] = f'https://www.rcsb.org/structure/{pdb_code}'
    row['histo_url'] = f'https://wwww.histo.fyi/structures/view/{pdb_code}'
    row['pdbe_url'] = f'https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_code}'
    return row


facets = {
    'alpha_chain_details':process_alpha_chain_details, 
    'assemblies':process_assemblies,
    'chronology':process_chronology,
    'complex_types':process_complex_types,
    'experiment':process_experiment,
    'peptide_details':process_peptide_details,
    'publication_details':process_publication_details,
    'species':process_species, 
    'tcr_details':process_tcr_details
}






def get_blank_data_dict():
    data_dict = {
        'pdb_code':None,
        'locus':None,
        'allele_slug':None,
        'allele_match': None,
        'species_slug':None,
        'peptide_sequence':None,
        'peptide_length':None,
        'peptide_features':None,
        'resolution':None,
        'deposition_date':None,
        'deposition_year':None,
        'release_date':None,
        'revision_date':None,
        'experimental_method':None,
        'publication_url':None,
        'open_access':False,
        'assembly_count':None,
        'complex_type':None,
        'components':None,
        'has_tcr':False,
        'histo_url':None,
        'rcsb_url':None,
        'pdbe_url':None
    }
    return data_dict




pdb_codes = ['1hhk']

pdb_codes = load_pdb_lists('class_i', config['WAREHOUSE_PATH'], console)

console.rule(f'[bold]Loading facets for structures')


rows = []
csv_rows = []
csv_filepath = f'{config["WAREHOUSE_PATH"]}/datacompilations/core.csv'



for pdb_code in pdb_codes:
    print (pdb_code)
    row = get_blank_data_dict()
    row['pdb_code'] = pdb_code
    for facet in facets:
        current_facet = load_facet(pdb_code, facet)
        row = facets[facet](row, current_facet)
    row = process_urls(row, pdb_code)
    rows.append(row)

i = 0

for row in rows:
    if i == 0:
        headers = [key for key in row]
        csv_rows.append(headers)
    csv_row = [row[key] for key in row]
    csv_rows.append(csv_row)
    i += 1

print (csv_rows)


with open(csv_filepath, 'w', newline='') as csvfile:
     csv_writer = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
     for csv_row in csv_rows:
        csv_writer.writerow(csv_row)
