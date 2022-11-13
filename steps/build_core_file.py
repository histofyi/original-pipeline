from typing import Dict, List, Tuple

from functions.files import load_facet, load_constants, write_json, write_facet
from functions.cli import load_config
from functions.pdb import load_pdb_lists
from functions.collation import build_core_data
from functions.titles import pick_best_title, clean_title, build_titles

from rich.progress import Progress
from rich.console import Console

import json


config = load_config()
console = Console()

peptide_lengths = load_constants('peptide_lengths')

def process_facet(core:Dict, facet:Dict, facet_type:str, facet_function) -> Dict:
    if facet is not None:
        core = facet_function(core, facet)
    else:
        if facet_type not in ['peptide_details','tcr_details']:
            errors[facet_type].append(core["pdb_code"])
    return core

def name_peptide_length(length:int) -> str:
    peptide_length_name = None
    for peptide_length in peptide_lengths:
        if peptide_lengths[peptide_length]['length'] == int(length):
            peptide_length_name = peptide_length
    return peptide_length_name


def process_alpha_chain_details(core:Dict, facet:Dict) -> Dict:
    core['allele']['alpha'] = facet    
    core['locus'] = facet['locus']
    return core


def process_species(core:Dict, facet:Dict) -> Dict:
    core['species'] = facet
    return core


def process_assigned_chains(core:Dict, facet:Dict) -> Dict:
    core['assigned_chains'] = facet
    return core


def process_chronology(core:Dict, facet:Dict) -> Dict:
    items = ['deposition','release','revision']
    for item in items:
        core['chronology'][f'{item}_date'] = facet[f'{item}_date']
        core['chronology'][f'{item}_year'] = int(facet[f'{item}_date'][:4])
    return core


def process_experiment(core:Dict, facet:Dict) -> Dict:
    core['experimental_method'] = facet['experimental_method']
    resolution = facet['resolution']
    if resolution is not None:
        core['resolution'] = f'{resolution:.2f}'
    return core


def process_complex_types(core:Dict, facet:Dict) -> Dict:
    complex_type_slug = facet['slug']
    if 'class_i' in complex_type_slug:
        core['classical'] = True
        core['class'] = 'class_i'
    else:
        core['classical'] = False
        core['class'] = 'class_i'
    core['complex_type'] = complex_type_slug
    core['unique_chain_count'] = facet['unique_chains']
    core['components'] = facet['components']
    return core


def process_peptide_details(core:Dict, facet:Dict) -> Dict:
    core['peptide'] = facet
    i = 0
    for chain in facet:
        if i == 0:
            core['peptide_sequence'] = facet[chain]['full_sequence']
            core['peptide_length'] = facet[chain]['full_length']
            core['peptide_length_name'] = name_peptide_length(facet[chain]['full_length'])
        i+=1
    return core


def process_tcr_details(core:Dict, facet:Dict) -> Dict:
    if facet is not None:
        core['tcr'] = facet
    return core


def process_publication_details(core:Dict, facet:Dict) -> Dict:
    core['publication'] = facet
    return core


def process_assemblies(core:Dict, facet:Dict) -> Dict:
    core['assemblies'] = facet
    return core


def process_title(core:Dict, facet:Dict) -> Dict:
    core['pdb_title'] = pick_best_title(core, facet)
    core['pdb_title'] = clean_title(core, facet)
    core['title'], core['page_title'] = build_titles(core)
    return core


def process_urls(core:Dict, pdb_code:str) -> Dict:
    core['same_as']['rcsb'] = {'url':f'https://www.rcsb.org/structure/{pdb_code}'}
    core['same_as']['pdbe'] = {'url':f'https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_code}'}
    if core['tcr'] is not None:
        core['same_as']['stcrdab'] = {'url':f'http://opig.stats.ox.ac.uk/webapps/stcrdab/StrViewer?pdb={pdb_code}'}
    return core


def make_listing(core:Dict) -> Dict:
    listing = {}
    listing['pdb_code'] = core['pdb_code']
    listing['pdb_title'] = core['pdb_title']
    listing['components'] = core['components']
    listing['species'] = core['species']['common_name']
    listing['title'] = core['title']
    listing['deposition_year'] = core['chronology']['deposition_year']
    if 'peptide_length' in core:
        listing['peptide_length'] = core['peptide_length']
        listing['peptide_length_name'] = core['peptide_length_name']
    return listing


facets = {
    'alpha_chain_details':process_alpha_chain_details, 
    'assemblies':process_assemblies,
    'assigned_chains':process_assigned_chains,
    'chronology':process_chronology,
    'complex_types':process_complex_types,
    'experiment':process_experiment,
    'peptide_details':process_peptide_details,
    'publication_details':process_publication_details,
    'species':process_species, 
    'tcr_details':process_tcr_details,
    'titles':process_title
}

errors = {}
for facet in facets:
    errors[facet] = []


cores = {}
listings = {}

pdb_codes = ['1hhk']

pdb_codes = load_pdb_lists('class_i', config['WAREHOUSE_PATH'], console)

console.rule(f'[bold]Loading core files for structures')


core_filepath = f'{config["WAREHOUSE_PATH"]}/datacompilations/core.json'
listings_filepath = f'{config["WAREHOUSE_PATH"]}/datacompilations/listings.json'


for pdb_code in pdb_codes:
    core = build_core_data(pdb_code)
    for facet in facets:
        current_facet = load_facet(pdb_code, facet)
        core = process_facet(core, current_facet, facet, facets[facet])
    core = process_urls(core, pdb_code)
    write_facet(pdb_code, 'hydrated', core)

    listing = make_listing(core)
    write_facet(pdb_code, 'listing', listing)
    listings[pdb_code] = listing
    cores[pdb_code] = core
    

#print (cores)
#print (listings)

console.rule(f'[bold]Errors')
for error in errors:
    this_error = errors[error]
    if len(this_error) > 0:
        print (f'{error} : {json.dumps(this_error)}\n')

write_json(core_filepath, cores, verbose=True)
write_json(listings_filepath, listings, verbose=True)
