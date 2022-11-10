
from typing import Dict, List, Tuple

from functions.files import load_facet, write_file, write_facet
from functions.cli import load_config
from functions.pdb import load_pdb_lists

from common.providers import httpProvider

from rich.progress import Progress
from rich.console import Console

import Bio.PDB
from Bio.PDB.mmcifio import MMCIFIO
from io import StringIO, TextIOWrapper
from Bio.PDB.MMCIFParser import MMCIFParser

import os
from datetime import datetime

config = load_config()
console = Console()

errors = []

def download_cif_coordinates(pdb_code:str, assembly_id:int) -> str:
    """
    This function
    """
    pdbe_url = f'https://www.ebi.ac.uk/pdbe/model-server/v1/{pdb_code}/assembly?name={assembly_id}&model_nums=1&encoding=cif&copy_all_categories=false&download=false'
    coordinates = httpProvider().get(pdbe_url, format='txt')
    return coordinates


def get_assembly_count(assigned_chains:Dict) -> int:
    chain_counts = []
    total_count = 0
    for chain in assigned_chains:
        if chain not in ['note','force_override']:
            chain_counts.append(len(assigned_chains[chain]['chains']))
            total_count += len(assigned_chains[chain]['chains'])

    processed_chain_counts = sorted([count for count in set(chain_counts)])

    if len(processed_chain_counts) == 1:
        assembly_count = processed_chain_counts[0]
    else:
        lowest_count = processed_chain_counts[0]
        successes = []
        for count in processed_chain_counts:
            if count % lowest_count == 0:
                successes.append(True)
        if len([success for success in set(successes)]) != 1:
            print ('MISMATCH')
            assembly_count = None
        else:
            assembly_count = lowest_count
    if total_count % assembly_count == 0:
        chain_count = int(total_count/assembly_count)
    else:
        print ('INCORRECT TOTAL CHAINS VS ASSEMBLIES')
        chain_count = None
    return assembly_count, chain_count 


def get_all_chain_letters(assigned_chains:Dict) -> List:
    chain_letters = []
    if assigned_chains is not None:
        for chain in assigned_chains:
            for individual_chain in assigned_chains[chain]['chains']:
                chain_letters.append(individual_chain)
    return chain_letters


def clean_assigned_chains(assigned_chains:Dict) -> Dict:
    if assigned_chains is not None:
        cleaned_assigned_chains = {}
        for chain in assigned_chains:
            if chain not in ['note','force_override']:
                cleaned_assigned_chains[chain] = assigned_chains[chain]
    else:
        cleaned_assigned_chains = None
    return cleaned_assigned_chains


pdb_codes = ['1hhk']

pdb_codes = load_pdb_lists('class_i', config['WAREHOUSE_PATH'], console)

console.rule(f'[bold]Loading raw files for structures')

all_structures_count = len(pdb_codes)

chain_number_mismatch = []
different_chain_sets = []
parsing_error = []
no_assigned_chains = []

with Progress() as progress:
        
    task = progress.add_task("[white]Processing...", total=all_structures_count)

    for pdb_code in pdb_codes:

        assigned_chains = clean_assigned_chains(load_facet(pdb_code, 'assigned_chains'))

        if assigned_chains is not None:
            to_write = False
            assemblies = {}

            assembly_count, chain_count = get_assembly_count(assigned_chains)
            chain_letters = get_all_chain_letters(assigned_chains)

            assembly_id = 1
            while assembly_id <= assembly_count:
                assembly_name = f'{pdb_code}_{assembly_id}'
                
                

                if chain_count is None:
                    print (pdb_code)
                    print (assembly_count)
                    print (chain_count)
                    errors.append(assembly_name)
                    filepath = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/with_solvent/errors/{assembly_name}.cif'
                else:
                    filepath = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/with_solvent/raw/{assembly_name}.cif'

                    if not os.path.exists(filepath):
                        coordinates = download_cif_coordinates(pdb_code, assembly_id)    
                        cif_file = StringIO(coordinates)
                        parser = MMCIFParser(QUIET=True)
                        try:
                            structure = parser.get_structure(assembly_name, cif_file)
                            chains = [chain.id for chain in structure.get_chains()]
                            print (chains)
                            write_file(filepath, coordinates)
                            assemblies[assembly_id] = {
                                'chains':chains,
                                'downloaded_at':datetime.now().isoformat(),
                                'assembly_name':assembly_name
                            }
                            to_write = True
                        except:
                            print ('PARSING ERROR')
                            errors.append(assembly_name)
                assembly_id += 1
            if to_write:
                print (assemblies)
                write_facet(pdb_code, 'assemblies', assemblies)


                



        progress.update(task, advance=1)


print (f'{errors=}')
print (f'{chain_number_mismatch=}')
print (f'{different_chain_sets=}')
print (f'{parsing_error=}')
print (f'{no_assigned_chains=}')
