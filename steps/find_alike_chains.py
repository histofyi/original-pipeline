from typing import Dict, List, Tuple

import os
from fuzzywuzzy import fuzz

from functions.cli import print_spacer
from functions.files import write_json, read_json, load_constants
from functions.cli import load_config, parse_args

from localpdb import PDB

from rich.progress import Progress
from rich.console import Console


console = Console()
config = load_config


fully_matched = []
partially_matched = []
only_one_assembly = []
override_used = []


def exact_match(sequence:str, chain:str, test:str, test_chain:str, verbose=False) -> bool:
    if sequence == test:
        if verbose:
            print (f'Exact match between {chain} and {test_chain}')
            print (sequence)
            print (test)  
            print ('')          
        return True
    else:
        return False


def fuzzy_match(sequence:str, chain:str, test:str, test_chain:str, verbose=False, cutoff=0.9) -> bool:
    length_diff = abs(len(sequence) - len(test))
    if length_diff > 10:
        return False
    else:
        if verbose:
            ratio = fuzz.ratio(test, sequence) / 100
            if ratio > cutoff:
                return True
            else:
                return False


def find_alike_chains(pdb_code:str, all_chains, overrides:Dict):
    structure_chains = [chain for chain in all_chains.loc[all_chains['pdb'] == pdb_code].index]
    structure_sequences = [all_chains.loc[chain]['sequence'] for chain in structure_chains]
    chain_labels = [chain.split('_')[1] for chain in structure_chains]
    i = 0
    matched_already = []
    chain_sets = {}
    unique_structure_sequences = [sequence for sequence in set(structure_sequences)]
    if pdb_code in overrides:
        chain_sets = overrides[pdb_code]
        override_used.append(pdb_code)
    else:
        if len(structure_chains) == len(unique_structure_sequences):
            only_one_assembly.append(pdb_code)
            fully_matched.append(pdb_code)
            for chain_label in chain_labels:
                chain_sets[chain_label] = [chain_label]
        else:
            for chain_label in chain_labels:
                sequence = structure_sequences[i]
                j = 0
                if not chain_labels[i] in matched_already:
                    for test_sequence in structure_sequences:
                        if i != j:
                            if exact_match(sequence, chain_label, structure_sequences[j], chain_labels[j], verbose=False):
                                if not chain_label in chain_sets:
                                    chain_sets[chain_label] = [chain_label]
                                chain_sets[chain_label].append(chain_labels[j])
                                if not chain_label in matched_already:
                                    matched_already.append(chain_label)
                                if not chain_labels[j] in matched_already:
                                    matched_already.append(chain_labels[j])
                            else:
                                if fuzzy_match(sequence, chain_label, structure_sequences[j], chain_labels[j], verbose=True):
                                    if not chain_label in chain_sets:
                                        chain_sets[chain_label] = [chain_label]
                                    chain_sets[chain_label].append(chain_labels[j])
                                    if not chain_label in matched_already:
                                        matched_already.append(chain_label)
                                    if not chain_labels[j] in matched_already:
                                        matched_already.append(chain_labels[j])
                                    
                        j += 1
                i += 1
            if len(matched_already) == len(structure_chains):
                fully_matched.append(pdb_code)
            else:
                partially_matched.append(pdb_code)
    print (pdb_code)
    print (chain_sets)
    return chain_sets
    


def find_alike_chains_for_list(pdb_codes:str, warehouse_path:str, all_chains, console, force:bool=False):
    """
    This function will build the core information for a set of pdb codes

    Args:
        pdb_codes (str): the pdb codes to be acted upon
        warehouse_path (str): the filepath to the local copy of the warehouse
        force (bool): whether the existing records should be wiped (use with care!)

    """
    record_count = len(pdb_codes)
    records_changed = []
    overrides = load_constants('alike_chains_overrides')

    console.rule('[bold]Matching alike chains')

    with Progress() as progress:
        
        task = progress.add_task("[white]Processing...", total=len(pdb_codes))

        for pdb_code in pdb_codes:
            filepath = f'{warehouse_path}/structures/info/public/alike_chains/{pdb_code}.json'
            # if we're forcing the update of this file (and removing all information from it)
            if force:
                alike_chains = find_alike_chains(pdb_code, all_chains, overrides)
            else:
                # check if file exists, if it does, don't do anything
                if not os.path.exists(filepath):
                    alike_chains = find_alike_chains(pdb_code, all_chains, overrides)
                else:
                    alike_chains = None
            if alike_chains:
                records_changed.append(pdb_code)
                write_json(filepath, alike_chains, verbose=True)
            
            progress.update(task, advance=1)
    
    print (f'{len(records_changed)} out of {record_count} have been changed')
    

        

def main():

    config = load_config()
    
    console.rule('[bold] Finding alike chains in structures')

    pdb_codes, mhc_class, set_slug, set_context, force = parse_args(console)

    localpdb_path = config['LOCALPDB_PATH']

    console.rule(f'[bold]Matching pdb_codes with localpdb')

    with console.status(f'Spinning up localpdb...'):
        lpdb = PDB(db_path=localpdb_path)
        all_chains = lpdb.chains

    find_alike_chains_for_list(pdb_codes, config['WAREHOUSE_PATH'], all_chains, console, force=force)

    console.rule(f'[bold]Statistics')

    print (f'Fully matched: {len(fully_matched)}')
    print (f'Partially matched: {len(partially_matched)}')
    print (f'Override used: {len(override_used)}')
    print (f'Success: {len(override_used) + len(fully_matched)} out of {len(pdb_codes)}')
    if len(partially_matched) > 0:
        print (partially_matched)
    
    console.print('[bold green]Done')

console = Console()
main()
