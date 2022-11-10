from typing import Dict, List, Tuple

import os
from fuzzywuzzy import fuzz

from functions.cli import print_spacer
from functions.files import write_json, read_json, load_constants, load_facet, write_facet
from functions.cli import load_config, parse_args

from localpdb import PDB

from rich.progress import Progress
from rich.console import Console

console = Console()
config = load_config


fully_matched = []
partially_matched = []
override_used = []

def assign_complex_type(pdb_code:str, complex_types:Dict) -> Dict:
    print (pdb_code)
    assigned_chains = load_facet(pdb_code, 'assigned_chains')

    unique_chains = [chain for chain in assigned_chains if chain != 'note']
    unique_chain_count = len(unique_chains)

    possible_matches = complex_types[str(unique_chain_count)]
    #except:
    #    print ('NO POSSIBLE MATCHES')
    #   print (unique_chains)
    #   print (unique_chain_count)
    #    possible_matches = None

    best_match = None
    if possible_matches:
        for match in possible_matches:
            match_count = 0
            matches = []

            for chain in unique_chains:
                if chain in match['components']:
                    match_count += 1
                    matches.append(chain)
            if match_count == unique_chain_count:
                best_match = match
                break
    if not best_match:
        print ('No exact match')
        print (unique_chains)
        print_spacer()
    else:
        print (best_match)
        print_spacer()
        write_facet(pdb_code, 'complex_types', best_match)
        fully_matched.append(pdb_code)
    return {}



def assign_complex_types_for_list(pdb_codes:List, warehouse_path:str, console, force:bool):
    """
    This function takes a set of PDB codes and assigns a complex type based upon the chain assignments

    """
    complex_types = load_constants('complex_types')

    with Progress() as progress:
        
        task = progress.add_task("[white]Processing...", total=len(pdb_codes))

        print ('')
        for pdb_code in pdb_codes:
            
            assign_complex_type(pdb_code, complex_types)

            progress.update(task, advance=1)



        


def main():

    config = load_config()
    
    console.rule('[bold] Assigning structures to complex types')

    pdb_codes, mhc_class, set_slug, set_context, force = parse_args(console)

    print (pdb_codes)

    assign_complex_types_for_list(pdb_codes, config['WAREHOUSE_PATH'], console, force=force)

    console.rule(f'[bold]Statistics')

    print (f'Fully matched: {len(fully_matched)}')
    print (f'Partially matched: {len(partially_matched)}')
    print (f'Override used: {len(override_used)}')
    print (f'Success: {len(override_used) + len(fully_matched)} out of {len(pdb_codes)}')
    console.print('[bold green]Done')

console = Console()
main()
