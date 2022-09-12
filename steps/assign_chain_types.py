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
fcrn = []
unmatched_sequences = {}


tcrs = load_constants('tcrs')
antibodies = load_constants('abs') 
chain_type_tests = load_constants('chain_type_tests')
assign_chain_types_overrides = load_constants('assign_chain_types_overrides')


def assign_other_types(chain_labels:List, sequence:str) -> Tuple[str, Dict]:
    best_ratio = 0
    best_match = None
    for chain_type in chain_type_tests:
        if chain_type_tests[chain_type]['sequences']:
            for test in chain_type_tests[chain_type]['sequences']:
                ratio = fuzz.ratio(test, sequence) / 100
                if ratio > best_ratio:
                    best_match = chain_type
                    best_ratio = ratio
    if best_ratio >= 0.7:
        return best_match, {'chains':chain_labels, 'sequence':sequence, 'match_type':'histo:assign_chains', 'score':best_ratio}
    else:
        return None, None


def assign_peptide(chain_labels:List, peptide_details:Dict, sequence:str) -> Tuple[str, Dict]:
    match = False
    for chain_label in peptide_details:
        if chain_label in chain_labels:
            match = True
    if match:
        return 'peptide', {'chains':chain_labels, 'sequence':sequence, 'match_type':'histo:classify_peptide'}
    else:
        return None, None


def assign_tcrs(chain_labels:List, tcr_details:Dict, sequence:str) -> Tuple[str, Dict]:
    tcr_chains = ['alpha','beta','gamma','delta']
    match = False
    match_info = None
    match_tcr_chain = None
    for tcr_chain in tcr_chains:
        if tcr_chain in tcr_details:
            for chain_label in chain_labels:
                if chain_label in tcr_details[tcr_chain]['chains']:
                    match = True
                    match_info = tcr_details[tcr_chain]
                    match_tcr_chain = f'tcr_{tcr_chain}'
                    match_info['sequence'] = sequence
                    match_info['match_type'] = 'stcrdab'
                    return match_tcr_chain, match_info                    
    return None, None


def assign_abs(chain_labels:List, ab_details:Dict, sequence:str) -> Tuple[str, Dict]:
    ab_chains = ['heavy','light']
    match = False
    match_info = None
    match_ab_chain = None
    for ab_chain in ab_chains:
        if ab_chain in ab_details:
            for chain_label in chain_labels:
                if chain_label in ab_details[ab_chain]['chains']:
                    match = True
                    match_info = ab_details[ab_chain]
                    match_ab_chain = f'ab_{ab_chain}'
                    match_info['sequence'] = sequence
                    match_info['match_type'] = 'sabdab'
                    return match_ab_chain, match_info                    
    return None, None


def assign_chain_type(chain_labels:List, peptide_details:Dict, tcr_details:Dict, ab_details:Dict, sequence:str) -> Tuple[str, Dict]:
    if peptide_details:
        chain_type, chain_details = assign_peptide(chain_labels, peptide_details, sequence)
    else:
        chain_type = None
    if not chain_type:        
        if tcr_details:
            chain_type, chain_details = assign_tcrs(chain_labels, tcr_details, sequence)
    if not chain_type:
        if ab_details:
            chain_type, chain_details = assign_abs(chain_labels, ab_details, sequence)
    if not chain_type:
        chain_type, chain_details = assign_other_types(chain_labels, sequence)
    return chain_type, chain_details



def assign_chain_types_for_list(pdb_codes:str, warehouse_path:str, all_chains, console, force:bool=False):
    console.rule('[bold]Assigning chains to chain types')



    with Progress() as progress:
        
        task = progress.add_task("[white]Processing...", total=len(pdb_codes))
    
        for pdb_code in pdb_codes:
            chain_set = {}

            alike_chains = load_facet(pdb_code, 'alike_chains')
            if 'note' in alike_chains:
                del alike_chains['note']
            if 'corrected' in alike_chains:
                alike_chains = alike_chains['corrected']

            peptide_details = load_facet(pdb_code, 'peptide_details')
            
            if pdb_code in tcrs:
                tcr_details = tcrs[pdb_code]
            else:
                tcr_details = None

            if pdb_code in antibodies:
                ab_details = antibodies[pdb_code]
            else:
                ab_details = None



            local_unmatched_sequences = {}
            for chain_label in alike_chains:
                chain = f'{pdb_code}_{chain_label}'
                sequence = all_chains.loc[chain]['sequence']
                
                chain_type, chain_details = assign_chain_type(alike_chains[chain_label], peptide_details, tcr_details, ab_details, sequence)
                
                if chain_type:
                    chain_set[chain_type] = chain_details
                else:
                    local_unmatched_sequences[chain_label] = sequence


            
            if len(chain_set) ==  len(alike_chains):
                fully_matched.append(pdb_code)
            else:
                if pdb_code in assign_chain_types_overrides:
                    print (pdb_code)
 
                    chain_overrides = assign_chain_types_overrides[pdb_code]

                    for chain in chain_overrides:
                        not_matched = False
                        if chain != 'note':
                            if chain in chain_set:
                                if chain_set[chain]['chains'] != chain_overrides[chain]['chains']:
                                    not_matched = True
                            else:
                                not_matched = True
                        if not_matched:
                            if pdb_code not in override_used:
                                override_used.append(pdb_code)
                            chain_set[chain] = chain_overrides[chain]
                            chain_set[chain]['match_type'] = 'histo:override'
                            chain_label = chain_overrides[chain]['chains'][0]
                            if not 'sequence' in chain_set[chain]:
                                if chain_label in local_unmatched_sequences:
                                    chain_set[chain]['sequence'] = local_unmatched_sequences[chain_label]
                                else:
                                    try:
                                        localpdb_chain = f'{pdb_code}_{chain_label}'
                                        chain_set[chain]['sequence'] = all_chains.loc[localpdb_chain]['sequence']
                                    except:
                                        print (f'sequence_retrieval_error for {localpdb_chain}')
                    if 'peptide_fragment1' in chain_set:
                        del chain_set['peptide']
                            
                    chain_set['note'] = chain_overrides['note']
                    print (f'Revised chain_set for {pdb_code}:')
                    for part in chain_set:
                        print (f'\n{part}:-')
                        if part != 'note':
                            for element in chain_set[part]:
                                print (f'{element} : {chain_set[part][element]}')
                        else:
                            print (chain_set[part])
                    print_spacer()

                partially_matched.append(pdb_code)

            write_facet(pdb_code, 'assigned_chains', chain_set)
            progress.update(task, advance=1)
    
    pass


def main():

    config = load_config()
    
    console.rule('[bold] Assigning chains to chain types')

    pdb_codes, mhc_class, set_slug, set_context, force = parse_args(console)

    localpdb_path = config['LOCALPDB_PATH']

    console.rule(f'[bold]Matching pdb_codes with localpdb')

    with console.status(f'Spinning up localpdb...'):
        lpdb = PDB(db_path=localpdb_path)
        all_chains = lpdb.chains

    assign_chain_types_for_list(pdb_codes, config['WAREHOUSE_PATH'], all_chains, console, force=force)

    console.rule(f'[bold]Statistics')

    print (f'Fully matched: {len(fully_matched)}')
    print (f'Partially matched: {len(partially_matched)}')
    print (f'Override used: {len(override_used)}')
    print (f'Success: {len(override_used) + len(fully_matched)} out of {len(pdb_codes)}')
    #if len(partially_matched) > 0:
    #   print (partially_matched)
    
    #for pdb_code in unmatched_sequences:
    #    print (pdb_code)
    #    print (unmatched_sequences[pdb_code])
    #    print ('')
    console.print('[bold green]Done')

console = Console()
main()
