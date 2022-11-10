from typing import Dict, List, Tuple

from localpdb import PDB
from fuzzywuzzy import fuzz

import toml
import json

from functions.files import load_constants, write_json, remove_file, read_json
from functions.cli import print_spacer

from rich.progress import Progress
from rich.console import Console

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--mhc_class")



console = Console()

class_i_starts = load_constants('mhc_starts')['class_i']['alpha']

def match_chain(sequence:str, match_sequences:List, match_labels:List) -> Dict:
    i = 0
    hit = False
    best_ratio = 0
    best_match = None
    
    start_index = None
    for class_i_start in class_i_starts:
        if class_i_start in sequence:
            start_index = sequence.index(class_i_start)

    if start_index:
        sequence = sequence[start_index:]
        if len(sequence) > 275:
            sequence = sequence[:275]

    for test in match_sequences:
        if len(test) > len(sequence):
            test = test[0:len(sequence) - 1]
        ratio = fuzz.ratio(test, sequence) / 100
        if ratio > 0.7:
            hit = True
        if ratio > best_ratio:
            best_ratio = round(ratio,3)
            best_match = match_labels[i]
        i += 1
    return {
        'hit':hit,
        'best_ratio':best_ratio,
        'best_match': best_match
    }


def checkpoint_save(query_info:Dict, searched:List):
    filepath = f'tmp/{query_info["label"]}_checkpoint.json'
    write_json(filepath, searched)


def checkpoint_match_save(query_info:Dict, possible_hits:Dict):
    filepath = f'tmp/{query_info["label"]}_match_checkpoint.json'
    write_json(filepath, possible_hits)


def run_query(lpdb, query_info:Dict, match_sequences:List, match_labels:List, hits_filepath:str, searched_filepath:str) -> Tuple[Dict,List]:
    low = query_info['low']
    high = query_info['high']

    previous_hits = read_json(hits_filepath)
    searched = read_json(searched_filepath)

    print (len(searched))

    possible_hits = {}
    new_hits = {}

    console.rule(f'[bold]Searching localpdb  - {query_info["label"]}')
    with console.status(f'Searching for length matches... {low} to {high} residues'):
        results = lpdb.chains.query(f'{low} < sequence.str.len() < {high}')


    results_count = len(results)
    print (f'\n{results_count} matching chain length of {low} to {high} residues\n')

    n = 0

    console.rule('[bold]Matching')

    with Progress() as progress:
        
        task = progress.add_task(f'[white]Processing... {query_info["label"]} [{results_count} length matches]', total=results_count)

        for structure in results.index:
            
            pdb_code = results.loc[structure]['pdb'].lower()
            chain = structure.split('_')[1]
                
            if pdb_code in previous_hits:
                if pdb_code not in possible_hits:
                    if chain in previous_hits[pdb_code]['chains']:
                        possible_hits[pdb_code] = previous_hits[pdb_code]
                        print_spacer()
                        print (structure)
                        print ('Existing match')
                        print (results.loc[structure])
                        print ('')
            else:
                sequence = results.loc[structure]['sequence']

                if structure not in searched:
                    match_stats = match_chain(sequence, match_sequences, match_labels)

                    if match_stats['hit']:
                        print_spacer()
                        print (structure)
                        print (f'Ratio: {match_stats["best_ratio"]}')
                        print (f'Test: {match_stats["best_match"]}')
                        print (results.loc[structure])
                        print ('')
                        if pdb_code not in new_hits:
                            new_hits[pdb_code] = {
                                'chains':[],
                                'best_ratio':0,
                                'best_match':''
                            }
                        new_hits[pdb_code]['chains'].append(chain)
                        if match_stats['best_ratio'] > new_hits[pdb_code]['best_ratio']:
                            new_hits[pdb_code]['best_ratio'] = match_stats['best_ratio']
                            new_hits[pdb_code]['best_match'] = match_stats['best_match']
                        
                        checkpoint_match_save(query_info, new_hits)


            progress.update(task, advance=1)
            searched.append(structure)
            n += 1
            if n == 10:
                n = 0
                checkpoint_save(query_info, searched)
    for pdb_code in new_hits:
        possible_hits[pdb_code] = new_hits[pdb_code]
    return possible_hits, new_hits, searched

    
def run_mhc_class(mhc_class:str, config:Dict):
    """
    This function runs the query on a variety of different query lengths. 

    This enables the detection of normal length sequences, truncated sequences, full length sequences and single chain constructs

    Args:
        mhc_class (str): the class of MHC molecule to search for e.g. class_i (only class_i currently supported)
        config (dict): the configuration data e.g the filepath of specific directories
    """
    chain_query_info = load_constants('chain_queries')[mhc_class]
    match_sequence_set = load_constants('search_sequences')[mhc_class]
    match_sequences = match_sequence_set['sequences']
    match_labels = match_sequence_set['labels']
    version_file = f'{config["OUTPUT_PATH"]}/current_version.json'

    version = read_json(version_file)['version']

    lpdb = PDB(db_path=config['LOCALPDB_PATH'])

    all_new = {}

    for search_set in chain_query_info:
        query_info = chain_query_info[search_set]

        hits_filepath = f'{config["OUTPUT_PATH"]}/queries/{query_info["label"]}_hits.json'
        searched_filepath = f'{config["OUTPUT_PATH"]}/queries/{query_info["label"]}_searched.json'

        print (f'\nQuerying localpdb for {search_set}\n')
   
        possible_hits, new_hits, searched = run_query(lpdb, query_info, match_sequences, match_labels, hits_filepath, searched_filepath)
        
        if len(new_hits) > 0:
            print (new_hits)
            new_hits_filepath = f'{config["OUTPUT_PATH"]}/queries/{query_info["label"]}_new_{version}.json'
            write_json(new_hits_filepath, new_hits, verbose=True, pretty=True)
            for pdb_code in new_hits:
                all_new[pdb_code] = new_hits[pdb_code]

        print('')
        console.rule('[bold]Stats and saving output')
        print (f'{len(possible_hits)} matches found from {len(searched)} sequence length search results\n')
        
        write_json(hits_filepath, possible_hits, verbose=True, pretty=True)
        
        write_json(searched_filepath, searched, verbose=True)
        
        console.rule('[bold]Tidying up')
        print ('Cleaning up tmp files\n')

        filepath = f'tmp/{query_info["label"]}_checkpoint.json'
        remove_file(filepath, verbose=True)

        filepath = f'tmp/{query_info["label"]}_match_checkpoint.json'
        remove_file(filepath, verbose=True)
    if len(all_new) > 0:
        all_new_filepath = f'{config["OUTPUT_PATH"]}/queries/{mhc_class}_new_{version}.json'
        write_json(all_new_filepath, all_new, verbose=True, pretty=True)


def main():
    args = parser.parse_args()

    mhc_class = args.mhc_class

    config = toml.load('config.toml')

    console.rule('[bold]Loading constants and localpdb')
    run_mhc_class(mhc_class, config)

    console.print('[bold green]Done')
    

main()