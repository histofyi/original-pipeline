from typing import Dict, List, Tuple

from localpdb import PDB
from fuzzywuzzy import fuzz

import toml
import json

from functions.files import load_constants, write_json, remove_file
from functions.cli import print_spacer

from rich.progress import Progress
from rich.console import Console

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--search_set")



console = Console()







def setup_query(search_set:str='class_i') -> Tuple[Dict, List, List]:
    """


    """
    
    chain_query_info = load_constants('chain_queries')[search_set]
    match_sequence_set = load_constants('search_sequences')[chain_query_info['sequences']]
    match_sequences = match_sequence_set['sequences']
    match_labels = match_sequence_set['labels']

    print (f'\nQuerying localpdb for {search_set}\n')
    return chain_query_info, match_sequences, match_labels


def match_chain(sequence:str, match_sequences:List, match_labels:List) -> Dict:
    i = 0
    hit = False
    best_ratio = 0
    best_match = None
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


def run_query(localpdb_path:str, query_info:Dict, match_sequences:List, match_labels:List) -> Tuple[Dict,List]:
    low = query_info['low']
    high = query_info['high']

    possible_hits = {}
    searched = []

    console.rule(f'[bold]Searching localpdb  - {query_info["label"]}')
    with console.status(f'Searching for length matches... {low} to {high} residues'):
        lpdb = PDB(db_path=localpdb_path)
        results = lpdb.chains.query(f'{low} < sequence.str.len() < {high}')


    results_count = len(results)
    print (f'\n{results_count} matching chain length of {low} to {high} residues\n')

    n = 0

    console.rule('[bold]Matching')

    with Progress() as progress:
        
        task = progress.add_task("[white]Processing...", total=results_count)

        for structure in results.index:
            
            progress.update(task, advance=1)

            pdb_code = results.loc[structure]['pdb']
            chain = structure.split('_')[1]
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
                    if pdb_code not in possible_hits:
                        possible_hits[pdb_code] = {
                            'chains':[],
                            'best_ratio':0,
                            'best_match':''
                        }
                    possible_hits[pdb_code]['chains'].append(chain)
                    if match_stats['best_ratio'] > possible_hits[pdb_code]['best_ratio']:
                        possible_hits[pdb_code]['best_ratio'] = match_stats['best_ratio']
                        possible_hits[pdb_code]['best_match'] = match_stats['best_match']
                    
                    checkpoint_match_save(query_info, possible_hits)

                searched.append(structure)
                n += 1
                if n == 10:
                    n = 0
                    checkpoint_save(query_info, searched)
    return possible_hits, searched

    




def main():
    args = parser.parse_args()

    search_set = args.search_set

    config = toml.load('config.toml')

 
    query_info, match_sequences, match_labels = setup_query(search_set)
    possible_hits, searched = run_query(config['LOCALPDB_PATH'], query_info, match_sequences, match_labels)
    
    print('')
    console.rule('[bold]Stats and saving output')
    print (f'{len(possible_hits)} matches found from {len(searched)} sequence length search results\n')
    
    write_json(f'{config["OUTPUT_PATH"]}/queries/{query_info["label"]}_hits.json', possible_hits, verbose=True, pretty=True)
    
    write_json(f'{config["OUTPUT_PATH"]}/queries/{query_info["label"]}_searched.json', searched, verbose=True)
    
    console.rule('[bold]Tidying up')
    print ('Cleaning up tmp files\n')

    filepath = f'tmp/{query_info["label"]}_checkpoint.json'
    remove_file(filepath, verbose=True)

    filepath = f'tmp/{query_info["label"]}_match_checkpoint.json'
    remove_file(filepath, verbose=True)

    console.print('[bold green]Done')
    

main()