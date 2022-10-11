from typing import List

import gzip
from Bio.PDB.MMCIFParser import MMCIFParser

parser = MMCIFParser(QUIET=True)

from .files import read_json, load_constants


def load_structure(pdb_code:str, pdb_filepath:str):
    print (pdb_filepath)
    try:
        with gzip.open(pdb_filepath, 'rt') as f:
            structure = parser.get_structure(pdb_code, f)[0]
    except:
        structure = None  
    return structure


def load_pdb_lists(mhc_class:str, warehouse_path:str, console) -> List:
    console.rule(f'[bold]Loading pdb code lists for - {mhc_class}')
    exclude = load_constants('exclude')
    pdb_codes = []
    for structure_type in ['class_i','truncated_class_i','full_length_class_i', 'single_chain_class_i']:
        print (f'Loading {structure_type}')
        filepath = f'{warehouse_path}/queries/{structure_type}_hits.json'
        json_data = read_json(filepath)
        pdb_codes += [pdb_code for pdb_code in json_data]
        pdb_codes = [pdb_code for pdb_code in set(pdb_codes)]
        print (len(pdb_codes))
    pdb_codes = [pdb_code for pdb_code in pdb_codes if pdb_code not in exclude]
    return pdb_codes
