from typing import Dict, List
import os

from functions.cli import load_config, parse_args
from functions.files import write_json
from functions.core import build_core_data

from rich.progress import Progress
from rich.console import Console




def build_core_for_list(pdb_codes:str, warehouse_path:str, console, force:bool=False):
    """
    This function will build the core information for a set of pdb codes

    Args:
        pdb_codes (str): the pdb codes to be acted upon
        warehouse_path (str): the filepath to the local copy of the warehouse
        force (bool): whether the existing records should be wiped (use with care!)

    """
    record_count = len(pdb_codes)
    records_changed = []
    for pdb_code in pdb_codes:
        filepath = f'{warehouse_path}/structures/info/public/core/{pdb_code}.json'
        # if we're forcing the update of this file (and removing all information from it)
        if force:
            core_info = build_core_data(pdb_code)
        else:
            # check if file exists, if it does, don't do anything
            if not os.path.exists(filepath):
                core_info = build_core_data(pdb_code)
            else:
                core_info = None
        if core_info:
            records_changed.append(pdb_code)
            write_json(filepath, core_info, verbose=True, pretty=True)
    print (f'{len(records_changed)} out of {record_count} have been changed')
    

        

def main():

    config = load_config()

    console.rule('[bold] Initialiseing records')
    pdb_codes, mhc_class, set_slug, set_context, force = parse_args(console)

    build_core_for_list(pdb_codes, config['WAREHOUSE_PATH'], console, force=force)


console = Console()
main()

    




