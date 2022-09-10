import toml
import argparse

from .pdb import load_pdb_lists


def print_spacer():
    print (f'{"-"*68}\n')


def load_config():
    config = toml.load('config.toml')
    return config


def parse_args(console):
    config = load_config()

    parser = argparse.ArgumentParser()
    
    #TODO create help text
    parser.add_argument("--pdb_code")
    parser.add_argument("--mhc_class")
    parser.add_argument("--set_slug")
    parser.add_argument("--set_context")
    parser.add_argument("--force")


    args = parser.parse_args()
    pdb_code = args.pdb_code
    mhc_class = args.mhc_class
    set_slug = args.set_slug
    set_context = args.set_context
    force = args.force

    pdb_codes = []

    if pdb_code:
        console.print(f'[bold] Single PDB code mode - {pdb_code}')
        pdb_codes = [pdb_code]
    else:
        if mhc_class:
            console.print(f'[bold] MHC class mode- {mhc_class}')
            pdb_codes = load_pdb_lists(mhc_class, config['WAREHOUSE_PATH'], console)
        elif set_slug and set_context:
            #TODO fill out this stub
            console.print(f'[bold] Set mode- {set_context}/{set_slug}')
            members = []
            pdb_codes = members

    return pdb_codes, mhc_class, set_slug, set_context, force