from typing import Dict
import pandas as pd


from functions.cli import load_config, parse_args
from functions.files import write_json


from rich.progress import Progress
from rich.console import Console



console = Console()
config = load_config()

ab_set = {}

error_list = []
ab_pdb_codes = []


def fetch_ab_info(pdb_code):
    url = f'http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/summary/{pdb_code}'
    raw_ab_df = pd.read_table(url, sep='\t')
    raw_ab_dict = pd.concat([
        raw_ab_df['Hchain'], 
        raw_ab_df['Lchain'], 
        raw_ab_df['heavy_species'], 
        raw_ab_df['light_species'], 
        raw_ab_df['heavy_subclass'], 
        raw_ab_df['light_subclass'], 
        raw_ab_df['light_ctype'],
        raw_ab_df['scfv']
    ], axis = 1).to_dict(orient='records')
    return raw_ab_dict


def process_ab_dict(raw_ab_dict:Dict):
    i = 0
    for row in raw_ab_dict:
        if i == 0:
            chain_info = {
                    'light':{
                        'sabdab':'Lchain',
                        'chains':[]
                    },
                    'heavy':{
                        'sabdab':'Hchain',
                        'chains':[]
                    }
                }
            for item in ['heavy_species','light_species', 'heavy_subclass','light_subclass','light_ctype','scfv']:
                chain_info[item] = row[item]
        for chain in [('heavy','Hchain'),('light','Lchain')]:
            if chain[1] in row:
                if row[chain[1]]:
                    chain_info[chain[0]]['chains'].append(row[chain[1]])

    print (chain_info)
    return chain_info



config = load_config()
    

pdb_codes, mhc_class, set_slug, set_context, force = parse_args(console)


console.rule('[bold] Fetching list of PDB codes of antibodies with an MHC molecule')

url = 'http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/search/?all=true'

df = pd.read_html(url)[0]

ab_pdb_codes = pd.concat([df['PDB']])

print (len(ab_pdb_codes))

i = 0

with Progress() as progress:
        
    task = progress.add_task("[white]Processing...", total=len(ab_pdb_codes))

    for pdb_code in ab_pdb_codes:
        if pdb_code in pdb_codes:
            print (pdb_code)

            raw_ab_dict = fetch_ab_info(pdb_code)
            chain_info = process_ab_dict(raw_ab_dict)
            ab_set[pdb_code] = chain_info

            i += 1
        progress.update(task, advance=1)

console.rule('[bold]Statistics')

console.rule('[bold]Tidying up and saving')
filepath = f'{config["CONSTANTS"]}/abs.json'
write_json(filepath, ab_set, verbose=True, pretty=True)
console.print('[bold green]Done')