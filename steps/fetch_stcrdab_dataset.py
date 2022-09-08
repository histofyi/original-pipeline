import pandas as pd


from functions.cli import load_config
from functions.files import write_json

from rich.progress import Progress
from rich.console import Console



console = Console()
config = load_config()

tcr_set = {}

chain_types = {
    'abTCR': [
        {'chain':'alpha','letter':'A'},
        {'chain':'beta','letter':'B'}
    ],
    'gdTCR': [
        {'chain':'gamma','letter':'G'},
        {'chain':'delta','letter':'D'}
    ]
}

error_list = []
pdb_code_list = []

def chain_set_template():
  return {
      'chains':[],
      'subgroup':''
  }


def process_tcr_data(raw_tcr_dict):
  for row in raw_tcr_dict:
    clean_tcr_dict = {}
    if str(row['TCRtype']) != 'nan':
      clean_tcr_dict['pdb_code'] = pdb_code
      chain_type = chain_types[row['TCRtype']]
      for chain_row in chain_type:
        if str(row['mhc_type']) != 'nan':
          if row['mhc_type'] == 'MH1':
            clean_tcr_dict['mhc_type'] = 'class_i'
          elif row['mhc_type'] == 'MH2':
            clean_tcr_dict['mhc_type'] = 'class_ii'
        if not chain_row['chain'] in clean_tcr_dict:
          clean_tcr_dict[chain_row['chain']] = chain_set_template()
        clean_tcr_dict[chain_row['chain']]['chains'].append(row[f'{chain_row["letter"]}chain'])
        if not pd.isna(row[f'{chain_row["chain"]}_subgroup']):
            clean_tcr_dict[chain_row['chain']]['subgroup'] = row[f'{chain_row["chain"]}_subgroup']
        else:
            clean_tcr_dict[chain_row['chain']]['subgroup'] = None
    else:
      clean_tcr_dict = None
    return clean_tcr_dict


def fetch_tcr_info(pdb_code):
  url = f'http://opig.stats.ox.ac.uk/webapps/stcrdab/summary/{pdb_code}'
  print (url)
  raw_tcr_df = pd.read_table(url, sep='\t')
  raw_tcr_dict = pd.concat([
    raw_tcr_df['Achain'], 
    raw_tcr_df['Bchain'], 
    raw_tcr_df['Gchain'], 
    raw_tcr_df['Dchain'], 
    raw_tcr_df['TCRtype'], 
    raw_tcr_df['mhc_type'], 
    raw_tcr_df['antigen_chain'], 
    raw_tcr_df['mhc_chain1'], 
    raw_tcr_df['mhc_chain2'],
    raw_tcr_df['alpha_subgroup'],
    raw_tcr_df['beta_subgroup'],
    raw_tcr_df['gamma_subgroup'],
    raw_tcr_df['delta_subgroup'],
  ], axis = 1).to_dict(orient='records')
  return raw_tcr_dict


console.rule('[bold] Fetching list of PDB codes of T-cell receptors from STCRdab')

df = pd.read_html('http://opig.stats.ox.ac.uk/webapps/stcrdab/Browser?all=true')[0]

pdb_codes = pd.concat([df['PDB']])

i = 0

with Progress() as progress:
        
    task = progress.add_task("[white]Processing...", total=len(pdb_codes))

    for pdb_code in pdb_codes:
        print (pdb_code)
        pdb_code_list.append(pdb_code)
        raw_tcr_dict = fetch_tcr_info(pdb_code)
        clean_tcr_dict = process_tcr_data(raw_tcr_dict)
        if clean_tcr_dict:
            tcr_set[pdb_code] = clean_tcr_dict
            filepath = f'{config["WAREHOUSE_PATH"]}/structures/info/public/tcr_details/{pdb_code}.json'
            write_json(filepath, clean_tcr_dict, verbose=True, pretty=True)
        else:
            error_list.append(pdb_code)
        print (clean_tcr_dict)
        i += 1
        progress.update(task, advance=1)

console.rule('[bold]Tidying up and saving')
filepath = f'{config["CONSTANTS"]}/tcrs.json'
write_json(filepath, tcr_set, verbose=True, pretty=True)
console.print('[bold green]Done')