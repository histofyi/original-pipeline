from typing import Dict, List

from common.providers import httpProvider
from rich.progress import Progress
from rich.console import Console

from functions.pdb import load_pdb_lists
from functions.cli import load_config
from functions.files import write_json

import json

config = load_config()
console = Console()

def rank_data_by_field(field:str, ascending:bool=True) -> List:
    if ascending:
        asc_desc = 'asc'
    else:
        asc_desc = 'desc'
    print (field)
    print (asc_desc)
    baseurl = 'http://127.0.0.1:8001/core.json?sql='
    select = 'select+pdb_code+from+core+'
    if field is not None:
        if field in ['deposition_date','pdb_code','resolution','revision_date','release_date']:
            
            order_by = f'order+by+{field}+{asc_desc}'
            url = f'{baseurl}{select}{order_by}'
            raw_return = httpProvider().get(url, 'json')
            pdb_codes = parse_datasette_return(raw_return['rows'])
            print (f'{len(pdb_codes)} in index')
            return pdb_codes
    else:
        print ('Field is None')
        return []


def parse_datasette_return(raw_return:Dict) -> List:
    pdb_codes = [pdb_code[0] for pdb_code in raw_return]
    return pdb_codes 



pdb_codes = load_pdb_lists('class_i', config['WAREHOUSE_PATH'], console)

print (f'{len(pdb_codes)} records in the current dataset')

indexes = {}

for field in ['pdb_code','resolution','deposition_date','revision_date','release_date']:
    index_name = f'{field}_asc'
    indexes[index_name] = rank_data_by_field(field)
    index_name = f'{field}_desc'
    indexes[index_name] = rank_data_by_field(field, ascending=False)


stringified_pdb_codes = spaced_stringify(pdb_codes)

for index in indexes:
    difference  = set(pdb_codes).difference(set(indexes[index]))
    if len(difference) > 0:
        print ('Reload the core.db with any updated data')
    else:
        print ('Matching sets')



index_filepath = f'{config["WAREHOUSE_PATH"]}/datacompilations/index.json'
write_json(index_filepath, indexes, verbose=True)




