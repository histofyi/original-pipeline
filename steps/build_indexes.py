from typing import Dict, List

from common.providers import httpProvider



def rank_data_by_field(field:str, ascending:bool=True) -> List:

    baseurl = 'http://127.0.0.1:8001/core.json?sql='
    select = 'select+pdb_code+from+core+'
    if field is not None:
        if field in ['deposition_date','pdb_code','resolution','revision_date','release_date']:
            if ascending:
                asc_desc = 'asc'
            else:
                asc_desc = 'desc'
            order_by = f'order+by+{field}+{asc_desc}'
            url = f'{baseurl}{select}{order_by}'
            print (url)
            raw_return = httpProvider().get(url, 'json')
            pdb_codes = parse_datasette_return(raw_return['rows'])
            return pdb_codes
    else:
        print ('Field is None')
        return []


def parse_datasette_return(raw_return:Dict) -> List:
    pdb_codes = [pdb_code[0] for pdb_code in raw_return]
    print (len(pdb_codes))
    return pdb_codes 

indexes = {}

for field in ['pdb_code','resolution','deposition_date','revision_date','release_date']:
    index_name = f'{field}_asc'
    indexes[index_name] = rank_data_by_field(field)
    index_name = f'{field}_desc'
    indexes[index_name] = rank_data_by_field(field, ascending=False)
    
    
print (indexes)




