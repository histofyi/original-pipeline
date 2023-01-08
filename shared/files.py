from typing import Dict, List, Tuple, Union
import json
import toml
import os



config = toml.load('config.toml')



def load_constants(constants_name:str) -> Dict:
    """
    This function returns a dictionary of constants contained in a constants file in the constants directory specified in the config
    
    Args:
        constants_name (str): the name of the constantts dictionary

    Returns:
        Dict: the dictionary of constants
    
    """
    filepath = f'{config["CONSTANTS"]}/{constants_name}.json'
    json_data = read_json(filepath)
    return json_data


def load_facet(pdb_code:str, facet_name:str) -> Dict:
    """
    This function returns a facet (such as 'peptide_details') for a specific pdb code

    Args:
        pdb_code (str): the pdb code in question e.g. 1hhk
        facet_name (str): the required facet e.g. alike_chains
    """
    filepath = f'{config["WAREHOUSE_PATH"]}/structures/info/public/{facet_name}/{pdb_code}.json'
    try:
        json_data = read_json(filepath)
    except:
        json_data = None
    return json_data


def write_facet(pdb_code:str, facet_name:str, contents:Dict):
    filepath = f'{config["WAREHOUSE_PATH"]}/structures/info/public/{facet_name}/{pdb_code}.json'
    write_json(filepath, contents, verbose=True, pretty=True)


def read_json(filepath:str) -> Union[Dict,List]:
    with open(filepath, 'r') as infile:
        json_data = json.load(infile)
    return json_data
    


def write_json(filepath:str, contents:Union[List,Dict], verbose:bool=False, pretty:bool=False):
    with open(filepath, 'w') as outfile:
        if pretty:
            outfile.write(json.dumps(contents, sort_keys=True, indent=4))
        else:
            outfile.write(json.dumps(contents, sort_keys=True))
    if verbose:
        print (f'JSON file \'{filepath}\' written')


def write_file(filepath:str, contents:str, verbose:bool=False):
    with open(filepath,'w') as outfile:
        outfile.write(contents)
    if verbose:
        print (f'String file \'{filepath}\' written')
    

def remove_file(filepath:str, verbose:bool=False):
    if os.path.exists(filepath):
        os.remove(filepath)
        if verbose:
            print (f'File {filepath} removed')
    else:
        if verbose:
            print (f'File {filepath} does not exist')

    

