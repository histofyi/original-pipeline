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
    with open(filepath) as f:
        json_data = json.load(f)
    return json_data


def write_json(filepath:str, contents:Union[List,Dict], verbose:bool=False, pretty:bool=False):
    with open(filepath, 'w') as outfile:
        if pretty:
            outfile.write(json.dumps(contents, sort_keys=True, indent=4))
        else:
            outfile.write(json.dumps(contents, sort_keys=True))
    if verbose:
        print (f'JSON file \'{filepath}\' written')


def remove_file(filepath:str, verbose:bool=False):
    if os.path.exists(filepath):
        os.remove(filepath)
        if verbose:
            print (f'File {filepath} removed')
    else:
        if verbose:
            print (f'File {filepath} does not exist')

    



