from typing import Dict, Optional
import toml



def load_config() -> Dict:
    """ 
    Loads the configuration file for the pipline and returns a dictionary of values

    Returns:
        dict: a dictionary of configuration variables 

    """
    config = toml.load('config.toml')
    return config


def read_lockfile(filename:str) -> Optional[str]:
    """
    Reads the content of a lockfile
    
    Args:
        filename (str) - the filename of the lockfile

    Returns:
        str: the contents of the 
    """
    try:
        with open(filename) as lockfile:
            test = lockfile.read()
    except:
        write_lockfile(filename, 'intialised')
        test = 'initialised'
    return test


def write_lockfile(filename:str, text:str):
    """
    Writes the content of a lockfile

    Args:
        filename (str) - the filename of the lockfile
        text (str) - the contents for the lockfile
    """
    with open(filename, 'w') as lockfile:
        lockfile.write(text)
    pass