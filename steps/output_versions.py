import watermark.watermark as watermark
from typing import List, Dict

from functions.helpers import slugify
from functions.files import write_json
from functions.cli import load_config

package_list = "numpy,localpdb,biopandas,pymol,pandas,scipy,fuzzywuzzy,algoliasearch,datasette"


def load_modules(package_list:str):
    """
    This function loads the molecules specified in the list
    """
    packages = package_list.split(',')

    for module in packages:
        module_obj = __import__(module)
        globals()[module] = module_obj



def get_info(package_list=None):
    if package_list:
        info = watermark(packages=package_list)
    else:
        info = watermark()
    info = output_to_dict(info)
    return info


def output_to_dict(info:str) -> Dict:
    dict_info = {}
    for row in info.split('\n'):
        if len(row) > 0:
            items = row.split(': ')
            dict_info[slugify(items[0].strip())] = items[1].strip()
    return dict_info



def main():

    config = load_config()
    

    load_modules(package_list)

    info = {}
    
    info['system'] =get_info()
    info['modules'] = get_info(package_list=package_list)

    info_file_path = f'{config["OUTPUT_PATH"]}/current_system.json'

    write_json(info_file_path, info, verbose=True, pretty=True)
    print (info)



main()