import os
import toml

from functions.files import read_json

import datetime 


def update_local_pdb(localpdb_path:str, warehouse_path:str):
    """
    This function updates the localpdb instance and any installed plugins

    Args:
        localpdb_path(str): the path to the localpdb instance
        warehouse_path(str): the path to the local copy of the histo.fyi warehouse
    """
    version_info = read_json(f'{warehouse_path}/current_version.json')
    current_version = version_info['version']


    dateformat = '%Y%m%d'

    current_version_date = datetime.datetime.strptime(current_version, dateformat)
    next_version_date = current_version_date + datetime.timedelta(days=7)
    next_version = next_version_date.strftime(dateformat)

    version_path = f'{localpdb_path}/data/{next_version}'
    version_file = f'{version_path}/obsolete.txt'

    if not os.path.exists(version_file):
        if not os.path.exists(version_path):
            os.makedirs(version_path)
        version_command = f'touch {version_file}'
        os.system(version_command)
        print (f'Creating version file for {next_version}')

    update_command = f'localpdb_setup -db_path {localpdb_path} --update'

    os.system(update_command)



def main():
    config = toml.load('config.toml')
    update_local_pdb(config['LOCALPDB_PATH'], config['WAREHOUSE_PATH'])


main()
