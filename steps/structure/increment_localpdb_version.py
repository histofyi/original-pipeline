from typing import Dict, List, Union, Optional

import datetime
import os

from shared.pipeline import load_config
from shared.files import read_json


def increment_localpdb_version(config:Dict) -> Union[List, List, int]:

    warehouse_path = config['WAREHOUSE_PATH']
    localpdb_path = config['LOCALPDB_PATH']

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
    else:
        print (f'Version file for {next_version} already created')
    pass

    
if __name__ == "__main__":
    config = load_config()
    increment_localpdb_version(config)