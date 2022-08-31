import os
import toml
import json


def get_localpdb_version(localpdb_path:str, output_path:str):
    """
    This function updates the localpdb version file

    Args:
        localpdb_path(str): the path to the localpdb instance
        output_path(str): the path to the output folder
    """

    status_log_path = f'{localpdb_path}/data/status.log'
    with open(status_log_path, "r") as infile:
        version_info = json.load(infile)
    for version in version_info:
        if version_info[version][0] == 'OK':
            current_version = {
                'version':version,
                'last_updated':version_info[version][1]
            }
    current_version_file_path = f'{output_path}/current_version.json'
    with open(current_version_file_path, "w") as outfile:
        outfile.write(json.dumps(current_version, sort_keys=True, indent=4))

    


def main():
    config = toml.load('config.toml')
    get_localpdb_version(config['LOCALPDB_PATH'],config['OUTPUT_PATH'])


main()
