import os
import toml


def update_local_pdb(localpdb_path:str):
    """
    This function updates the localpdb instance and any installed plugins

    Args:
        localpdb_path(str): the path to the localpdb instance
    """
    update_command = f'localpdb_setup -db_path {localpdb_path} --update'
    os.system(update_command)



def main():
    config = toml.load('config.toml')
    update_local_pdb(config['LOCALPDB_PATH'])


main()
