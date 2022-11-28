from typing import List
from pymol import cmd
import os, sys

import datetime

from functions.cli import load_config
from functions.files import load_facet, write_facet





def align_to_canonical(mhc_class:str, assembly_name:str, input_folder:str, output_folder:str, warehouse_path:str):
    # we get CIF files from pdbe
    cmd.load(f'{input_folder}/{assembly_name}.cif', quiet=1)
    # but we also want PDB files, so save a copy now
    cmd.save(f'{input_folder}/{assembly_name}.pdb')

    cmd.load(f'{warehouse_path}/structures/canonical/cannonical_class_i_1hhk.pdb', quiet=1)
    

    align = cmd.cealign('cannonical_class_i_1hhk',assembly_name)
    align['rmsd'] = align['RMSD']
    del align['RMSD']
    cmd.delete('cannonical_class_i_1hhk')
    print (align['rmsd'])
    file_types = ['pdb', 'cif']
    for file_type in file_types:
        file_name = f'{output_folder}/{assembly_name}.{file_type}'
        cmd.save(file_name)
    return align


def iterator_cleanup():
    cmd.delete('all')



def action_cleanup(pipeline_path:str):
    cif_files = [file_name for file_name in os.listdir(pipeline_path) if '.cif' in file_name]
    for cif_file in cif_files:
        os.remove(f'{pipeline_path}/{cif_file}')
        print (f'Local {cif_file} removed')
    pass


def perform_action(to_process:List, mhc_class:str, input_folder:str, output_folder:str, warehouse_path:str, pipeline_path:str):
    """
    Iterates through the list of items to process and performs the action
    """
    step_errors = []
    facet_name = 'alignment'
    for assembly_name in to_process:
        assembly_id = assembly_name.split('_')[1]
        pdb_code = assembly_name.split('_')[0]
        alignment_dict = load_facet(pdb_code, facet_name)
        if alignment_dict is None:
            alignment_dict = {}

        print ('-----')
        print (pdb_code)
        item_errors = []
        if not assembly_id in alignment_dict:
            if pdb_code:
                alignment = align_to_canonical(mhc_class, assembly_name, input_folder, output_folder, warehouse_path)
                aligned = {
                    'aligned_on': mhc_class,
                    'alignment': alignment,
                    'last_updated': datetime.datetime.now().isoformat()
                }
                alignment_dict[assembly_id] = aligned
            iterator_cleanup()
            print (alignment_dict)

            # write aligned facet
            write_facet(pdb_code, facet_name, alignment_dict)



        if len(item_errors) > 0:
            for error in item_errors:
                step_errors.append(error)
    action_cleanup(pipeline_path)
    return step_errors



def main():
    config = load_config()
    input_folders = ['raw','edited']
    for folder in input_folders:
        input_folder = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/with_solvent/{folder}'
        output_folder = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/with_solvent/aligned'
        structures = [file.split('.')[0] for file in os.listdir(input_folder)]
        print (structures)
        step_errors = perform_action(structures, 'class_i', input_folder, output_folder, config['WAREHOUSE_PATH'], config['PIPELINE_PATH'])
        print (step_errors)

    


main()