from typing import List
from pymol import cmd
import os

import datetime


folder = 'edited'


def align_to_canonical(mhc_class:str, pdb_code:str, assembly_id:int):
    print ('-----')
    print (pdb_code)
    cmd.load(f'input/{folder}/{pdb_code}_{assembly_id}.cif', quiet=1)
    cmd.save(f'output/{mhc_class}/with_solvent/raw/{pdb_code}_{assembly_id}.pdb')

    cmd.load('canonical/cannonical_class_i_1hhk.pdb', quiet=1)
    

    align = cmd.cealign('cannonical_class_i_1hhk',pdb_code)
    align['rmsd'] = align['RMSD']
    del align['RMSD']
    cmd.delete('cannonical_class_i_1hhk')
    print (align['rmsd'])
    file_types = ['pdb', 'cif']
    for file_type in file_types:
        file_name = f'output/{mhc_class}/with_solvent/aligned/{pdb_code}_{assembly_id}.{file_type}'
        cmd.save(file_name)
    return align


def iterator_cleanup():
    cmd.delete('all')



def action_cleanup():
    pass


def perform_action(to_process:List, mhc_class:str):
    """
    Iterates through the list of items to process and performs the action
    """
    step_errors = []
    for assembly_name in to_process:
        try:
            assembly_id = assembly_name.split('_')[1]
            pdb_code = assembly_name.split('_')[0]
        except:
            assembly_id = None
            pdb_code = None
            item_errors.append({'assembly_name':assembly_name, 'error':'split_error'})
        try:
            alignment_dict = {}
            item_errors = []
            if pdb_code:
                alignment = align_to_canonical(mhc_class, pdb_code, assembly_id)
                aligned = {
                    'aligned_on': mhc_class,
                    'alignment': alignment,
                    'last_updated': datetime.datetime.now().isoformat()
                }
                alignment_dict[assembly_id] = aligned
        except:
            item_errors.append({'assembly_name':assembly_name, 'error':'alignment_error'})
        iterator_cleanup()
        print (alignment_dict)

        # write aligned facet

        if len(item_errors) > 0:
            for error in item_errors:
                step_errors.append(error)
    action_cleanup()
    return step_errors



def main():
    structures = [file.split('.')[0] for file in os.listdir(f'input/{folder}')]
    #pdb_codes = ['1hhk_1', '1hhk_2']
    print (structures)
    step_errors = perform_action(structures, 'class_i')
    print (step_errors)


main()