from typing import List
from pymol import cmd
import os, sys

import datetime

from functions.cli import load_config
from functions.files import load_facet, write_facet, load_constants

from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import one_to_three

pocket_definitions = load_constants('pockets')['class_i']

def map_pockets(alpha_chains:List, species:str, assembly_name:str, input_folder:str, warehouse_path:str):
    direct_mapping = ['homo_sapiens','macaca_mulatta','equus_caballus','felis_catus','ailuropoda_melanoleuca','bos_taurus','oryctolagus_cuniculus','sus_scrofa','mus_musculus','rattus_norvegicus']
    
    print (pocket_definitions)
    if species in direct_mapping:
        pockets = {}
        for pocket_letter in pocket_definitions:
            pockets[pocket_letter] = {}
            for position in pocket_definitions[pocket_letter]:
                number = int(position.replace('a',''))
                index = number - 1
                one_letter = alpha_chains['sequence'][index] 
                pockets[pocket_letter][position] = {
                    'number':number,
                    'one_letter':one_letter,
                    'three_letter':one_to_three(one_letter)
                }
        print (pockets)
    else:
        pockets = None
        print (f'{species} needs special handling')
    return pockets


def iterator_cleanup():
    pass



def action_cleanup(pipeline_path:str):
    pass


def perform_action(to_process:List, mhc_class:str, input_folder:str, warehouse_path:str, pipeline_path:str):
    """
    Iterates through the list of items to process and performs the action
    """
    step_errors = []
    for assembly_name in to_process:
        print ('------')
        print (assembly_name)
        
        assembly_id = assembly_name.split('_')[1]
        pdb_code = assembly_name.split('_')[0]

        facet_name = 'pockets'

        pockets = load_facet(pdb_code, facet_name)
        assigned_chains = load_facet(pdb_code, 'assigned_chains')
        species = load_facet(pdb_code, 'species')
        print (pdb_code)

        if assigned_chains:
            if not pockets:
                pockets = {}
            class_i_alpha = None
            class_i_alphas = ['class_i_alpha','truncated_class_i_alpha']
            for option in class_i_alphas:
                if option in assigned_chains:
                    class_i_alpha = option
            if class_i_alpha:
                alpha_chains = assigned_chains[class_i_alpha]
                print (species)
                this_pockets = map_pockets(alpha_chains, species['slug'], assembly_name, input_folder, warehouse_path)

                if this_pockets:
                    pockets[assembly_id] = this_pockets
                    write_facet(pdb_code, facet_name, pockets)
            else:
                print ('not classical Class I')
                if assembly_name not in step_errors:
                    step_errors.append(assembly_name)
    return step_errors



def main():
    config = load_config()
    
    input_folder = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/without_solvent/aligned'
    structures = [file.split('.')[0] for file in os.listdir(input_folder) if len(file.split('.')[0]) > 0]
    #structures = ['1hhk_1','1hhk_2']
    print (structures)
    step_errors = perform_action(structures, 'class_i', input_folder, config['WAREHOUSE_PATH'], config['PIPELINE_PATH'])
    print (step_errors)

    


main()