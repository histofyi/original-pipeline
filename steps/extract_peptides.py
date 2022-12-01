from typing import List
from pymol import cmd
import os, sys

import datetime

from functions.cli import load_config
from functions.files import load_facet, write_facet

from Bio.PDB import PDBParser, MMCIFParser, Select
from Bio.PDB.mmcifio import MMCIFIO

class SelectChains(Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        if chain.get_id() == self.chain:
            return 1
        else:          
            return 0


class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0


def extract_peptides(peptide_chains:str, assembly_name:str, input_folder:str, warehouse_path:str):
    solvated_folder = f'{warehouse_path}/structures/coordinates/public/class_i/with_solvent/peptide'
    solvent_free_folder = f'{warehouse_path}/structures/coordinates/public/class_i/without_solvent/peptide'

    peptide_structure_details = {}

    input_filename = f'{input_folder}/{assembly_name}.pdb'
    
    structure = PDBParser(PERMISSIVE=1).get_structure(assembly_name, input_filename)

    coordinate_server_url = 'https://coordinates.histo.fyi/structures/view/class_i/'
    filetype_filename = f'{assembly_name}_peptide.pdb'

    #https://coordinates.histo.fyi/structures/view/class_i/with_solvent/1hhk_1_aligned.pdb


    for chain in structure.get_chains():
        if chain.id in peptide_chains:
            chain_id = chain.id
            print (chain)

            filename = f'{solvated_folder}/{assembly_name}.cif'
            
            io = MMCIFIO()
            io.set_structure(structure)
            io.save(filename, SelectChains(chain_id))

            structure = MMCIFParser(QUIET=True).get_structure(assembly_name, filename)
            
            filename = f'{solvent_free_folder}/{assembly_name}.cif'
            
            io.set_structure(structure)
            io.save(filename, NonHetSelect())

            peptide_structure_details['chain'] = chain_id
            peptide_structure_details['with_solvent'] = f''
            peptide_structure_details['without_solvent'] = f''

    print (input_filename)
    
    


    return {}


def iterator_cleanup():
    cmd.delete('all')



def action_cleanup(pipeline_path:str):
    cif_files = [file_name for file_name in os.listdir(pipeline_path) if '.cif' in file_name]
    for cif_file in cif_files:
        os.remove(f'{pipeline_path}/{cif_file}')
        print (f'Local {cif_file} removed')
    pass


def perform_action(to_process:List, mhc_class:str, input_folder:str, warehouse_path:str, pipeline_path:str):
    """
    Iterates through the list of items to process and performs the action
    """
    step_errors = []
    facet_name = 'peptide_structures'
    for assembly_name in to_process:
        assembly_id = assembly_name.split('_')[1]
        pdb_code = assembly_name.split('_')[0]
        peptide_structures = load_facet(pdb_code, facet_name)
        peptide_details = load_facet(pdb_code, 'peptide_details')
        print ('------')
        print (pdb_code)

        if peptide_details:
            peptide_chains = [chain for chain in peptide_details]
            print (peptide_chains)
            if peptide_structures is None:
                peptide_structures = extract_peptides(peptide_chains, assembly_name, input_folder, warehouse_path)

            else:
                print ('ALREADY SPLIT')
        else:
            print ('NO PEPTIDE')

    return step_errors



def main():
    config = load_config()
    
    input_folder = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/with_solvent/aligned'
    structures = [file.split('.')[0] for file in os.listdir(input_folder)]
    structures = ['1hhk_1','1hhk_2']
    print (structures)
    step_errors = perform_action(structures, 'class_i', input_folder, config['WAREHOUSE_PATH'], config['PIPELINE_PATH'])
    print (step_errors)

    


main()