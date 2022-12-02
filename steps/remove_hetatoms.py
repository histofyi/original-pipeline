from typing import List
from pymol import cmd
import os, sys

import datetime

from functions.cli import load_config
from functions.files import load_facet, write_facet

from Bio.PDB import PDBParser, MMCIFParser, Select
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import PDBIO


class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0


def remove_hetatoms(assembly_name:str, input_folder:str, warehouse_path:str):
    
    aligned_folder = f'{warehouse_path}/structures/coordinates/public/class_i/without_solvent/aligned'

    input_filename = f'{input_folder}/{assembly_name}.pdb'
    
    structure = PDBParser(PERMISSIVE=1).get_structure(assembly_name, input_filename)

    mmcif_io = MMCIFIO()
    pdb_io = PDBIO()

    export_formats = {
        'cif':mmcif_io,
        'pdb':pdb_io
    }

    for format in export_formats:
        filename = f'{aligned_folder}/{assembly_name}.{format}'
        io = export_formats[format]
        io.set_structure(structure)
        io.save(filename, NonHetSelect())

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
    facet_name = 'antigen_binding_domain_structures'
    for assembly_name in to_process:
        print ('------')
        print (assembly_name)
        
        assembly_id = assembly_name.split('_')[1]
        pdb_code = assembly_name.split('_')[0]

        remove_hetatoms(assembly_name, input_folder, warehouse_path)

    return step_errors



def main():
    config = load_config()
    
    input_folder = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/with_solvent/aligned'
    structures = [file.split('.')[0] for file in os.listdir(input_folder) if len(file.split('.')[0]) > 0]
    #structures = ['1hhk_1','1hhk_2']
    print (structures)
    step_errors = perform_action(structures, 'class_i', input_folder, config['WAREHOUSE_PATH'], config['PIPELINE_PATH'])
    print (step_errors)

    


main()