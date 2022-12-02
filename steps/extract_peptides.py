from typing import List
from pymol import cmd
import os, sys

import datetime

from functions.cli import load_config
from functions.files import load_facet, write_facet

from Bio.PDB import PDBParser, MMCIFParser, Select
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import PDBIO



class SelectChains(Select):
    """ Only accept the specified chains when saving."""
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
    filetype = 'peptide'
    
    solvated_folder = f'{warehouse_path}/structures/coordinates/public/class_i/with_solvent/{filetype}'
    solvent_free_folder = f'{warehouse_path}/structures/coordinates/public/class_i/without_solvent/{filetype}'

    peptide_structure_details = {}

    input_filename = f'{input_folder}/{assembly_name}.pdb'
    
    structure = PDBParser(PERMISSIVE=1).get_structure(assembly_name, input_filename)

    coordinate_server_url = 'https://coordinates.histo.fyi/structures/view/class_i'
    filetype_filename = f'{assembly_name}_{filetype}'


    for chain in structure.get_chains():
        if chain.id in peptide_chains:
            chain_id = chain.id
            print (chain)

            
            mmcif_io = MMCIFIO()
            pdb_io = PDBIO()

            export_formats = {
                'cif':mmcif_io,
                'pdb':pdb_io
            }

            for format in export_formats:
                filename = f'{solvated_folder}/{assembly_name}.{format}'
                io = export_formats[format]
                io.set_structure(structure)
                io.save(filename, SelectChains(chain_id))

            peptide_only_filename = f'{solvated_folder}/{assembly_name}.cif'
            
            try:
                structure = MMCIFParser(QUIET=True).get_structure(assembly_name, peptide_only_filename)
            
                for format in export_formats:
                    filename = f'{solvent_free_folder}/{assembly_name}.{format}'
                    io = export_formats[format]
                    io.set_structure(structure)
                    io.save(filename, NonHetSelect())

                peptide_structure_details = {
                    'chain':chain_id,
                    'formats':{
                        'cif':{},
                        'pdb':{}
                    }
                }
                for format in ['cif','pdb']:
                    for solvent in ['with_solvent', 'without_solvent']:
                        peptide_structure_details['formats'][format][solvent] = f'{coordinate_server_url}/{solvent}/{filetype_filename}.{format}'
            except:
                peptide_structure_details = {}
    if len(peptide_structure_details) == 0:
        peptide_structure_details = None

    return peptide_structure_details


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
        print ('------')
        print (assembly_name)
        
        assembly_id = assembly_name.split('_')[1]
        pdb_code = assembly_name.split('_')[0]

        peptide_structures = load_facet(pdb_code, facet_name)
        peptide_details = load_facet(pdb_code, 'peptide_details')
        
        print (pdb_code)

        if peptide_details:
            peptide_chains = [chain for chain in peptide_details]
            if peptide_structures is None:
                peptide_structures = {}
            if not assembly_id in peptide_structures:
                peptide_assembly_info = extract_peptides(peptide_chains, assembly_name, input_folder, warehouse_path)

                if peptide_assembly_info is not None:
                    peptide_structures[assembly_id] = peptide_assembly_info

                print (peptide_structures)
                if len(peptide_structures) > 0:
                    write_facet(pdb_code, facet_name, peptide_structures)
                else:
                    if assembly_name not in step_errors:
                        step_errors.append(assembly_name)
                
            else:
                print ('ALREADY SPLIT')
        else:
            print ('NO PEPTIDE')

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