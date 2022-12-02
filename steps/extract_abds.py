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


class SelectResidues(Select):
    """ Only accept the specified residues when saving. """
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def accept_residue(self, res):
        if res.id[1] >= self.start and res.id[1] <= self.end:
            return True
        else:
            return False


def extract_abds(alpha_chains:str, assembly_name:str, input_folder:str, warehouse_path:str):
    
    alpha_chain_folder = f'{warehouse_path}/structures/coordinates/public/class_i/without_solvent/alpha_chains'
    abd_only_folder = f'{warehouse_path}/structures/coordinates/public/class_i/without_solvent/antigen_binding_domains'

    alpha_chain_structure_details = {}
    abd_only_structure_details = {}

    input_filename = f'{input_folder}/{assembly_name}.pdb'
    
    structure = PDBParser(PERMISSIVE=1).get_structure(assembly_name, input_filename)

    coordinate_server_url = 'https://coordinates.histo.fyi/structures/view/class_i'


    for chain in structure.get_chains():
        if chain.id in alpha_chains:
            chain_id = chain.id
            print (chain)

            mmcif_io = MMCIFIO()
            pdb_io = PDBIO()

            export_formats = {
                'cif':mmcif_io,
                'pdb':pdb_io
            }

            for format in export_formats:
                filename = f'{alpha_chain_folder}/{assembly_name}.{format}'
                io = export_formats[format]
                io.set_structure(structure)
                io.save(filename, SelectChains(chain_id))

            alpha_only_filename = f'{alpha_chain_folder}/{assembly_name}.cif'
            
            try:
                structure = MMCIFParser(QUIET=True).get_structure(assembly_name, alpha_only_filename)
            
                for format in export_formats:
                    filename = f'{abd_only_folder}/{assembly_name}.{format}'
                    io = export_formats[format]
                    io.set_structure(structure)
                    io.save(filename, SelectResidues(1,181))

                alpha_chain_structure_details = {
                    'formats':{
                        'cif':{},
                        'pdb':{}
                    }
                }
                abd_only_structure_details = {
                    'formats':{
                        'cif':{},
                        'pdb':{}
                    }
                }
                for format in ['cif','pdb']:
                    alpha_chain_structure_details['formats'][format]['without_solvent'] = f'{coordinate_server_url}/without_solvent/{assembly_name}_alpha_chain.{format}'
                    abd_only_structure_details['formats'][format]['without_solvent'] = f'{coordinate_server_url}/without_solvent/{assembly_name}_abd.{format}'
            except:
                alpha_chain_structure_details = {}
                abd_only_structure_details = {}
    if len(alpha_chain_structure_details) == 0:
        alpha_chain_structure_details = None
        abd_only_structure_details = None

    return alpha_chain_structure_details, abd_only_structure_details


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
    for assembly_name in to_process:
        print ('------')
        print (assembly_name)
        
        assembly_id = assembly_name.split('_')[1]
        pdb_code = assembly_name.split('_')[0]

        alpha_chain_structure_details = load_facet(pdb_code, 'alpha_chain_structures')
        abd_only_structure_details = load_facet(pdb_code, 'antigen_binding_domain_structures')
        assigned_chains = load_facet(pdb_code, 'assigned_chains')
        
        print (pdb_code)

        if assigned_chains:
            if not alpha_chain_structure_details:
                alpha_chain_structure_details = {}
            if not abd_only_structure_details:
                abd_only_structure_details = {}
            class_i_alpha = None
            class_i_alphas = ['class_i_alpha','truncated_class_i_alpha','mr1','cd1a','cd1b','cd1c','cd1d','fcrn','mica','micb']
            for option in class_i_alphas:
                if option in assigned_chains:
                    class_i_alpha = option
            if class_i_alpha:
                alpha_chains = assigned_chains[class_i_alpha]['chains']
                print (alpha_chains)
                this_alpha_chain, this_abd = extract_abds(alpha_chains, assembly_name, input_folder, warehouse_path)
                if this_alpha_chain and this_abd:
                    alpha_chain_structure_details[assembly_id] = this_alpha_chain
                    abd_only_structure_details[assembly_id] = this_abd

                    write_facet(pdb_code, 'alpha_chain_structures', alpha_chain_structure_details)
                    write_facet(pdb_code, 'antigen_binding_domain_structures', abd_only_structure_details)
            else:
                print ('not Class I')
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