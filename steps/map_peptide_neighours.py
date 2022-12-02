from typing import List
from pymol import cmd
import os, sys

import datetime

from functions.cli import load_config
from functions.files import load_facet, write_facet

from Bio.PDB import PDBParser, NeighborSearch


def map_peptide_neighbours(alpha_chains:List, peptide_chains:List, assembly_name:str, input_folder:str, warehouse_path:str):
    
    input_filename = f'{input_folder}/{assembly_name}.pdb'
    peptide_neighbours = {}

    alpha_chain = None
    peptide_chain = None
    structure = PDBParser(PERMISSIVE=1).get_structure(assembly_name, input_filename)
    
    all_atoms = []
    for chain in structure.get_chains():
        to_process = False
        if chain.id in alpha_chains:
            alpha_chain = chain.id
            to_process = True  
        elif chain.id in peptide_chains:
            peptide_chain = chain.id
            to_process = True
        if to_process:
            peptide_neighbours[chain.id] = {}
            for residue in chain:
                if residue.id[0] == ' ':
                    for atom in residue:
                        all_atoms.append(atom)
        
    print (alpha_chain)
    print (peptide_chain)
    if alpha_chain and peptide_chain:
        #perform a neighbour search on all alpha and peptide, residue level, cutoff 5A
        neighbours = NeighborSearch(all_atoms).search_all(5, level='R')
        for residue_pair in neighbours:
            residue_1 = residue_pair[0]
            residue_2 = residue_pair[1]
            if residue_1.get_parent().id != residue_2.get_parent().id:
                chain_pair = [residue_1.get_parent().id, residue_2.get_parent().id]
                if peptide_chain in chain_pair and alpha_chain in chain_pair:
                    if residue_1.get_parent().id == alpha_chain:
                        class_i_details = {'residue':residue_1.resname, 'position':residue_1.get_id()[1]}
                        peptide_details = {'residue':residue_2.resname, 'position':residue_2.get_id()[1]}
                    else:
                        peptide_details = {'residue':residue_1.resname, 'position':residue_1.get_id()[1]}
                        class_i_details = {'residue':residue_2.resname, 'position':residue_2.get_id()[1]}


                    #TODO offset the peptide and MHC residue ids if needed
                    
                    class_i_residue_id = class_i_details['position']
                    peptide_residue_id = peptide_details['position']


                    if class_i_residue_id not in peptide_neighbours[alpha_chain]:
                        peptide_neighbours[alpha_chain][class_i_residue_id] = {'position':class_i_residue_id, 'residue':class_i_details['residue'], 'neighbours':[] }
                    peptide_neighbours[alpha_chain][class_i_residue_id]['neighbours'].append(peptide_details)

                    if peptide_residue_id not in peptide_neighbours[peptide_chain]:
                        peptide_neighbours[peptide_chain][peptide_residue_id] = {'position':peptide_residue_id, 'residue':peptide_details['residue'], 'neighbours':[] }
                    if class_i_details not in peptide_neighbours[peptide_chain][peptide_residue_id]['neighbours']:
                        peptide_neighbours[peptide_chain][peptide_residue_id]['neighbours'].append(class_i_details)
    else:
        peptide_neighbours = None
    print (peptide_neighbours)
    return peptide_neighbours


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

        peptide_neighbours = load_facet(pdb_code, 'peptide_neighbours')
        assigned_chains = load_facet(pdb_code, 'assigned_chains')
        
        print (pdb_code)

        if assigned_chains:
            if not peptide_neighbours:
                peptide_neighbours = {}
            class_i_alpha = None
            class_i_alphas = ['class_i_alpha','truncated_class_i_alpha']
            for option in class_i_alphas:
                if option in assigned_chains:
                    class_i_alpha = option
            if class_i_alpha and 'peptide' in assigned_chains:
                alpha_chains = assigned_chains[class_i_alpha]['chains']
                peptide_chains = assigned_chains['peptide']['chains']
                print (alpha_chains)
                print (peptide_chains)
                this_neighbours = map_peptide_neighbours(alpha_chains, peptide_chains, assembly_name, input_folder, warehouse_path)
                if this_neighbours is not None:
                    peptide_neighbours[assembly_id] = this_neighbours
                    write_facet(pdb_code, 'peptide_neighbours', peptide_neighbours)

            else:
                print ('not Class I or no peptide')
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