from typing import Dict, List
import os

from functions.cli import load_config, parse_args
from functions.files import write_json

from rich.progress import Progress
from rich.console import Console



def build_core_data(pdb_code: str) -> Dict:
    """
    This function returns the core metadata dictionary which is filled out by the different methods on the structure pipeline

    Args:
        pdb_code (str): the pdb code of the structure

    Returns:
        Dict: the default prototype data dictionary for structure metadata
    """
    return {
        'pdb_code':pdb_code, # the pdb code
        'assemblies':[], # the assemblies (copies of the complex) in the structure
        'organism':{}, # information on the organism 
        'class':None, # MHC class e.g. class_i
        'classical': None, # whether the MHC molecule is classical (e.g. HLA-A) or non-classical (e.g. CD1a)
        'complex_type':None, # a slug for the type of complex e.g. class_i_with_peptide_and_alpha_beta_tcr
        'locus':None, # the locus e.g. HLA-A
        'allele':{
            'alpha':None, # the allele name of the alpha chain e.g. HLA-A*68:01
            'beta':None # the allele name of the beta chain e.g humanB2m, HLA-DRB*04:01
        },
        'peptide':{
            'full_sequence':None, # the sequence of the peptide as specified in the pdb_file
            'actual_sequence':None, # the sequence of the peptide in the structure co-ordinates
            'gapped_sequence': None, # the sequence including gaps for missing residues in the structure
            'epitope_info':{}, # the organism, protein, from and to information about the peptides, and mutations
            'gap_info':{}, # information about any gaps in the peptide
            'length':{
                'numeric':None, # the peptide length e.g. 9
                'text':None # the peptide length name e.g nonamer
            },
            'unnatural_amino_acids':[], # information on any unnatural amino acids in the peptide
            'features':[] # features of the peptide e.g.
        },
        'ligands':[], # other ligands, or the ligands of non-peptide binding molecules e.g. CD1a binding lipids
        'accessory_molecules':{}, # a dictionary of accessory molecules and receptors
        'resolution':None, # the resolution of the structure e.g 1.9
        'methodology':None, # the methodology of the structure e.g. diffraction, cryoem
        'assembly_count':None, # how many assemblies there are
        'chain_count':None, # how many chains there are
        'unique_chain_count': None, # how many unique chains there are
        'chronology':{
            'deposition_date':None, # when the structure was deposited
            'release_date':None, # when the structure was released
            'update_date':None # when the structure was updated
        },
        'missing_residues':[], # a list of missing residues e.g. {"chain":"A","position":40,"residue_name":"ALA"}
        'pdb_title':None, # the title of the structure (from the PDB)
        'publication':{}, # details on the publication
        'doi':{'doi':None, 'url':None}, # the DOI and resolved URL for the publication
        'open_access':False, # if the paper is available under Open Access
        'manually_edited': {}, # whether the information has been manually edited or corrected
        'facets':{} # which facets are present for this structure e.g. peptide_neighbours, peptide_backbone
    }


def build_core_for_list(pdb_codes:str, warehouse_path:str, console, force:bool=False):
    """
    This function will build the core information for a set of pdb codes

    Args:
        pdb_codes (str): the pdb codes to be acted upon
        warehouse_path (str): the filepath to the local copy of the warehouse
        force (bool): whether the existing records should be wiped (use with care!)

    """
    record_count = len(pdb_codes)
    records_changed = []
    for pdb_code in pdb_codes:
        filepath = f'{warehouse_path}/structures/info/public/core/{pdb_code}.json'
        # if we're forcing the update of this file (and removing all information from it)
        if force:
            core_info = build_core_data(pdb_code)
        else:
            # check if file exists, if it does, don't do anything
            if not os.path.exists(filepath):
                core_info = build_core_data(pdb_code)
            else:
                core_info = None
        if core_info:
            records_changed.append(pdb_code)
            write_json(filepath, core_info, verbose=True, pretty=True)
    print (f'{len(records_changed)} out of {record_count} have been changed')
    

        

def main():

    config = load_config()

    console.rule('[bold] Initialiseing records')
    pdb_codes, mhc_class, set_slug, set_context, force = parse_args(console)

    build_core_for_list(pdb_codes, config['WAREHOUSE_PATH'], console, force=force)


console = Console()
main()

    




