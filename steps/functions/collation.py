from typing import List, Dict, Tuple



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
        'title': None, # the constructed title sentance for the structure
        'assemblies':{}, # the assemblies (copies of the complex) in the structure
        'assigned_chains':{}, # the chain assignements for the complex
        'species':{}, # information on the organism 
        'class':None, # MHC class e.g. class_i
        'classical': None, # whether the MHC molecule is classical (e.g. HLA-A) or non-classical (e.g. CD1a)
        'complex_type':None, # a slug for the type of complex e.g. class_i_with_peptide_and_alpha_beta_tcr
        'locus':None, # the locus e.g. HLA-A
        'allele':{
            'alpha':{}, # the allele information of the alpha chain e.g. HLA-A*68:01
            'beta':{} # the allele name of the beta chain e.g humanB2m, HLA-DRB*04:01
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
        'tcr':None, #details of the TCR in the structure (from STCRDab)
        'ligands':[], # other ligands, or the ligands of non-peptide binding molecules e.g. CD1a binding lipids
        'components':[], # a list of the components in the structure
        'resolution':None, # the resolution of the structure e.g 1.9
        'experimental_method':None, # the experimental method for solution of the structure e.g. diffraction, cryoem
        'assembly_count':None, # how many assemblies there are
        'unique_chain_count': None, # how many unique chains there are
        'chronology':{
            'deposition_date':None, # when the structure was deposited
            'deposition_year':None, # the year the structure was deposited
            'release_date':None, # when the structure was released
            'release_year':None, # the year the structure was deposited
            'revision_date':None, # when the structure was revised
            'revision_year':None # the year the structure was deposited
        },
        'missing_residues':[], # a list of missing residues e.g. {"chain":"A","position":40,"residue_name":"ALA"}
        'pdb_title':None, # the title of the structure (from the PDB)
        'publication':{}, # details on the publication
        'manually_edited': {}, # whether the information has been manually edited or corrected
        'facets':{}, # which facets are present for this structure e.g. peptide_neighbours, peptide_backbone
        'same_as':{} # other representations of this struture
    }





def build_list_data(pdb_code:str) -> Dict:
    """
    This function returns the listing metadata dictionary which is filled out by the different methods on the structure pipeline

    Args:
        pdb_code (str): the pdb code of the structure

    Returns:
        Dict: the default prototype data dictionary for structure listing metadata
    """
    return {
        'pdb_code':pdb_code,
        'title': None,
        'species':None,
        'deposition_year':None,
        'components':[],
        'pdb_title':None
    }


