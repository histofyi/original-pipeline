from typing import Dict, List, Tuple

import toml
import json

from functions.cli import print_spacer
from functions.files import write_json, read_json
from functions.pdb import load_structure

from localpdb import PDB

from Bio.PDB.Polypeptide import three_to_one, is_aa

from rich.progress import Progress
from rich.console import Console

config = toml.load('config.toml')
console = Console()


localpdb_path = config['LOCALPDB_PATH']


set_list = [
    'matching_peptide',
    'matching_natural_peptide',
    'matching_unnatural_peptide',
    'modified_residues',
    'disordered',
    'missing_residues',
    'different_sequences',
    'missing_start',
    'missing_end',
    'incorrect_start',
    'structure_errors',
    'short_peptides',
    'normal_length_peptides',
    'no_peptide',
    'hetatoms',
    'has_peptide'
]

sets = {}

for item in set_list:
    sets[item] = []





def deduplicate_list(pdb_code_list:List):
    deduplicated_list = [pdb_code for pdb_code in set(pdb_code_list)]
    sorted_list = sorted(deduplicated_list)
    return sorted_list


def sequence_is_peptide(sequence:str, chain_id) -> Dict:
    """
    This function checks if the sequence is a peptide or not (based on length)
    It also looks for the presence of unnatural amino acids

    Args:
        sequence (str): the sequence of the chain e.g. SIINFEKL
        chain_id (str): the chain_id of the chain e.g. A


    Returns:
        Dict: returns a dictionary of information about the peptide, or None if it's not likely to be a peptide
    """
    chain_peptide_info = None
    if len(sequence) < 25:
        if 'X' in sequence:
            contains_unnatural = True
        else:
            contains_unnatural = False
        chain_peptide_info = {
            'chain_id': chain_id,
            'full_sequence': sequence,
            'full_length': len(sequence),
            'actual_sequence': None,
            'actual_length': None,
            'features': []
        }
        if contains_unnatural:
            chain_peptide_info['features'].append('contains_unnatural')
        return chain_peptide_info
    else:
        return None


def residue_checks(pdb_code, residue):
    is_disordered = residue_is_disordered(pdb_code, residue)
    is_amino_acid, is_natural, residue_name = residue_is_natural(pdb_code, residue)
    return is_amino_acid, is_natural, residue_name    


def residue_is_disordered(pdb_code:str, residue):
    if residue.is_disordered():
        sets['disordered'].append((pdb_code, residue.id[1], residue.resname))
        return True
    else:
        return False


def residue_is_natural(pdb_code:str, residue) -> Tuple[bool, bool, str]:
    residue_name = residue.resname
    if is_aa(residue_name):
        if is_aa(residue_name, standard=True):
            return True, True, residue_name
        else:
            sets['modified_residues'].append((pdb_code, residue.id[1], residue.resname))
            return True, False, residue_name
    else:
        if residue_name not in sets['hetatoms']:
            sets['hetatoms'].append(residue_name)
        return False, False, None


def load_pdb_lists(mhc_class:str, warehouse_path:str) -> List:
    console.rule(f'[bold]Loading pdb code lists for - {mhc_class}')
    pdb_codes = []
    for structure_type in ['class_i','truncated_class_i']:
        print (f'Loading {structure_type}')
        filepath = f'{warehouse_path}/queries/{structure_type}_hits.json'
        json_data = read_json(filepath)
        pdb_codes += [pdb_code for pdb_code in json_data]
    return pdb_codes






pdb_codes = load_pdb_lists('class_i', config['WAREHOUSE_PATH'])


console.rule(f'[bold]Matching pdb_codes with localpdb')

lpdb = PDB(db_path=localpdb_path)


all_chains = lpdb.chains
all_structures_count = len(pdb_codes)



with Progress() as progress:
        
    task = progress.add_task("[white]Processing...", total=all_structures_count)

    for pdb_code in pdb_codes:
        # create a list of all the chains in a structure
        structure_chains = [chain for chain in all_chains.loc[all_chains['pdb'] == pdb_code].index]
        # initialise the peptide chains dictionary
        peptide_chains = {}
        
        print (pdb_code)
        progress.update(task, advance=1)

        # iterate through the chains
        for chain in structure_chains:
            sequence = all_chains.loc[chain]['sequence']
            chain_id = chain.split('_')[1]
            chain_peptide_info = sequence_is_peptide(sequence, chain_id)
            if chain_peptide_info:
                peptide_chains[chain_id] = chain_peptide_info
        
        # if there are peptide chains detected, then see if the peptide within the structure matches the sequence
        # places where it may not match:
        # - disorder in a residue
        # - missing residues at the N- or C-terminus
        # - missing residues in the centre for longer peptides
        # - incorrect sequence in PDB file
        # also check for short peptides (peptide fragments)
        
        if len(peptide_chains) > 0:
            # load the structure from the local pdb copy
            structure = load_structure(pdb_code, lpdb.entries.loc[pdb_code, 'mmCIF_fn'])
            if structure:    
                # iterate through the chains in the structure
                for chain in structure.get_chains():
                    # operate only on the chains we've marked as peptides 
                    chain_id = chain.id
                    if chain_id in peptide_chains:
                        structure_sequence = ''
                        peptide_positions = []
                        # iterate through the residues in the chain
                        for residue in chain:
                            residue_id = residue.id[1]
                            # first check if it's an amino acid and if it's a natural one
                            is_amino_acid, is_natural, residue_name = residue_checks(pdb_code, residue)
                            if is_amino_acid:
                                # then append the single letter code to the sequence of the actual structure
                                if is_natural:
                                    structure_sequence += three_to_one(residue_name)
                                else:
                                    structure_sequence += 'X'
                                # and append the position, this may not be the position in the peptide, some are numbered as position in whole antigen
                                peptide_positions.append(residue_id)

                        # add the structure derived sequence and length to the peptide chains dictionary            
                        peptide_chains[chain_id]['actual_sequence'] = structure_sequence
                        peptide_chains[chain_id]['actual_length'] = len(structure_sequence)

                        if structure_sequence == peptide_chains[chain_id]['full_sequence']:
                            # the sequence in the pdb file and the sequence of the structure match
                            sets['matching_peptide'].append(pdb_code)
                            peptide_chains[chain_id]['features'].append('correct_sequence_and_length')
                            # check whether all amino acids are natural
                            if 'X' not in structure_sequence:
                                sets['matching_natural_peptide'].append(pdb_code)
                            else:
                                sets['matching_unnatural_peptide'].append(pdb_code)
                                peptide_chains[chain_id]['features'].append('unnatural_amino_acids')
                        else:
                            # the sequence in the pdb file and the sequence of the structure don't match
                            # first of all we'll test if the structure sequence is contained within the pdb sequence
                            if structure_sequence in peptide_chains[chain_id]['full_sequence']:
                                structure_sequence_index = peptide_chains[chain_id]['full_sequence'].index(structure_sequence)
                                # we'll test to see if the structure sequence is at the start of the pdb sequence
                                if structure_sequence_index == 0:
                                    # it's most likely a C-terminal extension that's disordered
                                    sets['missing_end'].append(pdb_code)
                                    peptide_chains[chain_id]['features'].append('missing_end')
                                else:
                                    # it's most likely an N-terminal section that's disordered
                                    sets['missing_start'].append(pdb_code)
                                    peptide_chains[chain_id]['features'].append('missing_start')
                            else:
                                # if the structure sequence is not contained in the pdb sequence it is because of mismatch or disordered/missing residues
                                missing_residue_count = peptide_chains[chain_id]['full_length'] - peptide_chains[chain_id]['actual_length']
                                if missing_residue_count == 0:
                                    # it's a sequence mismatch
                                    sets['different_sequences'].append((pdb_code, missing_residue_count, peptide_chains[chain_id]['full_sequence'], peptide_chains[chain_id]['actual_sequence']))
                                    peptide_chains[chain_id]['features'].append('sequence_mismatch')
                                else:
                                    # it's a missing chunk of peptide 
                                    sets['missing_residues'].append((pdb_code, missing_residue_count, peptide_chains[chain_id]['full_sequence'], peptide_chains[chain_id]['actual_sequence']))
                                    previous_position = 0
                                    previous_aa = ''
                                    missing_residues = ''
                                    pre_gap = ''
                                    post_gap = ''
                                    gap_start = 0
                                    gap_end = 0
                                    gap_length = 0
                                    has_gap = False
                                    # iterate through the chain
                                    for residue in chain:
                                        current_id = residue.id[1]
                                        residue_name = residue.resname
                                        # only look for gaps involving amino acids, both natural and unnatural
                                        if is_aa(residue_name):
                                            try:
                                                current_aa = three_to_one(residue_name)
                                            except:
                                                current_aa = 'X'
                                            # check if there is a gap. incrememnt of greater than one indicates a gap
                                            if current_id - previous_position != 1:
                                                has_gap = True
                                                peptide_chains[chain_id]['gap_start_after'] = previous_position
                                                peptide_chains[chain_id]['gap_ends_at'] = current_id
                                                gap_length = current_id - 1 - previous_position
                                                peptide_chains[chain_id]['gap_length'] = gap_length
                                                peptide_chains[chain_id]['gap_sequence'] = peptide_chains[chain_id]['full_sequence'][previous_position: current_id - 1]
                                                print (f'GAP : P{previous_position}{previous_aa} to P{current_id}{current_aa}')
                                            # build pre_gap and post_gap strings, these will be combined with multiple dashes to form the gapped sequence
                                            if not has_gap:
                                                pre_gap += current_aa
                                            else:
                                                post_gap += current_aa
                                            previous_position = current_id
                                            previous_aa = current_aa
                                    # at the end of the loop, append the gapped sequence and add a feature
                                    peptide_chains[chain_id]['gapped_sequence'] = f'{pre_gap}{"-"*gap_length}{post_gap}'
                                    peptide_chains[chain_id]['features'].append('central_gap')
                        # next we'll test whether the peptide is correctly numbered
                        if peptide_positions[0] != 1:
                            print (f'Chain {chain_id} starts at {peptide_positions[0]}')
                            sets['incorrect_start'].append((pdb_code, peptide_positions[0]))
                            peptide_chains[chain_id]['features'].append('incorrect_start')
                        # finally we'll test if it's a short peptide
                        if len(structure_sequence)  < 8:
                            sets['short_peptides'].append((pdb_code, len(structure_sequence), structure_sequence))
                            peptide_chains[chain_id]['features'].append('short_peptide')
                        else:
                            sets['normal_length_peptides'].append(pdb_code)
            else:
                sets['structure_errors'].append(pdb_code)
        else:
            # or there is no peptide present, as in some non-classical and some early structures
            sets['no_peptide'].append(pdb_code)
            peptide_chains = None

        if peptide_chains:
            sets['has_peptide'].append(pdb_code)
            print (peptide_chains)
            filepath = f'{config["WAREHOUSE_PATH"]}/structures/info/public/peptide_details/{pdb_code}.json'
            write_json(filepath, peptide_chains, verbose=True, pretty=True)
        print_spacer()




console.rule('[bold]Writing files and tidying up')

for item in set_list:
    
    filepath = f'tmp/{item}.json'
    print (item)
    
    sets[item] = deduplicate_list(sets[item])
    print (f'{len(sets[item])} items')
    write_json(filepath, sets[item], verbose=True, pretty=True)
    print_spacer()



print (f'{all_structures_count} processed')

console.print('[bold green]Done')
