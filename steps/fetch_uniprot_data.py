from typing import Dict, List, Tuple


from common.providers import rcsbProvider, PDBeProvider, httpProvider
from Bio import SeqIO
from io import StringIO

from fuzzywuzzy import fuzz


from functions.pdb import load_pdb_lists
from functions.cli import load_config, print_spacer
from functions.files import load_constants, load_facet, write_facet
from functions.helpers import slugify

from rich.progress import Progress
from rich.console import Console
from rich import print_json

console = Console()
config = load_config()

completed = []
failed = []

missing_chain_assignments = []

mhc_class = 'class_i'

pdb_codes = load_pdb_lists('class_i', config['WAREHOUSE_PATH'], console)

pdb_codes = ['1hhk','6eny','6mpp', '6o9c','6o9b']

http = httpProvider()


def fasta_reader(fasta_string:str):
    fasta_io = StringIO(fasta_string) 
    records = SeqIO.parse(fasta_io, "fasta")
    record = [record for record in records][0]
    fasta_io.close() 
    return record


def trim_sequence(sequence:str):
    if sequence[0] == 'M':
        sequence = sequence[1:]
    return sequence


def match_peptide(peptide:str, fasta_record):
    match = {}
    match['structure_sequence'] = peptide
    if peptide in fasta_record.seq:
        # we have an exact match

        sequence_starts = fasta_record.seq.index(peptide)
        sequence_ends = sequence_starts + len(peptide)

        match['match_type'] = 'exact'

        
        print (f'MATCH STARTS AT {sequence_starts}')
        print ('PERFECT MATCH\n')
    else:
        # no exact match, so we're going to run a 3 amino acid long window along the peptide sequence and see where it first matches
        window_length = 3
        peptide_chops = [position for position in range(1, len(peptide) - window_length + 1)]
        match_start = 0
        match_fragment = None
        print (peptide_chops)
        for chop in peptide_chops:
            if len(peptide) >= chop + 3:
                chopped = peptide[chop - 1:chop + window_length - 1]
                if match_start == 0:
                    if chopped in fasta_record.seq:
                        print (chop)
                        print(f'PARTIAL MATCH starting with {chopped} at P{chop}')
                        match_start = chop
                        match_fragment = chopped
        if not match_start:
            match['type'] = 'mismatch'
        else:
            match['type'] = 'partial'
            differences = {}
            match_offset = match_start - 1
            sequence_starts = fasta_record.seq.index(match_fragment) - match_offset
            sequence_ends = sequence_starts + len(peptide)
            uniprot_sequence = fasta_record.seq[sequence_starts: sequence_ends]
            changes = 0
            j = 0
            for uniprot_aa in uniprot_sequence:
                if uniprot_aa != peptide[j]:
                    structure_aa = peptide[j]
                    sequence_location = sequence_starts + j + 1
                    peptide_position = f'P{j + 1}'
                    print (f'MUTATION at {peptide_position} from {uniprot_aa} to {structure_aa} ({ sequence_location })')
                    differences[f'P{j+1}'] = {
                        'uniprot_aa':uniprot_aa,
                        'structure_aa': structure_aa,
                        'sequence_location': sequence_location
                    }
                    changes += 1
                j += 1
            print (f'MATCH STARTS AT {sequence_starts}')
            print (f'MATCH ENDS AT {sequence_ends}')
            print (f'MATCHED SEQUENCE {uniprot_sequence}')
            match['differences'] = differences
            match['changes'] = changes
        
    if sequence_starts:
        match['sequence_starts'] = sequence_starts
        match['sequence_ends'] = sequence_ends
        match['uniprot_sequence'] = str(fasta_record.seq[sequence_starts:sequence_ends])
        
    print ('MATCH BELOW')
    print (match)
    return match


def match_chain(squence:str, fasta_record):
    return {}


def fetch_uniprot_fasta(uniprot_id:str):
    fasta_record = None
    if uniprot_id:
        url = f'https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta'
        uniprot_record = http.get(url, 'txt')
        if uniprot_record:
            fasta_record = fasta_reader(uniprot_record)
            print (fasta_record.description)
    return fasta_record


def get_rcsb_info(pdb_code:str, molecule:int) -> str:
    rcsb_data, success, errors = rcsbProvider(pdb_code).fetch_uniprot(i)
    rcsb_info = None
    if len(rcsb_data) > 0:
        rcsb_data = rcsb_data[0]
        rcsb_info = {
            'uniprot_id': rcsb_data['rcsb_id'],
            'source_organism': rcsb_data['rcsb_uniprot_protein']['source_organism'],
            'protein_name': rcsb_data['rcsb_uniprot_protein']['name']
        }
    return rcsb_info




with Progress() as progress:
        
    task = progress.add_task("[white]Processing...", total=len(pdb_codes))

    print ('')
    for pdb_code in pdb_codes:

        print (f'{pdb_code}\n')
        
        uniprot_info = {}

        molecules_info, success, errors = PDBeProvider(pdb_code).fetch_molecules()
        
        assigned_chains = load_facet(pdb_code, 'assigned_chains')

        if not assigned_chains:
            print ('NO ASSIGNED CHAINS?')
        else:
            i = 1
            for chain in molecules_info:
                uniprot_id = None
                fasta_record = None
                if 'molecule' not in chain:
                    if 'length' in chain:
                        rcsb_info = get_rcsb_info(pdb_code, i)
                        if rcsb_info:                        
                            for assigned_chain in assigned_chains:
                                if assigned_chains[assigned_chain]['chains'] == chain['in_chains']:
                                    print (f'{assigned_chain}:\n')

                                    uniprot_id = rcsb_info['uniprot_id']
                                    uniprot_info[assigned_chain] = rcsb_info
                                    fasta_record = fetch_uniprot_fasta(uniprot_id)
                                    if fasta_record:
                                        if assigned_chain in ['peptide', 'peptide_fragment1']:
                                            match_details = match_peptide(assigned_chains[assigned_chain]['sequence'], fasta_record)
                                        else:
                                            if trim_sequence(assigned_chains[assigned_chain]['sequence']) in fasta_record.seq:
                                                print ('MATCH IN SEQUENCE\n')
                                            else:
                                                print ('MISMATCH IN SEQUENCE\n')
                i += 1
                    

        print (uniprot_info)




