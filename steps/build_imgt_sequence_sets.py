from typing import List

from functions.files import write_file, load_constants, write_json
from functions.fasta import fasta_reader, parse_allele_description
from functions.helpers import slugify
from functions.cli import load_config


input_file = 'tmp/sequences/imgt_hla.fasta'

loci = {}
sequences = {}

class_i_loci = ['hla_a','hla_b','hla_c','hla_e','hla_f','hla_g','hla_mica','hla_micb','hla_hfe']

classical = ['hla_a','hla_b','hla_c','hla_e','hla_f','hla_g']

class_i_starts = load_constants('mhc_starts')['class_i']['alpha']


def test_class_i_start(sequence:str, class_i_starts:List) -> bool:
    for start in class_i_starts:
        if start[1:] in sequence:
            if start in sequence:
                return start, True
            else:
                return start[1:], False
    return None, False
     


def trim_class_i_sequences(sequence:str, start:str, first_res:bool) -> str:
    leader_removed_sequence = f'{start}{sequence.split(start)[1]}'
    if first_res:
        end = 275
    else:
        end = 274
    trimmed_sequence = leader_removed_sequence[0:end] 
    return trimmed_sequence


def build_sequence_array(first_res:bool, trimmed_sequence:str) -> List:
    if not first_res:
        trimmed_sequence = '-' + trimmed_sequence
    sequence_array = [residue for residue in trimmed_sequence]
    return sequence_array


def show_loci(loci):
    for locus in loci:
        print (locus)
        this_locus = loci[locus]
        for allele in this_locus:
            this_allele = this_locus[allele]
            print (this_allele)


for entry in fasta_reader(input_file):
    allele_name, locus, identifier = parse_allele_description(entry.description, hla=True)
    locus_slug = slugify(locus)

    # first ensure that the locus slug is in those for Class I genes
    # TODO work on Class II genes
    if locus_slug in class_i_loci:
        # if the locus slug is a classical one (i.e. presents peptides)
        if locus_slug in classical:
            # create a locus in the loci dictionary
            if locus_slug not in loci:
                loci[slugify(locus)] = {}
            sequence = str(entry.seq)
            # ignore allele names with low expression "Q" and null alleles "N"
            if 'Q' not in allele_name and 'N' not in allele_name:
                # test whether the sequence is longer than the standard length for MHC Class I chains
                if len(sequence) > 275:
                    # now find if the sequence starts with one of the standard MHC Class I starts
                    start, first_res = test_class_i_start(sequence, class_i_starts)
                    if start:                    
                        allele_name_slug = slugify(allele_name)
                        if allele_name_slug not in loci[locus_slug]:
                            trimmed_sequence = trim_class_i_sequences(sequence, start, first_res)
                            loci[locus_slug][allele_name_slug] = {
                                'name':allele_name,
                                'locus':locus,
                                'allele_group':allele_name.split(':')[0],
                                'slug': allele_name_slug,
                                'identifier':f'ipd-imgt:{identifier}',
                                'species_stem':'hla',
                                'sequence':trimmed_sequence
                            }
                            sequence_array = build_sequence_array(first_res, trimmed_sequence)
                            sequences[allele_name_slug] = sequence_array




config = load_config()


sequences_path = f'{config["WAREHOUSE_PATH"]}/datacompilations/ipd_hla.json'

write_json(sequences_path, sequences, verbose=True)












