from Bio.SeqIO.FastaIO import FastaIterator

def fasta_reader(filename):
    with open(filename) as handle:
        for record in FastaIterator(handle):
            yield record





def parse_allele_description(description, hla=False):
    description_elements = description.split(' ')
    if hla:
        prefix = 'HLA-'

    raw_allele = description_elements[1]
    allele_name = f'HLA-{raw_allele}'
    locus = f'HLA-{raw_allele.split("*")[0]}'

    allele_name = ':'.join(allele_name.split(":")[0:2])
    identifier = description_elements[0].split(':')[1]
    return allele_name, locus, identifier


