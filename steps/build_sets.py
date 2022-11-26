from typing import Dict, List, Tuple


from functions.files import write_json, read_json, load_constants
from functions.helpers import slugify
from functions.cli import load_config

import datetime

config = load_config()

index_filepath = f'{config["WAREHOUSE_PATH"]}/datacompilations/index.json'

core_filepath = f'{config["WAREHOUSE_PATH"]}/datacompilations/core.json'

sets_filepath = f'{config["WAREHOUSE_PATH"]}/datacompilations/sets.json'

complex_types_raw = load_constants('complex_types')
peptide_feature_labels = load_constants('peptide_features')

complex_types = {}

for component_number in complex_types_raw:
    for complex_type in complex_types_raw[component_number]:
        slug = complex_type['slug']
        if not slug in complex_types:
            complex_types[slug] = complex_type
print (complex_types)


indexes = read_json(index_filepath)
cores = read_json(core_filepath)


default_ordering = indexes['release_date_desc']


def build_last_updated_date():
    return datetime.datetime.now().isoformat()


def build_blank_set():
    return {
        'title':'',
        'description':'',
        'members':[],
        'count': 0,
        'slug':''
    }

sets = {
    "website":{}
}


def check_classical_class_i(core):
    if core['class'] == 'class_i' and core['classical']:
        return True
    else:
        return False


def build_species_sets():
    species_sets = {}
    for pdb_code in default_ordering:
        core = cores[pdb_code]
        species_slug = slugify(core['species']['common_name'])
        if not species_slug in species_sets:
            species_sets[species_slug] = build_blank_set()
            species_sets[species_slug]['slug'] = species_slug
            species_sets[species_slug]['title'] = f'{core["species"]["common_name"]} Class I structures'
        species_sets[species_slug]['members'].append(pdb_code)
        species_sets[species_slug]['count'] += 1
    return species_sets


def build_loci_sets():
    loci_sets = {}
    for pdb_code in default_ordering:
        core = cores[pdb_code]
        if check_classical_class_i(core):
            locus = slugify(core['allele']['alpha']['locus'])
            if not locus in loci_sets:
                loci_sets[locus] = build_blank_set()
            loci_sets[locus]['slug'] = locus
            loci_sets[locus]['title'] = f'{core["allele"]["alpha"]["locus"].upper()} structures'
            loci_sets[locus]['members'].append(pdb_code)
            loci_sets[locus]['count'] += 1
    return loci_sets


def build_allele_group_sets():
    allele_group_sets = {}
    for pdb_code in default_ordering:
        core = cores[pdb_code]
        if check_classical_class_i(core):
            if ':' in core['allele']['alpha']['name']:
                allele_group = core['allele']['alpha']['name'].split(':')[0]
            else:
                allele_group = core['allele']['alpha']['name']
            allele_group_slug = slugify(allele_group)
            if not allele_group_slug in allele_group_sets:
                allele_group_sets[allele_group_slug] = build_blank_set()
            allele_group_sets[allele_group_slug]['slug'] = allele_group_slug
            allele_group_sets[allele_group_slug]['title'] = f'{allele_group} structures'
            allele_group_sets[allele_group_slug]['members'].append(pdb_code)
            allele_group_sets[allele_group_slug]['count'] += 1
    return allele_group_sets


def build_allele_sets():
    allele_sets = {}
    for pdb_code in default_ordering:
        core = cores[pdb_code]
        if check_classical_class_i(core):
            allele_slug = core['allele']['alpha']['slug']
            if not allele_slug in allele_sets:
                allele_sets[allele_slug] = build_blank_set()
            allele_sets[allele_slug]['slug'] = allele_slug
            allele_sets[allele_slug]['title'] = f'{core["allele"]["alpha"]["name"]} structures'
            allele_sets[allele_slug]['members'].append(pdb_code)
            allele_sets[allele_slug]['count'] += 1
    return allele_sets
        

def build_peptide_sequence_sets():
    peptide_sequence_sets = {}
    for pdb_code in default_ordering:
        core = cores[pdb_code]
        if 'peptide_sequence' in core:
            if core['peptide_sequence'] is not None:
                peptide_sequence = core['peptide_sequence'].lower()
                if not peptide_sequence in peptide_sequence_sets:
                    peptide_sequence_sets[peptide_sequence] = build_blank_set()
                    peptide_sequence_sets[peptide_sequence]['slug'] = peptide_sequence.lower()
                    peptide_sequence_sets[peptide_sequence]['title'] = f'Structures containing "{peptide_sequence.upper()}"'
                peptide_sequence_sets[peptide_sequence]['members'].append(pdb_code)
                peptide_sequence_sets[peptide_sequence]['count'] += 1
    return peptide_sequence_sets


def build_peptide_length_sets():
    peptide_length_sets = {}
    for pdb_code in default_ordering:
        core = cores[pdb_code]
        if 'peptide_length_name' in core:
            if core['peptide_length_name'] is not None:
                peptide_length = core['peptide_length_name'].lower()
                if not peptide_length in peptide_length_sets:
                    peptide_length_sets[peptide_length] = build_blank_set()
                    peptide_length_sets[peptide_length]['slug'] = peptide_length
                    peptide_length_sets[peptide_length]['title'] = f'{peptide_length.capitalize()} peptide structures'
                peptide_length_sets[peptide_length]['members'].append(pdb_code)
                peptide_length_sets[peptide_length]['count'] += 1
    return peptide_length_sets

all_peptide_features = []


def build_peptide_feature_sets():
    peptide_feature_sets = {}
    for pdb_code in default_ordering:
        core = cores[pdb_code]
        if len(core['peptide']) > 1:
            i = 0
            for peptide in core['peptide']:
                if core['peptide'][peptide] is not None:
                    if 'features' in core['peptide'][peptide] and i == 0:
                        peptide_features = core['peptide'][peptide]['features']
                        for peptide_feature in peptide_features:
                            if peptide_feature not in all_peptide_features:
                                all_peptide_features.append(peptide_feature)
                            if peptide_feature not in peptide_feature_sets:
                                peptide_feature_sets[peptide_feature] = build_blank_set()
                                peptide_feature_sets[peptide_feature]['slug'] = peptide_feature
                                peptide_feature_sets[peptide_feature]['title'] = peptide_feature_labels[peptide_feature]
                            peptide_feature_sets[peptide_feature]['members'].append(pdb_code)
                            peptide_feature_sets[peptide_feature]['count'] += 1
                i += 1
    return peptide_feature_sets


def build_complex_type_sets():
    complex_type_sets = {}
    for pdb_code in default_ordering:
        core = cores[pdb_code]
        if 'complex_type' in core:
            if core['complex_type'] is not None:
                complex_type = core['complex_type']
                if not complex_type in complex_type_sets:
                    complex_type_sets[complex_type] = build_blank_set()
                    complex_type_sets[complex_type]['slug'] = complex_type
                    complex_type_sets[complex_type]['title'] = complex_types[complex_type]['label']                 
                complex_type_sets[complex_type]['members'].append(pdb_code)
                complex_type_sets[complex_type]['count'] += 1
    return complex_type_sets


def build_chronology_sets(chronology_type, chronology_title):
    chronology_type_sets = {}
    for pdb_code in default_ordering:
        core = cores[pdb_code]
        chronology_year = str(core['chronology'][f'{chronology_type}_year'])
        if not chronology_year in chronology_type_sets:
            chronology_type_sets[chronology_year] = build_blank_set()
            chronology_type_sets[chronology_year]['slug'] = str(chronology_year)
            chronology_type_sets[chronology_year]['title'] = f'Structures {chronology_title} in {chronology_year}'                
        chronology_type_sets[chronology_year]['members'].append(pdb_code)
        chronology_type_sets[chronology_year]['count'] += 1
    return chronology_type_sets


def build_resolution_sets():
    resolution_sets = {}
    for pdb_code in default_ordering:
        core = cores[pdb_code]
        if core['resolution'] is not None:
            resolution_low = int(core['resolution'].split('.')[0])
            resolution_high = resolution_low + 0.99
            resolution = f'{resolution_low}_{resolution_high}'         
            if not resolution in resolution_sets:
                resolution_sets[resolution] = build_blank_set()
                resolution_sets[resolution]['slug'] = resolution
                resolution_sets[resolution]['title'] = f'Structures between {resolution_low} and {resolution_high} &#8491;resolution'                
            resolution_sets[resolution]['members'].append(pdb_code)
            resolution_sets[resolution]['count'] += 1
    return resolution_sets


def build_latest_set():
    this_set = {
        'title': 'Latest structures',
        'description':f'Most recently released structures. Last updated {build_last_updated_date()}',
        'members':default_ordering[:10],
        'count':10,
        'slug':'latest'
    }
    return this_set


def build_all_set():
    this_set = {
        'title': 'All structures',
        'description':f'All released structures. Last updated {build_last_updated_date()}',
        'members':default_ordering,
        'count':len(default_ordering),
        'slug':'all'
    }
    return this_set



sets['website']['latest'] = build_latest_set()
sets['website']['all'] = build_all_set()

sets['species'] = build_species_sets()
sets['loci'] = build_loci_sets()
sets['allele_groups'] = build_allele_group_sets()
sets['alleles'] = build_allele_sets()

sets['peptide_lengths'] = build_peptide_length_sets()
sets['peptide_sequences'] = build_peptide_sequence_sets()
sets['peptide_features'] = build_peptide_feature_sets()

sets['complex_types'] = build_complex_type_sets()

sets['deposited'] = build_chronology_sets('deposition', 'deposited')
sets['revised'] = build_chronology_sets('revision', 'revised')
sets['released'] = build_chronology_sets('release', 'released')

sets['resolutions'] = build_resolution_sets()

write_json(sets_filepath, sets, verbose=True)
