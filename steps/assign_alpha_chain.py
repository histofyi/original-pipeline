from typing import Dict, List, Tuple

from common.providers import PDBeProvider

from functions.pdb import load_pdb_lists
from functions.cli import load_config, print_spacer
from functions.files import load_constants, load_facet, write_facet, read_json
from functions.helpers import slugify, deduplicate_list

from fuzzywuzzy import fuzz

from rich.progress import Progress
from rich.console import Console



console = Console()
config = load_config()


mhc_class = 'class_i'


species = load_constants('species')
class_i_starts = load_constants('mhc_starts')[mhc_class]['alpha']
assign_alpha_chain_overrides = load_constants('assign_alpha_chain_overrides')
class_i_firsts = [aa for aa in set([aa[0] for aa in class_i_starts])]

pdb_codes = load_pdb_lists('class_i', config['WAREHOUSE_PATH'], console)

#pdb_codes = ['6ilc','6ile','6j2d', '1hhk']
#pdb_codes = ["1fo0", "1fzj", "1fzk", "1fzm", "1fzo", "1kj2", "1kj3", "1mwa", "1n59", "1nam", "1nan", "1rjy", "1rjz", "1t0m", "1t0n", "1wbz", "2clv", "2clz", "2ol3", "2zsv", "2zsw", "3fol", "3fom", "3fon", "3p4m", "3p4n", "3p4o", "3p9l", "3p9m", "3pab", "3rgv", "3rol", "3roo", "3tid", "3tie", "4hkj", "4hs3", "4pv8", "4pv9", "6gb6", "6ilc", "6ilf", "6ilg", "6j2d", "6j2e", "6j2f", "6j2g", "6j2h", "6j2i", "6j2j", "6k7t", "6wl2", "6wl3", "6wl4", "7ji2"]



overrides = []
exact_matched = []
fuzzy_matched = []
unmatched = []

not_classical = []
classical = []

missing_loci = []
missing_species = []

unaasignable = {
    'poor_match':[],
    'other':[]
}

facet_name = 'allele'


def fetch_species_nomenclature(species_slug:str, mhc_class:str, species:Dict) -> Tuple[str, List]:
    if species_slug in species:
        return species[species_slug]['stem'], species[species_slug][mhc_class]
    else:
        return None, None


def prep_match(sequence:str, test:str):
    if len(sequence) > 275:
        # long sequence, longer than normal
        sequence = sequence[:275]
    elif len(sequence) < 200:
        # truncated, only alpha1 and alpha2
        test = test[0:len(sequence)]
    elif len(sequence) <= 274:
        # often missing first amino acid
        if sequence[0] not in class_i_firsts:
            test = test[1:]
        else:
            # there's a truncation
            test = test[0:len(sequence)]
    return sequence, test


def do_exact_match(sequence:str, test:str):
    sequence, test = prep_match(sequence, test)
    if sequence == test:
        return True
    else:
        return False


def do_fuzzy_match(sequence:str, test:str):
    sequence, test = prep_match(sequence, test)
    ratio = fuzz.ratio(test, sequence) / 100
    return ratio

def trim_start(sequence:str):
    found = False
    for class_i_start in class_i_starts:
        if not found:
            if class_i_start in sequence:
                found = True
                start_index = sequence.index(class_i_start)
                if start_index > 0:
                    sequence = sequence[start_index:]
                break
    return sequence


cached_class_i_loci ={}

with Progress() as progress:
        
    task = progress.add_task("[white]Processing...", total=len(pdb_codes))

    for pdb_code in pdb_codes:
        print_spacer()
        print (pdb_code)
        assigned_chains = load_facet(pdb_code, 'assigned_chains')
        alpha_chain = f'{mhc_class}_alpha'
        if not alpha_chain in assigned_chains:
            not_classical.append(pdb_code)
            print ('NOT CLASS I')
        else:
            classical.append(pdb_code)
            species_info = load_facet(pdb_code, 'species')
            sequence = assigned_chains[alpha_chain]['sequence']
            sequence = trim_start(sequence)
            species_slug = species_info['slug']
            species_stem, loci_letters = fetch_species_nomenclature(species_slug, mhc_class, species)
            if not species_stem:
                print (f'MISSING SPECIES - {species_info["slug"]}')
                missing_species.append(species_info)
            else:
                loci = [f'{species_stem}_{locus}' for locus in loci_letters]
                matched = False
                exact_match = None
                best_match = None
                best_score = 0
                for locus in loci:
                    allele = None
                    if not matched:
                        if locus in cached_class_i_loci:
                            allele_sequences = cached_class_i_loci[locus]
                        else:
                            locuspath = f'{config["WAREHOUSE_PATH"]}/sequences/{locus}.json'
                            try:
                                allele_sequences = read_json(locuspath)
                                cached_class_i_loci[locus] = allele_sequences
                            except:
                                allele_sequences = None
                                missing_loci.append(locus)
                        if allele_sequences:
                            if not matched:
                                for test in allele_sequences:
                                    if not matched:
                                        if do_exact_match(sequence, allele_sequences[test]['sequence']):
                                            best_match = allele_sequences[test]
                                            matched = True
                                            exact_matched.append(pdb_code)
                            if not matched:
                                for test in allele_sequences:
                                    if not matched:
                                        score = do_fuzzy_match(sequence, allele_sequences[test]['sequence'])
                                        if score:
                                            if score > best_score:
                                                best_match = allele_sequences[test]
                                                best_score = score
                error_type = None
                match_type = None
                if not matched:
                    if best_match:
                        if best_score >= 0.95:
                            print (f'Good score: {best_match} : {best_score}')
                            fuzzy_matched.append(pdb_code)
                            match_type = 'fuzzy'
                        else:
                            print (f'Low score: {best_match} : {best_score}')
                            if pdb_code in assign_alpha_chain_overrides:
                                best_match = assign_alpha_chain_overrides[pdb_code]
                                overrides.append(pdb_code)
                                print (f'Override used')
                                match_type = 'override'
                            else:
                                error_type = 'poor_match'
                                unmatched.append(pdb_code)
                            
                    else:
                        if pdb_code in assign_alpha_chain_overrides:
                            best_match = assign_alpha_chain_overrides[pdb_code]
                            overrides.append(pdb_code)
                            print (f'Override used')
                            match_type = 'override'
                        else:
                            error_type = 'other'
                            print ('Unable to assign due to issues')
                            unmatched.append(pdb_code)
                else:
                    print (f'Exact match : {exact_match}')
                    match_type = 'exact'

                if error_type:
                    unaasignable[error_type].append({
                        'species':species_stem,
                        'loci':loci,
                        'sequence':sequence,
                        'pdb_code':pdb_code,
                        'best_match': best_match,
                        'best_score': best_score
            })
            
        if not error_type:
            best_match['match_type'] = f'histo:{match_type}'
            match_data = {}
            for key in best_match:
                if key != 'sequence':
                    match_data[key] = best_match[key]
            write_facet(pdb_code, 'alpha_chain_details', match_data)

        progress.update(task, advance=1)


        

print (f'Total {len(pdb_codes)}')

print (f'{len(deduplicate_list(not_classical))} not classical class I')
print (f'{len(deduplicate_list(classical))} classical class I')


print (f'{len(deduplicate_list(exact_matched))} assigned using exact match')
print (f'{len(deduplicate_list(fuzzy_matched))} assigned using fuzzy match')
print (f'{len(deduplicate_list(unmatched))} unassigned')
print (f'{len(overrides)} overrides used')


for error_type in ['other', 'poor_match']:


    for structure in unaasignable[error_type]:
        print (f'{structure["pdb_code"]}')
        print (f'{structure["species"]}')
        print (f'{structure["loci"]}')
        if structure["best_match"]:
            print (f'{structure["best_match"]} : {structure["best_score"]}')
        else:
            print ('Unable to match')
        print (f'{structure["sequence"]}')
        print_spacer()

print (deduplicate_list(unmatched))
