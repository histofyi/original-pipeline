from typing import Dict, List, Tuple

from common.providers import PDBeProvider

from functions.pdb import load_pdb_lists
from functions.cli import load_config, print_spacer
from functions.files import load_constants, load_facet, write_facet

from functions.helpers import slugify

from rich.progress import Progress
from rich.console import Console



console = Console()
config = load_config()

def get_species_info(organism_scientific:str, species:Dict) -> Dict:
    organism_slug = slugify(organism_scientific)
    if organism_slug in species:
        organism_info = {
            'scientific_name':species[organism_slug]['scientific_name'],
            'common_name':species[organism_slug]['common_name'],
            'slug':organism_slug
        }
    else:
        organism_info = None
    return organism_info

alpha_chains = ['class_i_alpha', 'mr1', 'cd1a', 'cd1b', 'cd1d', 'fcrn', 'mica', 'micb', 'hfe2', 'h2-t22','zag']


species = load_constants('species')
species_overrides = load_constants('species_overrides')

mhc_class = 'class_i'

pdb_codes = load_pdb_lists('class_i', config['WAREHOUSE_PATH'], console)


overrides = []
assigned = []
missing = []

facet_name = 'species'

with Progress() as progress:
        
    task = progress.add_task("[white]Processing...", total=len(pdb_codes))

    for pdb_code in pdb_codes:
        species_info = load_facet(pdb_code, facet_name)
        if not species_info:
            species_info = None
            if pdb_code in species_overrides:
                overrides.append(pdb_code)
                species_info = get_species_info(species_overrides[pdb_code]['organism_scientific_name'], species)
                species_info['match_type'] = 'histo:override'
            else:
                assigned_chains = load_facet(pdb_code, 'assigned_chains')
                print (pdb_code)
                for chain in assigned_chains:
                    if chain in alpha_chains:
                        molecules_info, success, errors = PDBeProvider(pdb_code).fetch_molecules()
                        for molecule in molecules_info:
                            if assigned_chains[chain]['chains'][0] in molecule['in_chains']:
                                if 'source' in molecule:
                                    assigned.append(pdb_code)
                                    species_info = get_species_info(molecule['source'][0]['organism_scientific_name'], species)
                                    species_info['match_type'] = 'histo:assign_species'
            if not species_info:
                missing.append(pdb_code)
            else:
                write_facet(pdb_code, facet_name, species_info)
        else:
            assigned.append(pdb_code)
        print (species_info)

        progress.update(task, advance=1)

print (f'{len(assigned)} assigned automatically')
print (f'{len(overrides)} assigned using overrides')
print (f'{len(missing)} unassigned')