from typing import Dict, List, Tuple


from common.providers import PDBeProvider

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

mhc_class = 'class_i'

pdb_codes = load_pdb_lists('class_i', config['WAREHOUSE_PATH'], console)

#pdb_codes = ['1hhk','6eny','6mpp']

with Progress() as progress:
        
    task = progress.add_task("[white]Processing...", total=len(pdb_codes))

    for pdb_code in pdb_codes:

        print (pdb_code)

        experiment_info, success, errors = PDBeProvider(pdb_code).fetch_experiment()

        abridged_info = {}
        
        if experiment_info:
            if slugify(experiment_info['experimental_method']) == 'x_ray_diffraction':
                abridged_info['resolution'] = experiment_info['resolution_high']
                abridged_info['spacegroup'] = experiment_info['spacegroup']
                abridged_info['cell'] = experiment_info['cell']
            elif slugify(experiment_info['experimental_method']) == 'electron_microscopy':
                abridged_info['resolution'] = experiment_info['resolution']
            else:
                abridged_info['resolution'] = None
                #abridged_info['nmr_spectrum_type'] = experiment_info['spectrum_type']
            abridged_info['experimental_method'] = experiment_info['experimental_method']

        print (abridged_info)
        if len(abridged_info) > 0:
            write_facet(pdb_code, 'experiment', abridged_info)
            completed.append(pdb_code)
        else:
            failed.append(pdb_code)
        
        print_spacer()

        progress.update(task, advance=1)

print (len(pdb_codes))
print (len(completed))
print (len(failed))