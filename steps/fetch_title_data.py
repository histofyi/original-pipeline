from typing import Dict, List, Tuple


from common.providers import PDBeProvider

from functions.pdb import load_pdb_lists
from functions.cli import load_config, print_spacer
from functions.files import load_constants, load_facet, write_facet
from functions.helpers import parse_date_to_isoformat

from rich.progress import Progress
from rich.console import Console
from rich import print_json

console = Console()
config = load_config()

completed = []
failed = []

mhc_class = 'class_i'

pdb_codes = load_pdb_lists('class_i', config['WAREHOUSE_PATH'], console)

#pdb_codes = ['1hhk']

with Progress() as progress:
        
    task = progress.add_task("[white]Processing...", total=len(pdb_codes))

    for pdb_code in pdb_codes:

        chronology_data = load_facet(pdb_code, 'titles')

        if not chronology_data:

            print (pdb_code)

            summary, success, errors = PDBeProvider(pdb_code).fetch_summary()

            title_info = {}

            if summary:
                
                title_info['pdb_title_original'] = summary['title']
                title_info['pdb_title_capitalized'] = summary['title'].capitalize()
                
                
                print (title_info) 

                if len(title_info) > 0:
                    write_facet(pdb_code, 'titles', title_info)
                    completed.append(pdb_code)
                else:
                    failed.append(pdb_code)
            
                print_spacer()
        else:
            completed.append(pdb_code)

            

        progress.update(task, advance=1)

print (len(pdb_codes))
print (len(completed))
print (len(failed))