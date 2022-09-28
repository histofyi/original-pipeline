from typing import Dict, List, Tuple


from common.providers import rcsbProvider, PDBeProvider, httpProvider
from Bio import SeqIO
from io import StringIO


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

print (len(pdb_codes))