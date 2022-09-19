from localpdb import PDB
import toml

from functions.files import write_json
from functions.pdb import load_pdb_lists
from rich.console import Console


console = Console()


config = toml.load('config.toml')

class_i = load_pdb_lists('class_i', config['WAREHOUSE_PATH'], console)


lpdb = PDB(db_path=config['LOCALPDB_PATH'])

alpha_fold_sets = {}

alpha_fold_date_cutoff = 20180430
alpha_fold_resolution_cutoff = 8

pdb_after_cutoff = lpdb.entries.query(f'deposition_date >= {alpha_fold_date_cutoff} & resolution <= 8')
pdb_before_cutoff = lpdb.entries.query(f'deposition_date < {alpha_fold_date_cutoff} & resolution <= 8')

print (len(pdb_after_cutoff))
print (len(pdb_before_cutoff))

mhc_after_cutoff = [pdb_code for pdb_code in pdb_after_cutoff.index if pdb_code in class_i]
mhc_before_cutoff = [pdb_code for pdb_code in pdb_before_cutoff.index if pdb_code in class_i]

print (len(class_i))

print (len(mhc_after_cutoff))
print (len(mhc_before_cutoff))





all_under_8 = set(mhc_after_cutoff) | set(mhc_before_cutoff)

print (len(all_under_8))

excluded = set(class_i) - all_under_8

print (len(excluded))

print (excluded)

alpha_fold_sets['before_date_cutoff'] = mhc_before_cutoff
alpha_fold_sets['after_date_cutoff'] = mhc_after_cutoff
alpha_fold_sets['over_8_angstroms'] = [pdb_code for pdb_code in excluded]

print (alpha_fold_sets)


filepath = f'{config["CONSTANTS"]}/alpha_fold_sets.json'
write_json(filepath, alpha_fold_sets, verbose=True, pretty=True)

