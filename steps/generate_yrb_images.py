from typing import List
from pymol import cmd

import os

from functions.files import write_file
from functions.cli import load_config



def yrb_image(assembly_name, input_folder, output_folder):
    image_string = ""

    load_string = f"cmd.load('{input_folder}/{assembly_name}.pdb')"

    yrb_string = """
cmd.set_color('yellow',[0.950,0.78,0.0])
cmd.set_color('grey',[0.95,0.95,0.95])
cmd.set_color('red',[1.0,0.4,0.4])	
cmd.set_color('blue',[0.2,0.5,0.8])	

mapping = {}
mapping['arg'] = [ ('NE,NH2,NH1', 'blue'), ('CD,CZ', 'grey'), ('CG', 'yellow') ]
mapping['asn'] = [ ('CG,OD1,ND2', 'grey') ]
mapping['asp'] = [ ('CG', 'grey'), ('OD2,OD1', 'red')  ]
mapping['cys'] = [ ('SG', 'grey') ]	
mapping['gln'] = [ ('CG', 'yellow'), ('CD,OE1,NE2', 'grey') ]
mapping['glu'] = [ ('CG', 'yellow'), ('CD', 'grey'), ('OE1,OE2', 'red') ]
mapping['his'] = [ ('CG,CD2,ND1,NE2,CE1', 'grey') ]	
mapping['ile'] = [ ('CG1,CG2,CD1', 'yellow') ]
mapping['leu'] = [ ('CG,CD1,CD2', 'yellow') ]
mapping['lys'] = [ ('CG,CD', 'yellow'), ('CE', 'grey'), ('NZ', 'blue') ]
mapping['met'] = [ ('CG,CE', 'yellow'), ('SD', 'grey') ]
mapping['phe'] = [ ('CG,CD1,CE1,CZ,CE2,CD2', 'yellow') ]
mapping['pro'] = [ ('CG', 'yellow'), ('CD', 'grey') ]
mapping['ser'] = [ ('CB,OG', 'grey') ]
mapping['thr'] = [ ('CB,OG1', 'grey'), ('CG2', 'yellow') ]
mapping['trp'] = [ ('CG,CD2,CZ2,CH2,CZ3,CE3', 'yellow'), ('CD1,NE1,CE2', 'grey') ]
mapping['tyr'] = [ ('CG,CE1,CD1,CE2,CD2', 'yellow'), ('CZ,OH', 'grey') ]
mapping['val'] = [ ('CG1,CG2', 'yellow') ]

obj_list = cmd.get_names('objects')

for obj in obj_list:
    cmd.color('grey','(n. N,C,CA,O and ' + obj + ')')
    cmd.color('yellow','(n. CB and ' + obj + ')')
    
    for key in mapping:
        for (atom, color) in mapping[key]:
            cmd.color(color, '( n. ' + atom + ' and r. ' + key + ' and ' + obj + ' )')

cmd.show('surface', 'all')

cmd.set_view((\
    0.997432590,    0.070144698,    0.014454702,\
    -0.023742054,    0.133436739,    0.990772903,\
    0.067568570,   -0.988572299,    0.134759426,\
    0.000000000,    0.000000000, -167.991149902,\
    -42.288372040,   56.106903076,   63.375183105,\
    132.445495605,  203.536804199,  -20.000000000 ))
"""
    
    save_large = f"cmd.png('{output_folder}/{assembly_name}_yrb_large.png', width=2200, height=2200, dpi=300, ray=1, quiet=0)"
    save_medium = f"cmd.png('{output_folder}/{assembly_name}_yrb_medium.png', width=1000, height=1000, dpi=300, ray=1, quiet=0)"
    save_small = f"cmd.png('{output_folder}/{assembly_name}_yrb_small_.png', width=500, height=500, dpi=300, ray=1, quiet=0)"
    save_thumb = f"cmd.png('{output_folder}/{assembly_name}_yrb_thumb_.png', width=250, height=250, dpi=300, ray=1, quiet=0)"

    image_string = load_string + '\n'
    image_string += yrb_string
    image_string += save_large + '\n'
    image_string += save_medium + '\n'
    image_string += save_small + '\n'
    image_string += save_thumb + '\n'
    image_string += "cmd.delete('all')\n"


    return image_string


def iterator_cleanup():
    pass



def action_cleanup(pipeline_path:str):
    pass



def perform_action(to_process:List, mhc_class:str, input_folder:str, output_folder:str):
    """
    Iterates through the list of items to process and performs the action
    """
    step_errors = []
    
    all_image_string = ""

    for assembly_name in to_process:
        assembly_id = assembly_name.split('_')[1]
        pdb_code = assembly_name.split('_')[0]


        if assembly_id == '1':
            this_image_string = yrb_image(assembly_name, input_folder, output_folder)
            iterator_cleanup()
            all_image_string += this_image_string
    
        
    return all_image_string




def main():
    config = load_config()
    
    input_folder = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/without_solvent/antigen_binding_domains'
    output_folder = f'{config["IMAGES_PATH"]}/cleft/yrb'


    structures = [file.split('.')[0] for file in os.listdir(input_folder) if len(file.split('.')[0]) > 0]

    image_string = perform_action(structures, 'class_i', input_folder, output_folder)

    image_string+= 'cmd.quit()'

    write_file(f'{config["TMP_PATH"]}/yrb.py',image_string)

main()