from typing import List

import os

from functions.files import write_file
from functions.cli import load_config

pockets = {
        "a": {"color":"wheat", "residues":["5","59","63","66","159","163","167","171"]},
        "b": {"color":"lightpink", "residues":["7","9","24","25","33","34","45","60","67","70"]},
        "c": {"color":"palecyan", "residues":["73","74"]},
        "d": {"color":"palegreen", "residues":["99","114","155","156",]},
        "e": {"color":"lightblue", "residues":["97","114","147","152"]},
        "f": {"color":"lightorange", "residues":["77","80","81","84","95","116","123","143","146","147"]}
}

pymol_pockets = []
for pocket in pockets:
    selection_string = ''
    length = len(pockets[pocket]['residues'])
    i = 0
    for residue in pockets[pocket]['residues']:
        if i == length - 1:
            selection_string += f' resi {residue}'
        else:
            selection_string += f' resi {residue},'
        i += 1
    pymol_pockets.append({'letter':pocket.upper(), 'color': pockets[pocket]['color'],'selection':selection_string})


def pockets_image(assembly_name, input_folder, output_folder):
    image_string = ""

    load_string = f"cmd.load('{input_folder}/{assembly_name}.pdb')\n"

    setup_string = """
cmd.set_view((\
0.997432590,    0.070144698,    0.014454702,\
-0.023742054,    0.133436739,    0.990772903,\
0.067568570,   -0.988572299,    0.134759426,\
0.000000000,    0.000000000, -167.991149902,\
-42.288372040,   56.106903076,   63.375183105,\
132.445495605,  203.536804199,  -20.000000000 ))
cmd.color('grey80','all')
cmd.show('surface')
cmd.show('sticks', 'all and not (name c,n)')
cmd.set('transparency', 0.25, 'all')
cmd.bg_color('white')
"""

    pocket_items_string = ""
    for pocket in pymol_pockets:
        pocket_items_string += f'cmd.select("{pocket["letter"]}_pocket", "{pocket["selection"]}")\n'
        pocket_items_string += f'cmd.color("{pocket["color"]}", "{pocket["letter"]}_pocket")\n'
    pocket_items_string += 'cmd.select("none")\n'

    save_large = f"cmd.png('{output_folder}/{assembly_name}_pockets_large.png', width=2200, height=2200, dpi=300, ray=1, quiet=0)\n"
    save_medium = f"cmd.png('{output_folder}/{assembly_name}_pockets_medium.png', width=1000, height=1000, dpi=300, ray=1, quiet=0)\n"
    save_small = f"cmd.png('{output_folder}/{assembly_name}_pockets_small_.png', width=500, height=500, dpi=300, ray=1, quiet=0)\n"
    save_thumb = f"cmd.png('{output_folder}/{assembly_name}_pockets_thumb_.png', width=250, height=250, dpi=300, ray=1, quiet=0)\n"

    cleanup_string = "cmd.delete('all')\n"

    image_string = ""

    image_string += load_string
    image_string += setup_string
    image_string += pocket_items_string
    image_string += save_large
    image_string += save_medium
    image_string += save_small
    image_string += save_thumb
    image_string += cleanup_string

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
            this_image_string = pockets_image(assembly_name, input_folder, output_folder)
            iterator_cleanup()
            all_image_string += this_image_string

        
    return all_image_string




def main():
    config = load_config()
    
    input_folder = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/without_solvent/antigen_binding_domains'
    output_folder = f'{config["IMAGES_PATH"]}/cleft/pockets/unlabelled'


    structures = [file.split('.')[0] for file in os.listdir(input_folder) if len(file.split('.')[0]) > 0]
    structures = ['1hhg_1','1hhh_1','1hhi_1','1hhk_1']

    image_string = perform_action(structures, 'class_i', input_folder, output_folder)

    image_string+= 'cmd.quit()'

    write_file(f'{config["TMP_PATH"]}/pockets.py',image_string)

main()    

