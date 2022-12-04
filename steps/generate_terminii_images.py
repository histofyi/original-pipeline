from typing import List

import os

from functions.files import write_file, load_facet
from functions.cli import load_config

terminii = {
        "pn": {"color":"pn_col", "residues":["7","59","171"]},
        "p1": {"color":"p1_col", "residues":["159"]},
        "p2": {"color":"p2_col", "residues":["63"]},
        "pc-1": {"color":"pc-1_col", "residues":["147",]},
        "pc": {"color":"pc_col", "residues":["80","84","123","143","146"]}
}

groups = []
for group in terminii:
    selection_string = ''
    length = len(terminii[group]['residues'])
    i = 0
    for residue in terminii[group]['residues']:
        if i == length - 1:
            selection_string += f' resi {residue}'
        else:
            selection_string += f' resi {residue},'
        i += 1
    groups.append({'name':group, 'color': terminii[group]['color'],'selection':selection_string})

def terminii_image(assembly_name, has_peptide, abd_folder, peptide_folder, output_folder):
    image_string = ""

    abd_load_string = f"cmd.load('{abd_folder}/{assembly_name}.pdb', 'class_i')\n"
    peptide_load_string = f"cmd.load('{peptide_folder}/{assembly_name}.pdb', 'peptide')\n"
    


    setup_string = """
cmd.set_view((\
0.997432590,    0.070144698,    0.014454702,\
-0.023742054,    0.133436739,    0.990772903,\
0.067568570,   -0.988572299,    0.134759426,\
0.000000000,    0.000000000, -167.991149902,\
-42.288372040,   56.106903076,   63.375183105,\
132.445495605,  203.536804199,  -20.000000000 ))
cmd.color('grey80','all')

cmd.set_color('pn_col', [120,148,177])
cmd.set_color('p1_col', [137,141,160])
cmd.set_color('p2_col', [156,139,151])
cmd.set_color('pc-1_col', [190, 160, 164])
cmd.set_color('pc_col', [239, 176, 168])

cmd.bg_color('white')
"""

    terminii_items_string = ""
    for group in groups:
        terminii_items_string += f'cmd.select("{group["name"]}", "{group["selection"]}")\n'
        terminii_items_string += f'cmd.color("{group["color"]}", "{group["name"]}")\n'
        terminii_items_string += f'cmd.show("sticks", "{group["name"]} and not (name c,n)")\n'
    terminii_items_string += 'cmd.select("none")\n'


    peptide_setup_string = """
cmd.show(representation="sticks", selection="peptide")
cmd.hide('sticks', 'peptide and not (name c,n,ca,o)')
cmd.hide('cartoon', 'peptide')
cmd.remove("hydro")
cmd.color("grey20","peptide")
"""

    save_large = f"cmd.png('{output_folder}/{assembly_name}_terminii_large.png', width=2200, height=2200, dpi=300, ray=1, quiet=0)\n"
    save_medium = f"cmd.png('{output_folder}/{assembly_name}_terminii_medium.png', width=1000, height=1000, dpi=300, ray=1, quiet=0)\n"
    save_small = f"cmd.png('{output_folder}/{assembly_name}_terminii_small_.png', width=500, height=500, dpi=300, ray=1, quiet=0)\n"
    save_thumb = f"cmd.png('{output_folder}/{assembly_name}_terminii_thumb_.png', width=250, height=250, dpi=300, ray=1, quiet=0)\n"

    cleanup_string = "cmd.delete('all')\n"

    image_string = ""

    if has_peptide:
        image_string += peptide_load_string
    image_string += abd_load_string
    
    image_string += setup_string
    image_string += terminii_items_string
    if has_peptide:
        image_string += peptide_setup_string
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



def perform_action(to_process:List, mhc_class:str, abd_input_folder:str, peptide_input_folder:str, output_folder:str):
    """
    Iterates through the list of items to process and performs the action
    """
    step_errors = []
    
    all_image_string = ""

    for assembly_name in to_process:
        assembly_id = assembly_name.split('_')[1]
        pdb_code = assembly_name.split('_')[0]

        assigned_chains = load_facet(pdb_code, 'assigned_chains')

        if 'peptide' in assigned_chains:
            has_peptide = True

            if not os.path.exists(f'{peptide_input_folder}/{assembly_name}.pdb'):
                has_peptide = False

        else:
            has_peptide = False

        if assembly_id == '1':
            this_image_string = terminii_image(assembly_name, has_peptide, abd_input_folder, peptide_input_folder, output_folder)
            iterator_cleanup()
            all_image_string += this_image_string

        
    return all_image_string




def main():
    config = load_config()
    
    abd_input_folder = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/without_solvent/antigen_binding_domains'
    peptide_input_folder = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/without_solvent/peptide'
    output_folder = f'{config["IMAGES_PATH"]}/cleft/terminii/unlabelled'


    structures = [file.split('.')[0] for file in os.listdir(abd_input_folder) if len(file.split('.')[0]) > 0]
    #structures = ['1hhg_1','1hhh_1','1hhi_1','1hhk_1']

    image_string = perform_action(structures, 'class_i', abd_input_folder, peptide_input_folder, output_folder)

    image_string+= 'cmd.quit()'

    write_file(f'{config["TMP_PATH"]}/terminii.py',image_string)

main()    

