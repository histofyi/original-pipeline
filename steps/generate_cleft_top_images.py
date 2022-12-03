from typing import List

import os

from functions.files import write_file
from functions.cli import load_config

orientations = {
        'top':{
            'class_i':(\
                0.995811403,    0.028724836,    0.086771332,\
                -0.087024398,    0.007673016,    0.996177554,\
                0.027949072,   -0.999556422,    0.010141926,\
                -0.000007659,   -0.000004176, -150.956878662,\
                -41.813217163,   60.248783112,   63.533638000,\
                149.222396851,  152.691375732,  -20.000000000 ),
            'peptide': (\
                0.995811403,    0.028724836,    0.086771332,\
                -0.087024398,    0.007673016,    0.996177554,\
                0.027949072,   -0.999556422,    0.010141926,\
                -0.000007659,   -0.000004176, -150.956878662,\
                -41.813217163,   60.248783112,   63.533638000,\
                -1795.062377930, 2096.976074219,  -20.000000000 )
        },
        'side': {
            'class_i':(0.999964476,0.001130534,-0.008345760,-0.001050179,0.999953210,0.009623017,0.008356228,-0.009613929,0.999918759,-0.000013337,-0.000001206,-169.672134399,-42.279060364,53.602447510,63.292312622,168.412094116,170.932174683,-20.000000000),
            'peptide':(0.999972939,0.007359593,0.000000000,-0.007359593,0.999972939,0.000000000,0.000000000,0.000000000,1.000000000,-0.000002027,0.000003786,-169.672134399,-42.265258789,53.586551666,61.640075684,158.271789551,181.072418213,-20.000000000)
        }
    }


def cleft_top_image(assembly_name, has_peptide, abd_folder, peptide_folder, abd_output_folder, peptide_output_folder):
    image_string = ""

    orientation = 'top'

    abd_load_string = f"cmd.load('{abd_folder}/{assembly_name}.pdb', 'class_i')\n"
    peptide_load_string = f"cmd.load('{peptide_folder}/{assembly_name}.pdb', 'peptide')\n"

    abd_orientation_string = f'cmd.set_view({orientations[orientation]["class_i"]})'

    abd_setup_string = """
cmd.hide(representation="cartoon", selection="class_i")
cmd.hide(representation="sticks", selection="class_i")
cmd.hide(representation="spheres", selection="class_i")
cmd.show(representation="mesh", selection="class_i")
cmd.color("gray80","class_i")
"""

    abd_hide_peptide_string = 'cmd.hide(representation="cartoon", selection="peptide")'

    abd_save_large = f"cmd.png('{abd_output_folder}/{assembly_name}_large.png', width=2200, height=2200, dpi=300, ray=1, quiet=0)\n"
    abd_save_medium = f"cmd.png('{abd_output_folder}/{assembly_name}_medium.png', width=1000, height=1000, dpi=300, ray=1, quiet=0)\n"
    abd_save_small = f"cmd.png('{abd_output_folder}/{assembly_name}_small_.png', width=500, height=500, dpi=300, ray=1, quiet=0)\n"
    abd_save_thumb = f"cmd.png('{abd_output_folder}/{assembly_name}_thumb_.png', width=250, height=250, dpi=300, ray=1, quiet=0)\n"

    peptide_setup_string = """
cmd.delete('class_i')
cmd.show(representation="sticks", selection="peptide")
util.cbay("peptide")
cmd.remove("hydro")
"""

    peptide_orientation_string = f'cmd.set_view({orientations[orientation]["peptide"]})'

    peptide_save_large = f"cmd.png('{peptide_output_folder}/{assembly_name}_large.png', width=2200, height=2200, dpi=300, ray=1, quiet=0)\n"
    peptide_save_medium = f"cmd.png('{peptide_output_folder}/{assembly_name}_medium.png', width=1000, height=1000, dpi=300, ray=1, quiet=0)\n"
    peptide_save_small = f"cmd.png('{peptide_output_folder}/{assembly_name}_small_.png', width=500, height=500, dpi=300, ray=1, quiet=0)\n"
    peptide_save_thumb = f"cmd.png('{peptide_output_folder}/{assembly_name}_thumb_.png', width=250, height=250, dpi=300, ray=1, quiet=0)\n"


    cleanup_string = "cmd.delete('all')\n"

    image_string = ""

    image_string += abd_load_string
    if has_peptide:
        image_string += peptide_load_string

    image_string += abd_orientation_string
    image_string += abd_setup_string

    if has_peptide:
        image_string += abd_hide_peptide_string

    image_string += abd_save_large
    image_string += abd_save_medium
    image_string += abd_save_small
    image_string += abd_save_thumb

    if has_peptide:
        image_string += peptide_setup_string
        image_string += peptide_orientation_string
        image_string += peptide_save_large
        image_string += peptide_save_medium
        image_string += peptide_save_small
        image_string += peptide_save_thumb

    image_string += cleanup_string

    return image_string

def iterator_cleanup():
    pass



def action_cleanup(pipeline_path:str):
    pass



def perform_action(to_process:List, mhc_class:str, abd_input_folder:str, abd_output_folder:str, peptide_input_folder:str, peptide_output_folder:str ):
    """
    Iterates through the list of items to process and performs the action
    """
    step_errors = []
    
    all_image_string = ""

    for assembly_name in to_process:
        assembly_id = assembly_name.split('_')[1]
        pdb_code = assembly_name.split('_')[0]


        if assembly_id == '1':
            this_image_string = cleft_top_image(assembly_name, abd_input_folder, peptide_input_folder, abd_output_folder, peptide_output_folder)
            iterator_cleanup()
            all_image_string += this_image_string

        
    return all_image_string




def main():
    config = load_config()
    
    abd_input_folder = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/without_solvent/antigen_binding_domains'    
    abd_output_folder = f'{config["IMAGES_PATH"]}/cleft/top/cutaway'

    peptide_input_folder = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/without_solvent/peptide'
    peptide_output_folder = f'{config["IMAGES_PATH"]}/cleft/top/peptide'



    structures = [file.split('.')[0] for file in os.listdir(abd_input_folder) if len(file.split('.')[0]) > 0]
    structures = ['1hhg_1','1hhh_1','1hhi_1','1hhk_1']

    image_string = perform_action(structures, 'class_i', abd_input_folder, abd_output_folder, peptide_input_folder, peptide_output_folder)

    image_string+= 'cmd.quit()'

    write_file(f'{config["TMP_PATH"]}/pockets.py',image_string)

main()    

