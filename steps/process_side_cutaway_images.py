from PIL import Image
from typing import List

import os

from functions.files import write_file
from functions.cli import load_config


def process_pocket_images(assembly_name, image_input_folder, image_output_folder):
    print (assembly_name)
    errors = None
    for size in ['large','medium']:
        if size:
            cutaway_image = f'{image_input_folder}/cutaway/{assembly_name}_cutaway_{size}.png'
            peptide_image = f'{image_input_folder}/peptide/{assembly_name}_peptide_{size}.png'

            labelled = f'{image_output_folder}/{assembly_name}_combined_{size}.png'

            combined_img = Image.open(cutaway_image)
            try:
                peptide_img = Image.open(peptide_image)
                combined_img.paste(peptide_img, (0, 0), peptide_img)
            except:
                print ('no peptide')

            combined_img.save(labelled)
            print (size)
    print ('---')
    return errors



def perform_action(to_process:List, mhc_class:str, image_input_folder:str, image_output_folder:str):
    step_errors = []

    for assembly_name in to_process:
        assembly_id = assembly_name.split('_')[1]
        if assembly_id == '1':
            this_errors = process_pocket_images(assembly_name, image_input_folder, image_output_folder)
            if this_errors:
                step_errors.append(assembly_name)
    return step_errors


def main():
    config = load_config()
    
    input_folder = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/without_solvent/antigen_binding_domains'
    
    image_input_folder = f'{config["IMAGES_PATH"]}/cleft/side/'
    image_output_folder = f'{config["IMAGES_PATH"]}/cleft/side/combined'


    structures = [file.split('.')[0] for file in os.listdir(input_folder) if len(file.split('.')[0]) > 0]
    #structures = ['1hhg_1','1hhh_1','1hhi_1','1hhk_1']

    step_errors = perform_action(structures, 'class_i', image_input_folder, image_output_folder)

    print (step_errors)
main()    


