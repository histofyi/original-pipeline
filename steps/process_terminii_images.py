from PIL import Image
from typing import List

import os

from functions.files import write_file
from functions.cli import load_config


def process_pocket_images(assembly_name, image_input_folder, image_output_folder):
    print (assembly_name)
    errors = None
    for size in ['large','medium']:
        try:
            terminii_image = f'{image_input_folder}/{assembly_name}_terminii_{size}.png'
            terminii_labels = f'assets/labels/terminii_labels_{size}.png'
            
            labelled = f'{image_output_folder}/{assembly_name}_labelled_{size}.png'

            labels_img = Image.open(terminii_labels)
            combined_img = Image.open(terminii_image)

            combined_img.paste(labels_img, (0, 0), labels_img)
            combined_img.save(labelled)

            print (size)

        except:
            errors = 'no_files_for_'+ assembly_name
            print (errors)
    print ('---')
    return errors



def perform_action(to_process:List, mhc_class:str, image_input_folder:str, image_output_folder:str):
    step_errors = []

    for assembly_name in to_process:
        this_errors = process_pocket_images(assembly_name, image_input_folder, image_output_folder)
        if this_errors:
            step_errors.append(assembly_name)
    return step_errors

def main():
    config = load_config()
    
    input_folder = f'{config["WAREHOUSE_PATH"]}/structures/coordinates/public/class_i/without_solvent/antigen_binding_domains'
    
    image_input_folder = f'{config["IMAGES_PATH"]}/cleft/terminii/unlabelled'
    image_output_folder = f'{config["IMAGES_PATH"]}/cleft/terminii/labelled'


    structures = [file.split('.')[0] for file in os.listdir(input_folder) if len(file.split('.')[0]) > 0]
    #structures = ['1hhg_1','1hhh_1','1hhi_1','1hhk_1']

    step_errors = perform_action(structures, 'class_i', image_input_folder, image_output_folder)

    print (step_errors)
main()    


