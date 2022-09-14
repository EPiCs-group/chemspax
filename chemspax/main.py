# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
import sys
import os
import random
import glob
import argparse
import logging
from gooey import Gooey
from pathlib import Path
from chemspax.attach_substituent import Complex
from chemspax.utilities import copy_functionalization_list_xyz_2_mol
CHEMSPAX_HOME_DIR=os.environ["CHEMSPAX_HOME_DIR"]
if CHEMSPAX_HOME_DIR == None:
    print("none is set")
else:
    print(CHEMSPAX_HOME_DIR)
def initialize_complex(original_skeleton_name, source_data, substituent_name, path_to_database):
    """Main function that takes command line arguments and passes them to attach_substituent.py's Complex class
    """
    complex = Complex(original_skeleton_name, source_data, substituent_name, path_to_database)
    return complex
path_to_database = os.path.join(CHEMSPAX_HOME_DIR, "substituents_xyz", "manually_generated", "central_atom_centroid_database.csv")
path_to_substituents = os.path.join(CHEMSPAX_HOME_DIR,"substituents_xyz","manually_generated/")
path_to_output = os.path.join(CHEMSPAX_HOME_DIR,"functionalized_complexes/")
@Gooey
def main(path_to_database=path_to_database, path_to_substituents = path_to_substituents):
    parser = argparse.ArgumentParser(prog='chemspax', description='Attach substituents to a skeleton molecule')
    parser.add_argument('-s', '--substituent', help='Name of substituent to attach', required=False, action='append', default=None)
    parser.add_argument('-k', '--skeleton', help='Provide path for skeletons in .xyz format', required=False,action='store', default=None)
    parser.add_argument('-o', '--output', help="path to output files", required=False, action='store',default=path_to_output)
    args = parser.parse_args()
    if args.substituent is not None:
        substituent_list = []
        for arg in args.substituent:
            substituent_list.append(arg)
    else:
        substituent_list = None
    if args.skeleton is not None:
        list_of_skeleton_files = glob.glob(args.skeleton+"/*.xyz")
    else:
        list_of_skeleton_files = glob.glob(CHEMSPAX_HOME_DIR+"/skeletons/*.xyz")
    skeleton_list = list_of_skeleton_files 
    # initialize logger
    Path(path_to_output).mkdir(parents=True, exist_ok=True)  #create directory for functionalized output
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    # create a file handler
    handler = logging.FileHandler('chemspax.log')
    handler.setLevel(logging.INFO)
    # create a logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(handler)
    # create a console handler
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # add the handlers to the logger
    logger.addHandler(console)

    logger.log(logging.INFO, "Looks like you want to add some substituents to a metal-ligand complex & optimize 'em.\n")
    logger.log(logging.INFO, "Created by: Adarsh Kalikadien & Vivek Sinha")

    # main loop to attach substituents to the skeleton
    for skeleton in skeleton_list:
        skeleton_path = os.path.normpath(skeleton)
        skeleton_name = os.path.basename(skeleton_path).split(".")[0]
        logger.log(logging.INFO, f"Attaching substituents to skeleton: {skeleton_name}")
        # intialize the complex class with a random substituent to be able to find lenght of functionalization_list
        first_target_file = skeleton_name + "_func_" + '1'
#        path_to_substituent_database = os.path.join("substituents_xyz", "manually_generated", "central_atom_centroid_database.csv")    #to keep things consistent let us just use "path_to_database" as used in "Complex" class as well"
        complex = initialize_complex(skeleton_name, skeleton_path, "CH3", path_to_database)
        if substituent_list is not None:
            logger.log(logging.INFO, f"initializing complex by placing {substituent_list[0]} on {skeleton_name}")
        else:
            logger.log(logging.INFO, f"initializing complex by placing random substituents based on length of functionalization_list")
            # create substituent_list with correct length (number of substituents)
            len_functionalization_list = len(complex.functionalization_site_list)
            # select random carbon (C) substituents, same length as functionalization_list (number of functionalization sites)
            all_substituents = glob.glob(path_to_substituents+"/C*.xyz")
            # only get name of substituent file without extension
            all_substituents = [os.path.basename(os.path.normpath(substituent)).split(".")[0] for substituent in all_substituents]
            substituent_list = random.choices(all_substituents, k=len_functionalization_list)
            logger.log(logging.INFO, f"initializing complex by placing {substituent_list[0]} on {skeleton_name}")


        # if functionalization_list is empty, skip this skeleton
        if len(complex.functionalization_site_list) == 0:
            logger.log(logging.INFO, f"No functionalization sites found for {skeleton_name}")
            continue
        # if substituent_list is too large, just use the first n elements of the list
        if len(substituent_list) > len(complex.functionalization_site_list):
            logger.log(logging.INFO, f"Substituent list is larger than amount of specified functionalization sites, "
                                     f"only using first {len(complex.functionalization_site_list)} elements")
            substituent_list = substituent_list[:len(complex.functionalization_site_list)]
        # if substituent_list is too small, add random substituents to the list
        if len(substituent_list) < len(complex.functionalization_site_list):
            logger.log(logging.INFO, f"Substituent list is smaller than amount of specified functionalization sites, "
                                     f"adding random carbon (C) substituents to the list")
            # get all C substituents from the database
            all_substituents = glob.glob("substituents_xyz/manually_generated/C*.xyz")
#            all_substituents = glob.glob(path_to_substituents+"/C*.xyz".split(".")[0])
            # only get name of substituent file without extension
            all_substituents = [os.path.basename(os.path.normpath(substituent)).split(".")[0] for substituent in all_substituents]
            # add random substituents to the list
            for i in range(len(substituent_list), len(complex.functionalization_site_list)):
                substituent_list.append(random.choice(all_substituents))

        complex = initialize_complex(skeleton_name, skeleton_path, substituent_list[0], path_to_database)
        complex.generate_substituent_and_write_xyz(first_target_file, 1.54, False,path_to_output=path_to_output)

        # copy functionalization_list from xyz to molfile
        cwd = os.getcwd()
        copy_functionalization_list_xyz_2_mol(os.path.join(path_to_output, first_target_file + '.xyz'),
                                              os.path.join(path_to_output, first_target_file + '.mol'))
        # Continue with looping over substituents and attaching them to the skeleton
        for idx, substituent in enumerate(substituent_list[1:]):
            logger.log(logging.INFO, f"Attaching substituent: {substituent} iteration: {idx + 1}")
            # enumerate starts at 0, so we need to add +1 to get the correct index of previous functionalization
            new_skeleton_name = skeleton_name + f'_func_{idx+1}'
            new_skeleton_path = os.path.join(path_to_output, new_skeleton_name + '.xyz')
            target_name = skeleton_name + "_func_" + str(idx + 2)
            complex = initialize_complex(skeleton_name, new_skeleton_path, substituent, path_to_database)
            complex.generate_substituent_and_write_xyz(target_name, 1.54, False, path_to_output=path_to_output)
            copy_functionalization_list_xyz_2_mol(os.path.join(path_to_output, target_name + '.xyz'),
                                                  os.path.join(path_to_output, target_name + '.mol'))
            logger.log(logging.INFO, f"Attached substituent: {substituent} iteration: {idx + 1}")
        logger.log(logging.INFO, f"Attached substituents to skeleton: {skeleton_name}")
        logger.log(logging.INFO, f"\n")

if __name__ == "__main__":

    main(path_to_database=path_to_database, path_to_substituents = path_to_substituents)
