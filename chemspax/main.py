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
from pathlib import Path

from chemspax.attach_substituent import Complex
from chemspax.utilities import copy_functionalization_list_xyz_2_mol
from chemspax.data_preparation import prepare_data


def initialize_complex(original_skeleton_name, source_data, substituent_name, path_to_database, path_to_skeletons, path_to_substituents):
    """Main function that takes command line arguments and passes them to attach_substituent.py's Complex class
    """
    structure = Complex(original_skeleton_name, source_data, substituent_name, path_to_database, path_to_skeletons, path_to_substituents)
    return structure


def main(skeleton_list, substituent_list, path_to_database, path_to_substituents, path_to_skeletons, working_directory, path_to_output):
    # initialize logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    # create a file handler
    handler = logging.FileHandler('chemspax.log', mode='a')
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

    Path(path_to_output).mkdir(parents=True, exist_ok=True)  # create directory for functionalized output

    # check if working directory is set in environment variable, else set to current folder
    # if working_directory is None:
    #     logger.log(logging.INFO, "No ChemSpaX home directory specified, using default")
    #     working_directory = os.getcwd()
    # else:
    #     logger.log(logging.INFO, f"ChemSpaX working directory is set in environment variable to {working_directory}")
    logger.log(logging.INFO, f"ChemSpaX working directory is set to {working_directory}")

    # start ChemSpaX main loop to attach substituents to the skeleton
    logger.log(logging.INFO, "Looks like you want to add some substituents to a metal-ligand complex & optimize 'em.\n")
    logger.log(logging.INFO, "Created by: Adarsh Kalikadien & Vivek Sinha\n")

    if not skeleton_list:
        logger.log(logging.INFO, "No skeletons supplied, exiting")
        exit(1)
    else:
        logger.log(logging.INFO, f"{len(skeleton_list)} skeletons provided to functionalize")

    # list of all substituents for later use
    all_substituents = glob.glob(os.path.join(path_to_substituents, "C*.xyz"))
    # only get name of substituent file without extension
    all_substituents = [os.path.basename(os.path.normpath(substituent)).split(".")[0] for substituent in
                        all_substituents]

    for skeleton in skeleton_list:
        skeleton_path = os.path.normpath(skeleton)
        skeleton_name = os.path.basename(skeleton_path).split(".")[0]
        logger.log(logging.INFO, f"Attaching substituents to skeleton: {skeleton_name}")
        # intialize the complex class with a random substituent to be able to find lenght of functionalization_list
        first_target_file = skeleton_name + "_func_" + '1'
#        path_to_substituent_database = os.path.join("substituents_xyz", "manually_generated", "central_atom_centroid_database.csv")    #to keep things consistent let us just use "path_to_database" as used in "Complex" class as well"
        complex = initialize_complex(skeleton_name, skeleton_path, "CH3", path_to_database, path_to_skeletons, path_to_substituents)
        if substituent_list is not None:
            logger.log(logging.INFO, f"initializing complex by placing {substituent_list[0]} on {skeleton_name}")
        else:
            logger.log(logging.INFO, f"initializing complex by placing random substituents based on length of functionalization_list")
            # create substituent_list with correct length (number of substituents)
            len_functionalization_list = len(complex.functionalization_site_list)
            # select random carbon (C) substituents, same length as functionalization_list (number of functionalization sites)
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
            # add random substituents to the list
            for i in range(len(substituent_list), len(complex.functionalization_site_list)):
                substituent_list.append(random.choice(all_substituents))

        complex = initialize_complex(skeleton_name, skeleton_path, substituent_list[0], path_to_database, path_to_skeletons, path_to_substituents)
        complex.generate_substituent_and_write_xyz(target_filename=first_target_file, path_to_output=path_to_output,
                                                   length_skeleton_bonded_substituent_central=1.54,
                                                   use_xtb_script_after=False)

        # copy functionalization_list from xyz to molfile
        copy_functionalization_list_xyz_2_mol(os.path.join(path_to_output, first_target_file + '.xyz'),
                                              os.path.join(path_to_output, first_target_file + '.mol'))
        # Continue with looping over substituents and attaching them to the skeleton
        for idx, substituent in enumerate(substituent_list[1:]):
            logger.log(logging.INFO, f"Attaching substituent: {substituent} iteration: {idx + 1}")
            # enumerate starts at 0, so we need to add +1 to get the correct index of previous functionalization
            new_skeleton_name = skeleton_name + f'_func_{idx+1}'
            new_skeleton_path = os.path.join(path_to_output, new_skeleton_name + '.xyz')
            target_name = skeleton_name + "_func_" + str(idx + 2)
            complex = initialize_complex(skeleton_name, new_skeleton_path, substituent, path_to_database, path_to_skeletons, path_to_substituents)
            complex.generate_substituent_and_write_xyz(target_filename=target_name, path_to_output=path_to_output,
                                                       length_skeleton_bonded_substituent_central=1.54,
                                                       use_xtb_script_after=False)
            copy_functionalization_list_xyz_2_mol(os.path.join(path_to_output, target_name + '.xyz'),
                                                  os.path.join(path_to_output, target_name + '.mol'))
            logger.log(logging.INFO, f"Attached substituent: {substituent} iteration: {idx + 1}")
        logger.log(logging.INFO, f"Attached substituents to skeleton: {skeleton_name}")
        logger.log(logging.INFO, f"\n")


if __name__ == "__main__":
    # get current directory since this is where substituents and substituents CSV file is located
    current_directory = os.getcwd()
    try:
        chemspax_working_directory = os.environ["CHEMSPAX_HOME_DIR"]
    except KeyError:
        chemspax_working_directory = current_directory

    path_to_substituents = os.path.join(current_directory, "substituents_xyz", "manually_generated/")
    path_to_database = os.path.join(path_to_substituents, "central_atom_centroid_database.csv")
    path_to_output = os.path.join(current_directory, "substituents_xyz", "automatically_generated")
    path_to_skeletons = os.path.join(current_directory, "skeletons")

    # parse input arguments
    parser = argparse.ArgumentParser(prog='chemspax', description='Attach substituents to a skeleton molecule')
    parser.add_argument('-s', '--substituent', help='Name of substituent to attach', required=False, action='append',
                        default=None)
    parser.add_argument('-k', '--skeleton', help='Provide folder where skeletons are located', required=False,
                        action='store', default=None)
    parser.add_argument('-o', '--output', help="path to output files", required=False, action='store',
                        default=None)
    # if --prepare argument is there, the data_preparation script will run prior to the main file
    # users should still check their structures after this preparation is done
    parser.add_argument('-pd', '--prepare', help="do xyz and mol conversions, and add new substituents to substituent database", required=False, action='store_true', default=False)
    args = parser.parse_args()
    if args.substituent is not None:
        # make list of supplied substituents
        substituent_list = []
        for arg in args.substituent:
            substituent_list.append(arg)
    else:
        # if None the substituents will be chosen automatically in the main loop
        substituent_list = None
    if args.skeleton is not None:
        path_to_skeletons = args.skeleton
        list_of_skeleton_files = glob.glob(os.path.join(args.skeleton, "*.xyz"))
    else:
        list_of_skeleton_files = glob.glob(os.path.join(chemspax_working_directory, "skeletons", "*.xyz"))
    if args.output is not None:
        path_to_output = args.output
    if args.prepare:
        prepare_data(path_to_substituents=path_to_substituents, path_to_skeletons=path_to_skeletons, path_to_database=path_to_database)
    else:
        main(skeleton_list=list_of_skeleton_files, substituent_list=substituent_list, path_to_database=path_to_database, path_to_substituents=path_to_substituents,
             path_to_skeletons=path_to_skeletons, working_directory=chemspax_working_directory, path_to_output=path_to_output)
