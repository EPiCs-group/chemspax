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

from chemspax.attach_substituent import Complex
from chemspax.utilities import copy_functionalization_list_xyz_2_mol
"""Main class that takes command line arguments and passes them to attach_substituent.py's Complex class
"""


def initialize_complex(original_skeleton_name, source_data, substituent_name, path_to_database):
    """Main function that takes command line arguments and passes them to attach_substituent.py's Complex class
    """
    complex = Complex(original_skeleton_name, source_data, substituent_name, path_to_database)
    return complex


def main(skeleton_list, substituent_list):
    # initialize logger
    logger = logging.getLogger(__name__)
    logger.log(logging.INFO, "Looks like you want to add some substituents to a metal-ligand complex & optimize 'em.\n")
    logger.log(logging.INFO, "Created by: Adarsh Kalikadien & Vivek Sinha")

    # main loop to attach substituents to the skeleton
    for skeleton in skeleton_list:
        skeleton_path = os.path.normpath(skeleton)
        skeleton_name = os.path.basename(skeleton_path).split(".")[0]
        print("Attaching substituents to skeleton: " + skeleton_name)
        target_name = skeleton_name + "_func"
        # intialize the complex class with a random substituent to be able to find lenght of functionalization_list
        first_target_file = target_name + '_1'
        path_to_substituent_database = os.path.join("substituents_xyz", "manually_generated", "central_atom_centroid_database.csv")
        complex = initialize_complex(skeleton_name, skeleton_path, "CH3", path_to_substituent_database)
        if substituent_list is not None:
            print(f"initializing complex by placing {substituent_list[0]} on {skeleton_name}")
            complex = initialize_complex(skeleton_name, skeleton_path, substituent_list[0], path_to_substituent_database)
            complex.generate_substituent_and_write_xyz(first_target_file, 1.54, False)
        else:
            print("No substituent list found, using random substituents based on length of functionalization_list")
            # create substituent_list with correct length (number of substituents)
            len_functionalization_list = len(complex.functionalization_site_list)
            # select random substituents, same length as functionalization_list (number of functionalization sites)
            all_substituents = glob.glob("substituents_xyz/manually_generated/*.xyz")
            # only get name of substituent file without extension
            all_substituents = [os.path.basename(os.path.normpath(substituent)).split(".")[0] for substituent in all_substituents]
            substituent_list = random.choices(all_substituents, k=len_functionalization_list)
            print(f"initializing complex by placing {substituent_list[0]} on {skeleton_name}")
            complex = initialize_complex(skeleton_name, skeleton_path, substituent_list[0], path_to_substituent_database)
            complex.generate_substituent_and_write_xyz(first_target_file, 1.54, False)

        # # copy functionalization_list from xyz to molfile
        # copy_functionalization_list_xyz_2_mol(target_name + '_1.xyz', target_name + '_1.mol')
        # # Continue with looping over substituents and attaching them to the skeleton
        # for idx, substituent in enumerate(substituent_list[1:]):
        #     print("Attaching substituent: " + substituent)
        #     # enumerate starts at 0, so we need to add +2 to get the correct index of previous functionalization
        #     complex = initialize_complex()
        #     complex.generate_substituent_and_write_xyz()
        #     print("Attached substituent: " + substituent)
        # print("Attached substituents to skeleton: " + skeleton_name)


if __name__ == "__main__":


    parser = argparse.ArgumentParser(prog='chemspax', description='Attach substituents to a skeleton molecule')
    parser.add_argument('-s', '--substituent_list', help='List of substituents to attach', required=False, default=None)

    args = parser.parse_args()
    substituent_list = args.substituent_list



    list_of_skeleton_files = glob.glob("skeletons/*.xyz")
    # substituent_list = ["CH3", "CH3"]
    main(list_of_skeleton_files, substituent_list)

#
# def main():
#     original_skeleton_name = sys.argv[1]
#     source_data = sys.argv[2]
#     target_file = sys.argv[3]
#     substituent_molecule = sys.argv[4]
#     path_to_database = sys.argv[5]
#     length_skeleton_bonded_atom_substituent_central_atom = sys.argv[6]
#     use_xtb_script_after = True if sys.argv[7].lower() == 'true' else False
#
#     some_complex = Complex(original_skeleton_name, source_data, substituent_molecule, path_to_database)
#     some_complex.generate_substituent_and_write_xyz(target_file, length_skeleton_bonded_atom_substituent_central_atom,
#                                                     use_xtb_script_after)
#
#
# if __name__ == '__main__':
#     main()
