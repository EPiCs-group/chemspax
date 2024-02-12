# -*- coding: utf-8 -*-
#                                                     #
#  __author__ = Adarsh Kalikadien                     #
#  __institution__ = TU Delft                         #
#  __contact__ = a.v.kalikadien@tudelft.nl            #
import os
import glob
from chemspax.main import main

current_directory = os.getcwd()


path_to_substituents = os.path.join(current_directory, "substituents_xyz", "manually_generated/")  # should always point to the chemspax/substituents_xyz/manually_generated folder
path_to_database = os.path.join(path_to_substituents, "central_atom_centroid_database.csv")

substituent_list = ["CH3", "CH3", "CH3", "CH3"]
skeleton_list = glob.glob(os.path.join(current_directory, "skeletons_test", "*.xyz"))
path_to_skeletons = os.path.join(current_directory, "skeletons_test")
working_directory = current_directory
path_to_output = os.path.join(current_directory, "substituents_xyz", "automatically_generated")
main(skeleton_list, substituent_list, path_to_database, path_to_substituents, path_to_skeletons, working_directory, path_to_output)

