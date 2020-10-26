# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
from utilities import *
from attach_substituent_folder.attach_substituent import Substituent
import glob
import os
"""Run this script before running any of the main files and bash scripts to prepare the csv database and .mol files 
for attach_substituent.py 
"""


def prepare_data():
    substituent_folder = 'substituents_xyz/manually_generated/'
    skeleton_folder = 'skeletons/'
    csv_database_file = substituent_folder + 'central_atom_centroid_database.csv'

    # check if csv database exists and delete if that's the case
    if glob.glob(csv_database_file):
        os.remove(csv_database_file)

    # ToDo: compress these 4 functions into one
    # convert substituents .xyz files to .mol files
    for file in glob.glob(substituent_folder + '*.xyz'):
        # conversion is only necessary if the .mol file doesn't exist
        if not glob.glob(file[:-4]+'.mol'):
            convert_xyz_2_mol_file(file)

    # convert skeleton .xyz files to .mol files
    for file in glob.glob(skeleton_folder + '*.xyz'):
        # conversion is only necessary if the .mol file doesn't exist
        if not glob.glob(file[:-4]+'.mol'):
            convert_xyz_2_mol_file(file)

    # convert substituents .mol file to .xyz files
    for file in glob.glob(substituent_folder + '*.mol'):
        # conversion is only necessary if the .xyz file doesn't exist
        if not glob.glob(file[:-4]+'.xyz'):
            convert_mol_2_xyz_file(file)

    # convert substituents .mol file to .xyz files
    for file in glob.glob(skeleton_folder + '*.mol'):
        # conversion is only necessary if the .xyz file doesn't exist
        if not glob.glob(file[:-4] + '.xyz'):
            convert_mol_2_xyz_file(file)

    # create csv database for all substituents
    for file in glob.glob(substituent_folder + '*.xyz'):
        # this assumes that the central atom of the substituent is the first atom in the file!
        atom = Substituent(file[36:-4], 0, 2.0)
        atom.write_central_atom_and_centroid_to_csv('manually')


if __name__ == '__main__':
    prepare_data()
