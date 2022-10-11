# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
import glob
import os
import logging

from utilities import *
from chemspax.attach_substituent import Substituent
"""Run this script before running any of the main files and bash scripts to prepare the csv database and .mol files 
for attach_substituent.py 
"""


def convert_file(from_extension, to_extension, folder):
    for file_path in glob.glob(os.path.join(folder, f'*.{from_extension}')):
        # conversion is only necessary if output file doesn't exist
        if not glob.glob(file_path[:-4] + f'.{to_extension}'):
            if to_extension == 'mol':
                convert_xyz_2_mol_file(file_path)
            elif to_extension == 'xyz':
                convert_mol_2_xyz_file(file_path)
            else:
                raise SystemExit('File extension is not supported for conversion')


def prepare_data(path_to_substituents, path_to_skeletons, path_to_database):
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

    substituent_folder = path_to_substituents
    skeleton_folder = path_to_skeletons
    csv_database_file = path_to_database

    # ToDo: add logging
    logger.log(logging.INFO, 'Converting substituent files to .xyz format')
    # convert substituents .mol file to .xyz files
    convert_file('mol', 'xyz', substituent_folder)
    logger.log(logging.INFO, 'Converting substituent files to .mol format')
    # convert substituents .xyz files to .mol files
    convert_file('xyz', 'mol', substituent_folder)

    logger.log(logging.INFO, 'Converting skeleton files to .xyz format')
    # convert substituents .mol file to .xyz files
    convert_file('mol', 'xyz', skeleton_folder)
    logger.log(logging.INFO, 'Converting skeleton files to .mol format')
    # convert skeleton .xyz files to .mol files
    convert_file('xyz', 'mol', skeleton_folder)

    # check if csv database exists and delete if that's the case
    if glob.glob(csv_database_file):
        logger.log(logging.INFO, 'Deleting existing csv database')
        os.remove(csv_database_file)
    #
    # create csv database for all substituents
    logger.log(logging.INFO, 'Creating csv database for substituents')
    for file_path in glob.glob(substituent_folder + '*.xyz'):
        # this assumes that the central atom of the substituent is the first atom in the file!
        molecule_name = os.path.basename(os.path.normpath(file_path)).split(".")[0]
        atom = Substituent(molecule=molecule_name, path_to_substituents=substituent_folder, central_atom=0, bond_length=2.0)
        atom.write_central_atom_and_centroid_to_csv()


if __name__ == '__main__':
    # test functionality
    current_directory = os.getcwd()
    path_to_substituents = os.path.join(current_directory, "substituents_xyz", "manually_generated/")
    path_to_database = os.path.join(path_to_substituents, "central_atom_centroid_database.csv")
    path_to_skeletons = os.path.join(current_directory, "skeletons")

    prepare_data(path_to_substituents, path_to_skeletons, path_to_database)
