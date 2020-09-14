# -*- coding: utf-8 -*-
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
import os
import sys
import ast
import numpy as np
import pandas as pd
sys.path.append("..")
from utilities import *

'''class Substituent: Take ligands and find centroids that point in correct direction to be added to a skeleton
   class Complex: Functionalize a skeleton with a substituent using the generated .csv database 
'''


class Substituent:
    def __init__(self, molecule, central_atom=0, bond_length=1.2):
        folder = '../substituents_xyz/manually_generated/'
        extension = '.xyz'
        self.molecule = molecule
        self.path = folder + self.molecule + extension
        self.bond_length = bond_length  # bond_length for the newly formed bond
        self.data_matrix = pd.read_table(self.path, skiprows=2, delim_whitespace=True, names=['atom', 'x', 'y', 'z'])  # read standard .xyz file
        self.central_atom_index = central_atom
        self.central_atom = self.data_matrix.loc[self.central_atom_index, ['x', 'y', 'z']]  # get xyz coordinate of central atom

    @staticmethod
    def distance(a, b):
        d = a - b
        return np.sqrt((d[0]) ** 2 + (d[1]) ** 2 + (d[2]) ** 2)

    def scale(self, vector, central_atom):
        vector = vector/np.linalg.norm(vector)
        return central_atom + vector*self.bond_length  # length in angstrom

    def first_coordination(self):
        # find centroid around a central atom, it is assumed that central atom has 3 bonds
        # the atom can be symmetrical or asymmetrical, say: C-X C-Y C-Z and a free electron pair for bonding
        amount_atoms = len(self.data_matrix.index)  # amount of rows in data_matrix
        coordination = []
        r_critical = 2.5  # circle radius in which we want to search (P-I has largest bond length)
        for i in range(0, amount_atoms):
            current_atom = self.data_matrix.loc[i, ['x', 'y', 'z']]
            d = self.central_atom - current_atom
            # if coordinates are within radius there is a bond
            if ((d[0]) ** 2 + (d[1]) ** 2 + (
                    d[2]) ** 2) < r_critical**2 and i != self.central_atom_index:
                coordination.append([np.linalg.norm(d), current_atom[0], current_atom[1], current_atom[2]])
        # sort by norm(d) to find the 3 shortest bonds at bottom of array
        coordination = np.array(coordination)
        coordination = coordination[coordination[:, -1].argsort()]
        edges = coordination[-3:, 1:4]
        # scale bonds such that an hypothetical symmetrical molecule is created say C-X' C-Y' C-Z'
        for i in range(np.shape(edges)[0]):
            scale_vector(self.central_atom, (edges[i, :]-self.central_atom), self.bond_length)
        # calculate centroid of this hypothetical molecule, which will be similar to real molecule
        centroid = (edges[0, :]+edges[1, :]+edges[2, :])/3
        # get correct orientation of total group s.t. the centroid vector is pointing towards the bond to be made
        centroid = (centroid - self.central_atom)/np.linalg.norm(centroid - self.central_atom)
        return np.array(centroid)

    def write_central_atom_and_centroid_to_csv(self, manually_or_automatically_generated):
        folder = '../substituents_xyz/'
        filename = 'central_atom_centroid_database.csv'
        path_to_file = folder + manually_or_automatically_generated + '_generated/' + filename
        centroid = self.first_coordination()
        write_data = pd.DataFrame([[self.molecule, self.central_atom_index, centroid]],
                                  columns=["group_to_be_attached", "central_atom_index", "centroid"])\
            .set_index("group_to_be_attached")

        if not os.path.exists(path_to_file):
            write_data.to_csv(path_to_file, sep=',', header=True)
        else:
            write_data.to_csv(path_to_file, sep=',', header=False, mode='a')


class Complex:
    def __init__(self, source_data, substituent_to_be_attached, path_to_database):
        # skeleton data
        self.skeleton_path = source_data
        self.skeleton_xyz = pd.read_table(self.skeleton_path, skiprows=2, delim_whitespace=True,
                                         names=['atom', 'x', 'y', 'z'])  # read standard .xyz file
        # substituent data
        self.substituent_molecule = substituent_to_be_attached
        substituent_folder = '../substituents_xyz/manually_generated/'
        extension = '.xyz'
        self.substituent_path = substituent_folder + self.substituent_molecule + extension
        self.substituent_xyz = pd.read_table(self.substituent_path, skiprows=2, delim_whitespace=True,
                                         names=['atom', 'x', 'y', 'z'])  # read standard .xyz file
        self.database_df = pd.read_csv(path_to_database, delimiter=',')
                                       # ,converters={'centroid': convert_list_of_string_to_np_array})
        try:
            self.substituent_df = self.database_df.loc[self.database_df['group_to_be_attached']
                                                       == self.substituent_molecule]
        except KeyError:
            raise KeyError
        self.substituent_central_atom_index = int(self.substituent_df['central_atom_index'])
        self.substituent_centroid_vector = self.substituent_df['centroid'].values
        self.substituent_centroid_vector = convert_list_of_string_to_np_array(self.substituent_centroid_vector)
        self.substituent_central_atom_xyz = self.substituent_xyz.loc[
            self.substituent_central_atom_index, ['x', 'y', 'z']]

        # functionalization list from source file
        with open(self.skeleton_path) as f:
            lines = f.readlines()
            self.functionalization_site_list = lines[1]
            # print(self.functionalization_site_list)
            # f.close()
        # convert list from string to integer
        self.functionalization_site_list = ast.literal_eval(self.functionalization_site_list)
        if len(self.functionalization_site_list) != 0:
            # take indices from converted list and assign to correct variable
            self.skeleton_atom_to_be_functionalized_index = self.functionalization_site_list[
                0][0]  # index in .xyz file of atom to be functionalized
            self.skeleton_bonded_atom_index = self.functionalization_site_list[
                0][1]  # index in .xyz file of atom bonded to atom to be functionalized
            # remove first item of nested list for correct formatting later
            self.functionalization_site_list = self.functionalization_site_list[1:]
            # write to .xyz file in generate_and_write_xyz function
        else:
            print('No more indices left. Exiting program')
            sys.exit()

        self.skeleton_atom_to_be_functionalized_xyz = self.skeleton_xyz.loc[
            self.skeleton_atom_to_be_functionalized_index, [
                'x', 'y', 'z']]  # get xyz coordinate of atom to be functionalized - H
        self.skeleton_bonded_atom_xyz = self.skeleton_xyz.loc[
            self.skeleton_bonded_atom_index, [
                'x', 'y', 'z']]  # get xyz coordinate of bonded atom - C (in CH3: C= central atom)
        self.bond_length = self.skeleton_atom_to_be_functionalized_xyz \
                         - self.skeleton_bonded_atom_xyz  # vector with origin on C and points to H in xyz plane
        self.normalized_bond_vector = self.bond_length / np.linalg.norm(
            self.bond_length)  # real bond between C-H in xyz plane

    def generate_substituent_group_vector(self, length_skeleton_bonded_substituent_central=1.54):
        normal_vector = self.substituent_centroid_vector
        # normal_vector = normal_vector / np.linalg.norm(normal_vector)  # vector is already unit vector (redundant)

        # construct rotation matrix
        bond_length_norm = np.array(self.normalized_bond_vector.astype('float64'))

        # https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
        v = np.cross(normal_vector.T, bond_length_norm.T)  # v is perpendicular to normal vector and bond between C-H
        v_x = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        v_xsq = np.dot(v_x, v_x)
        c = np.dot(bond_length_norm.T, normal_vector.T)
        if c != -1.0:
            rotation_matrix = np.eye(3) + v_x + v_xsq * (1 / (1 + c))
        else:
            rotation_matrix = np.eye(3)
            # raise RotationMatrixError

        n_atoms, n_columns = len(self.substituent_xyz), len(self.substituent_xyz.columns)  # atoms in substituent group
        substituent_vectors = np.array(self.substituent_xyz.loc[:, ['x', 'y', 'z']].values)  # xyz only
        # calculate new position of substituent_central_atom, the other atoms will be placed around this
        new_position_substiuent = np.array(scale_vector(self.skeleton_bonded_atom_xyz,
                            (self.skeleton_atom_to_be_functionalized_xyz - self.skeleton_bonded_atom_xyz),
                            float(length_skeleton_bonded_substituent_central)))

        # do rotation first
        for i in range(n_atoms):
            substituent_vectors[i, :] = np.dot(rotation_matrix, substituent_vectors[i, :].T)
        # save copy of correctly rotated central atom of substituent
        substituent_central_atom = substituent_vectors[self.substituent_central_atom_index, :].copy()
        # do translation after
        for i in range(n_atoms):
            # ToDo: for O-CH3 do new_position_substituent - Oxygen, how to detect these cases?
            substituent_vectors[i, :] = substituent_vectors[i, :] + (new_position_substiuent - substituent_central_atom)
        return substituent_vectors

    def generate_substituent_and_write_xyz(self, target_filename, length_skeleton_bonded_substituent_central=1.54):
        folder = '../substituents_xyz/automatically_generated/'
        extension = '.xyz'
        target_path = folder + target_filename + extension

        # replace substituent x y z with newly calculated positions
        substituent_vectors = self.generate_substituent_group_vector(float(length_skeleton_bonded_substituent_central))
        substituents_new_data = self.substituent_xyz.copy()  # always copy a df to modify it!
        substituents_new_data.loc[:, ['x', 'y', 'z']] = substituent_vectors
        # replace skeleton_atom_to_be_functionalized with central atom of substituent group
        skeleton_new_data = self.skeleton_xyz.copy()
        skeleton_new_data.loc[self.skeleton_atom_to_be_functionalized_index, :] = \
            substituents_new_data.loc[self.substituent_central_atom_index, :]
        # remove central atom from dataframe of substituent group
        substituents_new_data = substituents_new_data.drop([self.substituent_central_atom_index])
        # concat both dataframes
        write_data = pd.concat([skeleton_new_data, substituents_new_data])
        write_data = write_data.astype({'x': float, 'y': float, 'z': float})  # correct types in df
        write_data = write_data.set_index('atom')
        # increase n_atoms of source file accordingly
        with open(self.skeleton_path) as f:
            n_atoms = int(f.readline())
        n_atoms += len(substituents_new_data)
        # write to file
        with open(target_path, 'w') as wr:
            wr.write(str(n_atoms)+'\n')
            wr.write(str(self.functionalization_site_list)+'\n')

        write_data.to_csv(target_path, sep=' ', header=False, mode='a')
        # remove last whiteline generated by pandas' to_csv function
        remove_last_line(target_path)


if __name__ == "__main__":
    # os.remove('substituents_xyz/manually_generated/central_atom_centroid_database.csv')
    # for file in glob.glob('substituents_xyz/manually_generated/*.xyz'):
    #     atom = Substituent(file[36:-4], 0, 2.0)
    #     atom.write_central_atom_and_centroid_to_csv('manually')
    # ethyl has central atom index=4 and needs to be done separately
    # ethyl = Substituent('CH2CH3', 4, 2.0)
    # ethyl.write_central_atom_and_centroid_to_csv('manually')
    if os.path.exists('../substituents_xyz/automatically_generated/something.xyz'):
        os.remove('../substituents_xyz/automatically_generated/something.xyz')
    some_complex = Complex('../skeletons/RuPNP_iPr_skl.xyz', 'CCl2F',
                           '../substituents_xyz/manually_generated/central_atom_centroid_database.csv')
    some_complex.generate_substituent_and_write_xyz('something', 1.54)
    other_complex = Complex('../substituents_xyz/automatically_generated/something.xyz', 'CCl2F',
                           '../substituents_xyz/manually_generated/central_atom_centroid_database.csv')
    other_complex.generate_substituent_and_write_xyz('something_1', 1.54)
