# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
"""A different approach to the same problem: The use case would be given a C-H bond to be functionalized:
generate 3 atoms A,B,C such that H,A,B,C make a tetrahedron
then H can be replaced/substituted by say another atoms Z (practical case would be C,N,P)
in the end this code will also allow substitution by O-CH3 type groups
"""
import numpy as np
import pandas as pd
from exceptions import *
from utilities import *
import os
import ast
import sys


class Complex:
    def __init__(self, path_to_source_data, functionalization_site_list=None, recursive_or_intial='recursive'):
        self.path = path_to_source_data
        self.functionalization_site_list = functionalization_site_list
        self.recursive_or_initial = recursive_or_intial
        self.data_matrix = pd.read_table(self.path, skiprows=2, delim_whitespace=True,
                                         names=['atom', 'x', 'y', 'z'])  # read standard .xyz file

        if self.recursive_or_initial.lower() == 'initial':
            self.functionalization_site_list = ast.literal_eval(self.functionalization_site_list)
            self.atom_to_be_functionalized_index = self.functionalization_site_list[0][0]  # index in .xyz file of atom to be functionalized
            self.bonded_atom_index = self.functionalization_site_list[0][1]  # index in .xyz file of atom bonded to atom to be functionalized
        elif self.recursive_or_initial.lower() == 'recursive' or functionalization_site_list is None:
            # read second line of .xyz file
            with open(path_to_source_data) as f:
                lines = f.readlines()
                self.functionalization_site_list = lines[1]
                # print(self.functionalization_site_list)
                # f.close()
            # convert list from string to integer
            self.functionalization_site_list = ast.literal_eval(self.functionalization_site_list)
            if len(self.functionalization_site_list) != 0:
                # take indices from converted list and assign to correct variable
                self.atom_to_be_functionalized_index = self.functionalization_site_list[0][0]  # index in .xyz file of atom to be functionalized
                self.bonded_atom_index = self.functionalization_site_list[0][1]  # index in .xyz file of atom bonded to atom to be functionalized
                # remove first item of nested list for correct formatting in generate_and_write_xyz function
                self.functionalization_site_list = self.functionalization_site_list[1:]
                # write to .xyz file in generate_and_write_xyz function
            else:
                print('No more indices left. Exiting program')
                sys.exit()
        else:
            raise InvalidRecursiveOrInitialArgumentError

        self.atom_to_be_functionalized_xyz = self.data_matrix.loc[
            self.atom_to_be_functionalized_index, ['x', 'y', 'z']]  # get xyz coordinate of atom to be functionalized - H
        self.bonded_atom_xyz = self.data_matrix.loc[
            self.bonded_atom_index, ['x', 'y', 'z']]  # get xyz coordinate of bonded atom - C (in CH3: C= central atom)
        self.bond_length = self.atom_to_be_functionalized_xyz - self.bonded_atom_xyz  # vector with origin on C and points to H in xyz plane
        self.normalized_bond_vector = self.bond_length / np.linalg.norm(self.bond_length)  # real bond between C-H in xyz plane
        self.equilateral_triangle = np.array([[0, 1/np.sqrt(3.0), 0],
                                              [-0.5, -0.5/np.sqrt(3.0), 0],
                                              [0.5, -0.5/np.sqrt(3.0), 0]])  # equilateral triangle with centroid at origin

    @staticmethod
    def distance(a, b):
        d = a - b
        return np.sqrt((d[0])**2 + (d[1])**2 + (d[2])**2)

    @staticmethod
    def find_angle(v1, v2):
        cos_th = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return np.arccos(cos_th) * 57.29577951308232

    def find_centroid(self):
        # find new centroid and find where equilateral triangle needs to be translated to
        b = np.linalg.norm(self.bond_length)  # bond to be functionalized -H
        b = b * (2.0 * np.sqrt(2.0 / 3.0))
        self.equilateral_triangle = b*self.equilateral_triangle  # make side lengths equal to tetrahedral bond lengthw
        centroid = self.atom_to_be_functionalized_xyz + (b/3.0) * self.normalized_bond_vector
        return centroid

    def generate_substituent_vectors(self):
        centroid = self.find_centroid()
        normal_vector = np.array([0, 0, 1])
        normal_vector = normal_vector / np.linalg.norm(normal_vector)  # make unit vector
        starting_coordinate = np.zeros(3)  # origin of original defined equilateral triangle
        # theta = np.arccos(np.dot(self.bond_length_norm.T, normal_vector.T))

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

        num_rows, num_columns = self.equilateral_triangle.shape
        substituent_vectors = np.zeros((num_rows, num_columns))
        # find path in correct direction & shift object to new position + rotate vectors s.t. they are aligned with C-H
        for i in range(num_rows):
            substituent_vectors[i, :] = (centroid - starting_coordinate) + \
                                        np.dot(rotation_matrix, self.equilateral_triangle[i, :].T)
        return substituent_vectors[0], substituent_vectors[1], substituent_vectors[2]

    def generate_substituent_and_write_xyz(self, filename, substituent_0_atom, substituent_1_atom, substituent_2_atom, substituent_3_atom,
                                           view_file=False, save_file=False):
        folder = 'substituents_xyz/automatically_generated/'
        extension = '.xyz'
        filename = folder + filename + extension
        v_1, v_2, v_3 = self.generate_substituent_vectors()

        # scale vectors with right bond length between atom_to_be_functionalized and new substituent
        v_1 = self.atom_to_be_functionalized_xyz + 1.1*(v_1-self.atom_to_be_functionalized_xyz)/(np.linalg.norm(v_1-self.atom_to_be_functionalized_xyz))
        v_2 = self.atom_to_be_functionalized_xyz + 1.1*(v_2-self.atom_to_be_functionalized_xyz)/(np.linalg.norm(self.atom_to_be_functionalized_xyz-v_2))
        v_3 = self.atom_to_be_functionalized_xyz + 1.1 * (v_3-self.atom_to_be_functionalized_xyz ) / (
            np.linalg.norm(self.atom_to_be_functionalized_xyz - v_3))

        atom_to_be_functionalized = self.data_matrix.loc[
            self.atom_to_be_functionalized_index, ['atom', 'x', 'y', 'z']]  # type = pandas.Series convert to df
        atom_to_be_functionalized = pd.DataFrame([[atom_to_be_functionalized[0], atom_to_be_functionalized[1],
                                                   atom_to_be_functionalized[2], atom_to_be_functionalized[3]]],
                                                 columns=['atom', 'x', 'y', 'z']).set_index('atom')
        bonded_atom = self.data_matrix.loc[
            self.bonded_atom_index, ['atom', 'x', 'y', 'z']]
        bonded_atom = pd.DataFrame([[bonded_atom[0], bonded_atom[1], bonded_atom[2], bonded_atom[3]]],
                                   columns=['atom', 'x', 'y', 'z']).set_index('atom')

        substituent_1 = pd.DataFrame([[substituent_1_atom.upper() if len(substituent_1_atom) == 1 else
                                       substituent_1_atom, v_1[0], v_1[1], v_1[2]]],
                                     columns=['atom', 'x', 'y', 'z']).set_index('atom')
        substituent_2 = pd.DataFrame([[substituent_2_atom.upper() if len(substituent_2_atom) == 1 else
                                       substituent_2_atom, v_2[0], v_2[1], v_2[2]]],
                                     columns=['atom', 'x', 'y', 'z']).set_index('atom')
        substituent_3 = pd.DataFrame([[substituent_3_atom.upper() if len(substituent_3_atom) == 1 else
                                       substituent_3_atom, v_3[0], v_3[1], v_3[2]]],
                                     columns=['atom', 'x', 'y', 'z']).set_index('atom')

        write_data = pd.concat([bonded_atom, atom_to_be_functionalized, substituent_1, substituent_2, substituent_3])
        write_data = write_data.astype({'x': float, 'y': float, 'z': float})  # correct types in df

        # write correct .xyz format
        if self.recursive_or_initial.lower() == 'initial':
            n_atoms = 5
            with open(filename, 'a') as file:
                file.write(str(n_atoms) + '\n')
                file.write(str(self.functionalization_site_list[1:]) + '\n')  # write list for recursive functionalizations
            file.close()
            write_data.to_csv(filename, sep=' ', header=False, mode='a')
            # change atom to be functionalized with substituent 0
            with open(filename) as f:
                lines = f.readlines()
                substituent_0_line = list(lines[self.atom_to_be_functionalized_index + 2])
                substituent_0_line[0] = substituent_0_atom
                lines[self.atom_to_be_functionalized_index + 2] = ''.join(
                    substituent_0_line)  # don't forget comment line and num_atom line
                with open(filename, 'w') as wr:
                    # write to target filename
                    wr.writelines(lines)
                wr.close()
            # remove last whiteline generated by pandas' to_csv function
            remove_last_line(filename)
        elif self.recursive_or_initial.lower() == 'recursive':
            with open(self.path) as f:
                # read source file and modify
                lines = f.readlines()
                lines[0] = int(lines[0]) + 3  # increase n_atoms correctly
                lines[0] = str(lines[0])+'\n'
                lines[1] = str(self.functionalization_site_list)+'\n'  # write modified functionalization list to file
                # replace atom_to_be_functionalized with substituent 0 (str type can't be assigned a value)
                substituent_0_line = list(lines[self.atom_to_be_functionalized_index+2])
                substituent_0_line[0] = substituent_0_atom
                lines[self.atom_to_be_functionalized_index+2] = ''.join(substituent_0_line) # don't forget comment line and num_atom line
                lines[-1] = lines[-1] + '\n'  # add newline to prevent .to_csv to write on same line
            with open(filename, 'w') as wr:
                # write to target filename
                wr.writelines(lines)
            wr.close()
            # remove atom_to_be_functionalized and bonded_atom from dataframe since they are already in .xyz file
            write_data = pd.concat([substituent_1, substituent_2, substituent_3])
            write_data = write_data.astype({'x': float, 'y': float, 'z': float})  # correct types in df

            write_data.to_csv(filename, sep=' ', header=False, mode='a')
            # remove last whiteline generated by pandas' to_csv function
            remove_last_line(filename)
        else:
            raise InvalidRecursiveOrInitialArgumentError

        # check whether file must be saved or viewed only
        if 'manually_generated' in filename and view_file and save_file:
            visualize_xyz_file(filename, save_file, True)
        elif 'manually_generated' in filename and view_file and not save_file:
            visualize_xyz_file(filename, save_file, True)
        elif 'automatically_generated' in filename and view_file and save_file:
            visualize_xyz_file(filename, save_file, False)
        elif 'automatically_generated' in filename and view_file and not save_file:
            visualize_xyz_file(filename, save_file, False)
        else:  # write the file only, no viewer
            return


if __name__ == '__main__':
    # example usage for ethyl
    # generates a CH3-CH3 'molecule' (orientation and bonds are of course incorrect)
    if os.path.exists('substituents_xyz/automatically_generated/CH4.xyz'):
        os.remove('substituents_xyz/automatically_generated/CH4.xyz')
        functionalization_list = '[[1, 0], [4, 0]]'
        methyl = Complex('substituents_xyz/manually_generated/CH3.xyz', functionalization_list, 'initial')
        methyl.generate_substituent_and_write_xyz('CH4', 'C', 'H', 'H', 'C', False, False)
        ethyl = Complex('substituents_xyz/automatically_generated/CH4.xyz', None, 'recursive')
        ethyl.generate_substituent_and_write_xyz('CH4', 'C', 'H', 'H', 'H', True, False)
    functionalization_list = '[[1, 0], [4, 0]]'
    methyl = Complex('substituents_xyz/manually_generated/CH3.xyz', functionalization_list, 'initial')
    methyl.generate_substituent_and_write_xyz('CH4', 'C', 'H', 'H', 'C', False, False)
    ethyl = Complex('substituents_xyz/automatically_generated/CH4.xyz', None, 'recursive')
    ethyl.generate_substituent_and_write_xyz('CH4', 'C', 'H', 'H', 'H', True, False)


