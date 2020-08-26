# -*- coding: utf-8 -*-
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
import numpy as np
import pandas as pd
'''Take ligands and find centroids'''


class Substituent:
    # ask user for central atom
    def __init__(self, path_to_data='substituents_xyz/CH3.xyz', bond_length=1.2):
        self.path = path_to_data
        self.bond_length = bond_length  # bond_length for the newly formed bond
        self.data_matrix = pd.read_table(self.path, skiprows=2, delim_whitespace=True, names=['atom', 'x', 'y', 'z'])  # read standard .xyz file
        self.central_atom = self.data_matrix.loc[0, ['x', 'y', 'z']]  # get xyz coordinate of central atom

    @staticmethod
    def distance(a, b):
        d = a - b
        return np.sqrt((d[0]) ** 2 + (d[1]) ** 2 + (d[2]) ** 2)

    def scale(self, vector, central_atom):
        return central_atom + vector*self.bond_length  # length in angstrom

    def first_coordination(self):
        indices = []
        centroid = np.zeros((3, 3))
        r_critical = 2.0  # circle radius in which we want to search
        for i in range(1, len(self.central_atom)+1):
            temp = self.data_matrix.loc[i, ['x', 'y', 'z']]
            d = self.central_atom - temp
            if ((d[0]) ** 2 + (d[1]) ** 2 + (
                    d[2]) ** 2) < r_critical**2:  # if coordinates are within radius append
                centroid[i-1, :] = np.array(temp[0:4])
                indices.append(i)
        return centroid, indices

    def check_vector(self):
        centroid, indices = self.first_coordination()
        vector = centroid - np.array(self.central_atom)  # diff in distance = vector
        assert_list = []
        for i in range(len(vector)):
            norm = np.linalg.norm(np.array(vector[i]))
            vector[i] = vector[i] / norm  # vector is now normalized
            vector_euclidian_norm = vector[i, 0]**2+vector[i,1]**2+vector[i,2]**2
            assert_list.append(vector_euclidian_norm)
        assert_list = [True for i, j in zip(assert_list, [1, 1, 1]) if np.floor(i) == j]  # check if every value == 1.0
        return vector, True if all(assert_list) else False

    def generate_substituent_vectors(self):
        centroid, indices = self.first_coordination()
        vector = centroid - np.array(self.central_atom)
        new_xyz_matrix = []
        for i in range(0, 3):
            norm = np.linalg.norm(np.array(vector[i]))
            vector[i] = vector[i] / norm
            scaled_vector = self.scale(vector[i], self.central_atom)
            new_xyz_matrix.append(scaled_vector)
        return new_xyz_matrix


if __name__ == "__main__":
    methyl = Substituent('substituents_xyz/CH3.xyz', 1.2)
    # print(methyl.data_matrix)
    # print(methyl.central_atom)
    # centroid , indices = methyl.first_coordination()
    # vector, assertion = methyl.check_vector()
    # print(vector)
    # print(assertion)
    # test_scale = methyl.scale(vector[0], methyl.central_atom, 1.2)
    # print(methyl.distance(methyl.central_atom, test_scale))
    c = methyl.generate_substituent_vectors()
    c = (c[0] + c[1] + c[2])/3
    print(c)
    print(methyl.distance(c, methyl.central_atom))  # should be approx 0.4 in correct equilateral triangle
