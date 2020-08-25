# -*- coding: utf-8 -*-
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
import numpy as np
import pandas as pd


class Substituent:
    def __init__(self, path_to_data):
        self.path = path_to_data
        self.centroid = np.zeros((3, 3))
        self.data_matrix = pd.read_table(self.path, skiprows=2, delim_whitespace=True, names=['atom', 'x', 'y', 'z'])  # read standard .xyz file
        self.central_atom = self.data_matrix.loc[0, ['x', 'y', 'z']]  # get xyz coordinate of central atom

    @staticmethod
    def scale(vector, central_atom, length):
        return central_atom + vector*length  # length in angstrom

    @staticmethod
    def distance(a, b):
        return np.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)

    def first_coordination(self):
        indices = []
        r_critical = 2.0  # circle radius in which we want to search
        for i in range(1, len(self.central_atom)+1):
            temp = self.data_matrix.loc[i, ['x', 'y', 'z']]
            if ((self.central_atom[0] - temp[0]) ** 2 + (self.central_atom[1] - temp[1]) ** 2 + (
                    self.central_atom[2] - temp[2]) ** 2) < r_critical**2:  # if coordinates are within radius append
                self.centroid[i-1, :] = np.array(temp[0:4])
                indices.append(i)
        return self.centroid, indices

    def check_vector(self):
        centroid, indices = self.first_coordination()
        vector = centroid - np.array(self.central_atom)  # diff in distance = vector
        assert_list = []
        for i in range(0, 3):
            norm = np.linalg.norm(np.array(vector[i]))
            vector[i] = vector[i] / norm  # vector is now normalized
            vector_euclidian_norm = vector[i, 0]**2+vector[i,1]**2+vector[i,2]**2
            assert_list.append(vector_euclidian_norm)
        assert_list = [True for i, j in zip(assert_list, [1, 1, 1]) if np.floor(i) == j]  # check if every value == 1.0
        return vector, True if all(assert_list) else False


if __name__ == "__main__":
    methyl = Substituent('substituents_xyz/methyl.xyz')
    # print(methyl.data_matrix)
    # print(methyl.central_atom)
    # centroid , indices = methyl.first_coordination()
    vector, assertion = methyl.check_vector()
    print(vector)
    print(assertion)
    test_scale = methyl.scale(vector[0], methyl.central_atom, 1.2)
    print(methyl.distance(methyl.central_atom, test_scale))