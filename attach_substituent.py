# -*- coding: utf-8 -*-
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
import numpy as np
import pandas as pd
from utilities import *
import os
import glob
'''Take ligands and find centroids that point in correct direction to be added to a skeleton'''


class Substituent:
    # ask user for central atom
    def __init__(self, molecule, central_atom=0, bond_length=1.2):
        folder = 'substituents_xyz/manually_generated/'
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
        edges[0, :] = self.scale(edges[0, :]-self.central_atom, self.central_atom)
        edges[1, :] = self.scale(edges[1, :]-self.central_atom, self.central_atom)
        edges[2, :] = self.scale(edges[2, :]-self.central_atom, self.central_atom)
        # calculate centroid of this hypothetical molecule, which will be similar to real molecule
        centroid = (edges[0, :]+edges[1, :]+edges[2, :])/3
        # get correct orientation of total group s.t. the centroid vector is pointing towards the bond to be made
        centroid = (centroid - self.central_atom)/np.linalg.norm(centroid - self.central_atom)
        return np.array(centroid)

    def write_central_atom_and_centroid_to_csv(self, manually_or_automatically_generated):
        folder = 'substituents_xyz/'
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


if __name__ == "__main__":
    os.remove('substituents_xyz/manually_generated/central_atom_centroid_database.csv')
    for file in glob.glob('substituents_xyz/manually_generated/*.xyz'):
        atom = Substituent(file[36:-4], 0, 2.0)
        atom.write_central_atom_and_centroid_to_csv('manually')
    # ethyl has central atom index=4 and needs to be done separately
    # ethyl = Substituent('CH2CH3', 4, 2.0)
    # ethyl.write_central_atom_and_centroid_to_csv('manually')

    # centroid = methyl.first_coordination()
    # print(centroid)
