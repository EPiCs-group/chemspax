# -*- coding: utf-8 -*-
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
import os
import sys
import ast
import glob
import numpy as np
import pandas as pd
sys.path.append("..")
from utilities import *
"""A substituent from the library can be attached to another molecule with the functions given in this file. 
The substituent is seen as a 'rigid block' that is rotated and translated. After placement of a new substituent, the 
new substituent is optimized with openbabel's FF methods. This is a constrained optimization since the skeleton 
atoms will be frozen.

This approach can be used in the functionalize_and_optimize_obabel.sh and functionalize_and_optimize_xtb.sh scripts.
"""


class Substituent:
    def __init__(self, molecule, central_atom=0, bond_length=1.2):
        """Holds information about substituents, used to write the central atom and centroid vector to the .csv database
        The information of the .csv database is then used to attach the substituents to a skeleton in
        attach_substituent.py

        :param molecule: path to substituent xyz file
        :param central_atom: index of central atom
        :param bond_length: bond length between new substituent and skeleton

        Example of writing adding a CH3 entry to the .csv database (this is done in data_preparation.py):
        >>> atom = Substituent(../CH3.xyz[36:-4], 0, 2.0)
        >>> atom.write_central_atom_and_centroid_to_csv('manually')
        """
        folder = 'substituents_xyz/manually_generated/'
        extension = '.xyz'
        self.molecule = molecule
        self.path = folder + self.molecule + extension
        self.bond_length = bond_length  # bond_length for the newly formed bond
        self.data_matrix = pd.read_table(self.path, skiprows=2, delim_whitespace=True, names=['atom', 'x', 'y', 'z'])  # read standard .xyz file
        if len(self.data_matrix) == 0:
            raise ValueError('Substituent .xyz is empty')
        self.central_atom_index = central_atom
        self.central_atom = self.data_matrix.loc[self.central_atom_index, ['x', 'y', 'z']]  # get xyz coordinate of central atom

    def scale(self, vector, central_atom):
        """Function used to scale a vector to the given input bond length

        :param vector: vector to be scaled
        :param central_atom: xyz coordinates of central atom
        :return: scaled vector
        """
        vector = vector/np.linalg.norm(vector)
        return central_atom + vector*self.bond_length  # length in angstrom

    def first_coordination(self):
        """Find centroid around a central atom, it is assumed that central atom has 3 bonds
        the atom can be symmetrical or asymmetrical, say: C-X C-Y C-Z and a free electron pair for bonding

        :return: centroid vector of substituent group
        """


        #OLD METHOD
        # amount_atoms = len(self.data_matrix.index)  # amount of rows in data_matrix
        # coordination = []
        # r_critical = 2.5  # circle radius in which we want to search (P-I has largest bond length)
        # for i in range(0, amount_atoms):
        #     current_atom = self.data_matrix.loc[i, ['x', 'y', 'z']]
        #     d = self.central_atom - current_atom
        #     # if coordinates are within radius there is a bond
        #     if ((d[0]) ** 2 + (d[1]) ** 2 + (
        #             d[2]) ** 2) < r_critical**2 and i != self.central_atom_index:
        #         coordination.append([np.linalg.norm(d), current_atom[0], current_atom[1], current_atom[2]])
        # # sort by norm(d) to find the 3 shortest bonds at top of array
        # coordination = np.array(coordination)
        # # sort in ascending order, first column is the distance
        # coordination = coordination[coordination[:, 0].argsort()]
        # # find only the xyz coordinates of shortest bonds
        # edges = coordination[0:3, 1:4]

        # find atoms bonded to central atom of substituent, use mol file since graph representation is more accurate
        edges = get_bonded_atoms(self.path[:-4]+'.mol', self.central_atom_index)
        # scale bonds such that an hypothetical symmetrical molecule is created say C-X' C-Y' C-Z'
        for i in range(np.shape(edges)[0]):
            scale_vector(self.central_atom, (edges[i, :]-self.central_atom), self.bond_length)
        # calculate centroid of this hypothetical molecule, which will be similar to real molecule
        centroid = np.sum(edges, axis=0)/edges.shape[0]  # sum over rows and divide by amount of atoms found
        # get correct orientation of total group s.t. the centroid vector is pointing towards the bond to be made
        centroid = (centroid - self.central_atom)/np.linalg.norm(centroid - self.central_atom)
        return np.array(centroid)

    def write_central_atom_and_centroid_to_csv(self, manually_or_automatically_generated):
        """Write the central atom and centroid data of the substituent to the csv database

        :param manually_or_automatically_generated: folder where substituent.xyz should be read from
        :return: new entry to csv database with molecule name, index of central atom, centroid vector for the s
        substituent
        """
        folder = 'substituents_xyz/'
        filename = 'central_atom_centroid_database.csv'
        path_to_file = folder + manually_or_automatically_generated + '_generated/' + filename
        # if there is only 1 atom to be attached there's no need to calculate a centroid, position atom == centroid
        if len(self.data_matrix) != 1:
            centroid = self.first_coordination()
        else:
            centroid = np.array(self.central_atom.values)
        write_data = pd.DataFrame([[self.molecule, int(self.central_atom_index), centroid]],
                                  columns=["group_to_be_attached", "central_atom_index", "centroid"])\
            .set_index("group_to_be_attached")

        if not os.path.exists(path_to_file):
            write_data.to_csv(path_to_file, sep=',', header=True)
        else:
            write_data.to_csv(path_to_file, sep=',', header=False, mode='a')


class Complex:
    def __init__(self, original_skeleton_name, source_data, substituent_to_be_attached, path_to_database):
        """Holds information about skeletons on which substituents will be placed,
        used to attach substituents to the skeleton using information from the .csv database.

        :param original_skeleton_name: name of the skeleton that needs to be functionalized
        :param source_data: path to skeleton xyz file
        :param substituent_to_be_attached: name of the substituent that will be attached
        :param path_to_database: path to .csv database

        Example of placing a CH3 substituent on a skeleton (this can be done in functionalize_and_optimize_obabel.sh and
        functionalize_and_optimize_xtb.sh):
        >>> some_complex = Complex('PCP-cy', '../skeletons/PCP-cy.xyz', 'CH3',
                           '../substituents_xyz/manually_generated/central_atom_centroid_database.csv')
        >>> some_complex.generate_substituent_and_write_xyz('PCP-cy_func_1', 1.54, False)
        """
        # original skeleton data (so starting point of each functionalization)
        skeleton_folder = 'skeletons/'
        extension = '.xyz'
        self.original_skeleton_path = skeleton_folder + original_skeleton_name +extension
        # skeleton data
        self.skeleton_path = source_data
        # for the first usage this is purely the skeleton, for recursive usage it's skeleton + prev. functionalization
        self.skeleton_xyz = pd.read_table(self.skeleton_path, skiprows=2, delim_whitespace=True,
                                         names=['atom', 'x', 'y', 'z'])  # read standard .xyz file
        if len(self.skeleton_xyz) == 0:
            raise ValueError('Skeleton .xyz is empty')
        # substituent data
        self.substituent_molecule = substituent_to_be_attached
        substituent_folder = 'substituents_xyz/manually_generated/'
        extension = '.xyz'
        self.substituent_path = substituent_folder + self.substituent_molecule + extension
        self.substituent_xyz = pd.read_table(self.substituent_path, skiprows=2, delim_whitespace=True,
                                         names=['atom', 'x', 'y', 'z'])  # read standard .xyz file
        if len(self.substituent_xyz) == 0:
            raise ValueError('Substituent .xyz is empty')
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

    def create_functionalization_list_all_hydrogens(self):
        # create functionalization list by finding all H atoms and atom bonded to it
        source_mol_file = self.skeleton_path
        # find all H
        search_this_atomic_num = 1
        # initalize openbabel classes
        obconversion = openbabel.OBConversion()
        # both xyz and mol can be used as input but mol contains an accurate graph representation
        if source_mol_file[-4:] == '.mol':
            obconversion.SetInFormat('mol')
        elif source_mol_file[-4:] == '.xyz':
            obconversion.SetInFormat('xyz')
        else:
            raise Exception('file type is incorrect, .mol and .xyz are supported, not', source_mol_file[-4:])

        mol = openbabel.OBMol()
        obconversion.ReadFile(mol, source_mol_file)

        functionalization_site_list = []
        # iterate over all atoms in structure
        for atom in openbabel.OBMolAtomIter(mol):
            atomic_number = atom.GetAtomicNum()
            # if a hydrogen is found
            if atomic_number == search_this_atomic_num:
                atom_to_be_functionalized_index = atom.GetIndex()  # indexing for this OB method starts at 0
                # if the hydrogen has 1 bond it can be functionalized and added to functionalization_list
                if atom.CountBondsOfOrder(1) == 1:
                    for neighbour_atom in openbabel.OBAtomAtomIter(atom):
                        bonded_atom_index = neighbour_atom.GetIndex()  # indexing for this OB method starts at 0
                        functionalization_site_list.append(
                            [int(atom_to_be_functionalized_index), int(bonded_atom_index)])

        return functionalization_site_list

    def generate_substituent_group_vector(self, length_skeleton_bonded_substituent_central=1.54):
        """Used to translate and rotate the substituent group for correct placement on the skeleton

        :param length_skeleton_bonded_substituent_central: bond length between the substituent group and skeleton
        :return: xyz matrix with correctly rotated and translated coordinates of substituent group's atoms
        """
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

    def write_connectivity_in_file(self, target_path, new_connectivity_data):
        """Used to write given connectivity data to a MDL molfile. This necessary to add correct bonding information
        from the input substituent and skeleton file to prevent weird bonds from being formed upon file conversions.
        This assumes that the MDL molfile for substituents and skeletons are correct.

        :param target_path: path to MDL molfile
        :param new_connectivity_data: dataframe of new connectivity data
        :return: same input MDL molfile but now with bonding information from given connectivity data
        """
        # convert connectivity data to ints
        new_connectivity_data = new_connectivity_data.astype(int)

        # save first part of file to write later
        lines = open(target_path).readlines()
        n_atoms = len(self.skeleton_xyz) - 1 + len(self.substituent_xyz)
        n_atoms_and_comments = n_atoms + 4

        first_part = lines[:n_atoms_and_comments]
        first_part[3] = print_mol_counts_block(first_part[3], n_atoms, len(new_connectivity_data))  # correct counts
        # read data and skip first 4 lines
        all_data = pd.read_table(target_path, skiprows=4, delim_whitespace=True, header=None)
        # fill NaN with space & save ending line to write at end of file
        end_line = all_data.fillna('').iloc[[-1]]

        # write new .mol file correctly
        os.remove(target_path)  # delete original file to prevent errors
        with open(target_path, 'w') as wr:
            wr.writelines(first_part)

        # connectivity is separated by 2 spaces (I thought), but this is not correct for the official .mol format
        # I couldn't find a way to directly write every digit in the correct format from dataframe to file
        # so now the dataframe is first written and string formatting is done afterwards
        f = open(target_path, 'ab')  # open file in binary to be able to append with np.savetxt
        np.savetxt(f, new_connectivity_data, delimiter='  ', fmt='%d')  # pd doesn't support '  ' as delimiter :(
        f.close()

        # add correct spacing of connectivity table for official .mol format, each number has 3 allocated spaces
        f = open(target_path, 'r')
        lines = f.readlines()
        connectivity_lines_only = lines[n_atoms_and_comments:]
        newlines = []
        for line in connectivity_lines_only:
            line = print_correct_connectivity_line(line)
            newlines.append(line)
        wr = open(target_path, 'w')
        lines_above_connectivity_lines = lines[:n_atoms_and_comments]
        wr.writelines(lines_above_connectivity_lines + newlines)
        wr.close()

        f = open(target_path, 'ab')  # open file in binary to be able to append with np.savetxt
        np.savetxt(f, end_line, delimiter='  ', fmt="%s")  # pd doesn't support '  ' as delimiter :(
        f.close()

    def generate_substituent_and_write_xyz(self, target_filename, length_skeleton_bonded_substituent_central=1.54,
                                           functionalize_all_hydrogens=True, use_xtb_script_after=True):
        """Used to generate translated and rotated vectors for substituents and attaching them to the skeleton

        :param target_filename: name of target xyz file that will be written
        :param length_skeleton_bonded_substituent_central: bond length between the skeleton and substituent group
        (does not really matter since force field optimization is used after placement of a new substituent)
        :param use_xtb_script_after: REDUNDANT parameter, was used to enable xtb optimization after a functionalization.
        This can still be done by uncommenting last part of this function if desired.
        :return: xyz file of substituent attached to skeleton and FF optimized MDL molfile of substituent attached to
        skeleton
        """
        folder = 'substituents_xyz/automatically_generated/'
        extension = '.xyz'
        target_path = folder + target_filename + extension

        # replace substituent x y z with newly calculated positions
        substituent_vectors = self.generate_substituent_group_vector(float(length_skeleton_bonded_substituent_central))
        substituents_new_data = self.substituent_xyz.copy()  # always copy a df to modify it!
        substituents_new_data.loc[:, ['x', 'y', 'z']] = substituent_vectors
        # remove skeleton_atom_to_be_functionalized and place substituent data at end of file
        skeleton_new_data = self.skeleton_xyz.copy()
        skeleton_new_data = skeleton_new_data.drop([self.skeleton_atom_to_be_functionalized_index])

        # old method: write substituent central atom at atom_to_be_functinoalized and paste rest of sub. to bottom
        # skeleton_new_data.loc[self.skeleton_atom_to_be_functionalized_index, :] = \
        #     substituents_new_data.loc[self.substituent_central_atom_index, :]
        # # remove central atom from dataframe of substituent group
        # substituents_new_data = substituents_new_data.drop([self.substituent_central_atom_index])

        # if we want to functionalize all hydrogens, the funtionalization_list needs to be created first
        # else just use the existing one
        if functionalize_all_hydrogens:
            self.functionalization_site_list = self.create_functionalization_list_all_hydrogens()
            # reassign atom_to_be_functionalized and bonded_atom based on new functionalization list
            if len(self.functionalization_site_list) != 0:
                # take indices from converted list and assign to correct variable
                self.skeleton_atom_to_be_functionalized_index = self.functionalization_site_list[
                    0][0]  # index in .xyz file of atom to be functionalized
                self.skeleton_bonded_atom_index = self.functionalization_site_list[
                    0][1]  # index in .xyz file of atom bonded to atom to be functionalized
                # remove first item of nested list for correct formatting later
                self.functionalization_site_list = self.functionalization_site_list[1:]
        
        # since atom_to_be_functionalized is dropped, indices in functionalization list need to shift
        # shift bonded_atom first
        self.skeleton_bonded_atom_index = self.skeleton_bonded_atom_index - 1 if self.skeleton_bonded_atom_index > self\
            .skeleton_atom_to_be_functionalized_index else self.skeleton_bonded_atom_index
        # make nested list as big as functionalization list
        new_functionalization_list = [[] for i in range(len(self.functionalization_site_list))]
        for i in range(len(self.functionalization_site_list)):
            some_list = self.functionalization_site_list[i]
            for j in range(len(some_list)):
                item = some_list[j]
                # shift is only needed if item contains index larger than atom to be functionalized
                item = item - 1 if item > self.skeleton_atom_to_be_functionalized_index else item
                some_list[j] = item
            new_functionalization_list[i] = some_list
        self.functionalization_site_list = new_functionalization_list

        # concat both dataframes and write to file
        write_data = pd.concat([skeleton_new_data, substituents_new_data])
        write_data = write_data.astype({'x': float, 'y': float, 'z': float})  # correct types in df
        write_data = write_data.set_index('atom')
        # # check if there is overlap between atoms and then ToDO: then what? Not necessary if ff_opt works
        # is_there_overlap, overlapping_atoms = check_overlap(write_data)
        # print(is_there_overlap, overlapping_atoms)
        # increase n_atoms of source file accordingly
        with open(self.skeleton_path) as f:
            n_atoms = int(f.readline())
        n_atoms += len(substituents_new_data) - 1  # atom_to_be_functionalized is dropped so 1 less atom to count
        # write to file
        with open(target_path, 'w') as wr:
            wr.write(str(n_atoms) + '\n')
            wr.write(str(self.functionalization_site_list) + '\n')
        write_data.to_csv(target_path, sep=' ', header=False, mode='a')
        # remove last whiteline generated by pandas' to_csv function
        remove_last_line(target_path)

        # fix bug with bonds being formed & weirdly broken
        # remember, indexing in .mol files starts from 1 for some reason...

        # convert skeleton to .mol if .mol file doesn't exist
        if not glob.glob(self.skeleton_path[:-4]+'.mol'):
            convert_xyz_2_mol_file(self.skeleton_path)

        # read connectivity of skeleton
        skeleton_connectivity = read_connectivity_from_mol_file(self.skeleton_path[:-4]+'.mol', len(self.skeleton_xyz))
        skeleton_connectivity = skeleton_connectivity.astype(int)
        # substituent is pasted below skeleton data so index of central atom += len(original skeleton)
        new_substituent_central_atom_index = self.substituent_central_atom_index + len(self.skeleton_xyz)

        # drop bonds with atom_to_be_functionalized
        skeleton_connectivity = skeleton_connectivity[skeleton_connectivity[0] !=
                                                      self.skeleton_atom_to_be_functionalized_index+1]
        skeleton_connectivity = skeleton_connectivity[skeleton_connectivity[1] !=
                                                      self.skeleton_atom_to_be_functionalized_index+1]

        # since atom_to_be_functionalized is dropped the indices of atoms below that need to be decreased by 1
        skeleton_connectivity[0] = skeleton_connectivity[0].apply(
            lambda x: x-1 if x > self.skeleton_atom_to_be_functionalized_index+1 else x)
        skeleton_connectivity[1] = skeleton_connectivity[1].apply(
            lambda x: x-1 if x > self.skeleton_atom_to_be_functionalized_index+1 else x)

        # append new bond between substituent_central_atom and skeleton_bonded_atom
        new_bond = pd.DataFrame({0: [new_substituent_central_atom_index], 1: [self.skeleton_bonded_atom_index+1],
                                 2: [1], 3: [0], 4: [0], 5: [0], 6: [0]})
        skeleton_connectivity = skeleton_connectivity.append(new_bond)
        skeleton_connectivity = skeleton_connectivity.astype(int)

        # convert substituent to .mol if .mol file doesn't exist
        if not glob.glob(self.substituent_path[:-4]+'.mol'):
            print('have you forgotten to run data_preparation.py? generating .mol for substituent...')
            convert_xyz_2_mol_file(self.substituent_path)

        # read connectivity of substituent
        substituent_connectivity = read_connectivity_from_mol_file(self.substituent_path[:-4]+'.mol', len(self.substituent_xyz))
        substituent_connectivity = substituent_connectivity.astype(int)
        # make indices correct, originally from substituent file indices start from 0, now: 0 + len(skeleton)
        # for initial case, for recursive case: 0 + len(skeleton) + len(substituent)
        substituent_connectivity[0] = substituent_connectivity[0] + len(skeleton_new_data)
        substituent_connectivity[1] = substituent_connectivity[1] + len(skeleton_new_data)
        # concat connectivities
        total_connectivities = pd.concat([skeleton_connectivity, substituent_connectivity])

        # convert functionalized skeleton to .mol
        # ToDo: C-O 13 and 14 are triple bonded, manually fix or does xtb/gaussian do that?
        convert_xyz_2_mol_file(target_path)
        self.write_connectivity_in_file(target_path[:-4]+'.mol', total_connectivities)
        # optimize .mol file
        # indices_to_freeze = None
        # since self.skeleton changes because of recursive functionalizations, the original skeleton atoms need to be
        # freezed before ff optimization to try to improve calculation efficiency

        # in functionalize_and_optimize scripts name after uff is skeleton + _func_i
        # and in xtb skeleton + _func + _i + _opt
        n_iteration = int(target_filename.split('_')[-1])

        n_atoms_original_skeleton = int(open(self.original_skeleton_path).readline()) - n_iteration
        try:
            indices_to_freeze = [i for i in range(n_atoms_original_skeleton)]
        except:
            print('No indices found to freeze')
            indices_to_freeze = None
        ff_optimize(target_path[:-4]+'.mol', 'gaff', indices_to_freeze)
        ff_optimize(target_path[:-4]+'.mol', 'uff', indices_to_freeze)

        # conversion from .mol to .xyz is taken care of in xtb bash script: xtbopt.mol --> xtbopt.xyz
        # REDUNDANT: xyz conversion can be done in batch after all functionalizations are done
        # set to True if the xtb bash script will be used
        # if not use_xtb_script_after:
        #     # convert .mol file back to xyz file
        #     convert_mol_2_xyz_file(target_path[:-4]+'.mol')
        #     # remove last white line
              # remove_last_line(target_path)


if __name__ == "__main__":
    # os.remove('../substituents_xyz/manually_generated/central_atom_centroid_database.csv')
    # for file in glob.glob('../substituents_xyz/manually_generated/*.xyz'):
    #     print(file)
    #     atom = Substituent(file[39:-4], 0, 2.0)
    #     atom.write_central_atom_and_centroid_to_csv('manually')
    # phenyl = Substituent('C6H6', 0, 2.0)
    # print(phenyl.first_coordination())
    # phenyl.write_central_atom_and_centroid_to_csv('manually')
    folder_name = '../substituents_xyz/automatically_generated/'
    folder = os.listdir('../substituents_xyz/automatically_generated/')
    for item in folder:
        if item.endswith(".xyz"):
            os.remove(os.path.join(folder_name, item))
        elif item.endswith(".mol"):
            os.remove(os.path.join(folder_name, item))

    some_complex = Complex('PCP-cy', '../skeletons/PCP-cy.xyz', 'CH3',
                           '../substituents_xyz/manually_generated/central_atom_centroid_database.csv')
    some_complex.generate_substituent_and_write_xyz('PCP-cy_func_1', 1.54, False)
    # some_complex.write_connectivity_in_file('../substituents_xyz/automatically_generated/something.mol', 'moh')
    other_complex = Complex('PCP-cy', '../substituents_xyz/automatically_generated/PCP-cy_func_1.xyz', 'CH3',
                            '../substituents_xyz/manually_generated/central_atom_centroid_database.csv')
    other_complex.generate_substituent_and_write_xyz('PCP-cy_func_2', 1.54, False)
    some_other_complex = Complex('PCP-cy', '../substituents_xyz/automatically_generated/PCP-cy_func_2.xyz', 'CH3',
                            '../substituents_xyz/manually_generated/central_atom_centroid_database.csv')
    some_other_complex.generate_substituent_and_write_xyz('PCP-cy_func_3', 1.54, False)
    # other_other_complex = Complex('../substituents_xyz/automatically_generated/something_2.xyz', 'CH4N2OH',
    #                         '../substituents_xyz/manually_generated/central_atom_centroid_database.csv')
    # other_other_complex.generate_substituent_and_write_xyz('something_3', 1.54)
