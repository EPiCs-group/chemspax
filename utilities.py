# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
import ase.io as io
from ase.visualize import view
import ase.build
import numpy as np
import pandas as pd
from openbabel import openbabel


def distance(a, b):
    d = a - b
    return np.sqrt((d[0]) ** 2 + (d[1]) ** 2 + (d[2]) ** 2)


def visualize_xyz_file(filename, save_picture=False, manually_generated=True):
    """Visualize an .xyz file, number_atoms at top of .xyz needs to be correct!

    :param filename:
    :param save_picture:
    :param manually_generated:
    :return: visualization of .xyz file
    """
    molecule = io.read(filename)
    view(molecule)
    if save_picture:  # ToDo: use matplotlib for nicer visualizations
        if manually_generated:
            io.write('substituents_xyz/visualizations/' + filename[35:-4] + '.png', molecule, rotation='45x,45y,0z')
        else:
            io.write('substituents_xyz/visualizations/' + filename[40:-4] + '.png', molecule, rotation='45x,45y,0z')
    else:
        return


def read_central_atom_index(filename):
    """Read second line of a .xyz file to get the central atom index

    :param filename:
    :return: index of central atom
    """
    with open(filename) as f:
        next(f)
        return int(next(f))


def find_distance(filename, index_of_atom_1, index_of_atom_2):
    """Find distance between two atoms based on their index in .xyz file

    :param filename:
    :param index_of_atom_1:
    :param index_of_atom_2:
    :return: distance between atoms
    """
    molecule = io.read(filename)
    return molecule.get_distance(index_of_atom_1, index_of_atom_2)


def remove_last_line(filename):
    with open(filename) as f:
        lines = f.readlines()
        last = len(lines) - 1
        lines[last] = lines[last].replace('\r', '').replace('\n', '')
    with open(filename, 'w') as wr:
        wr.writelines(lines)


def create_molecule_and_write_xyz(input_molecule, filename):
    """
    https://wiki.fysik.dtu.dk/ase/ase/build/build.html?highlight=ase%20build%20molecule#ase.build.molecule
    :param input_molecule: 
    :param filename: 
    :return: 
    """
    molecule = ase.build.molecule(input_molecule)
    molecule.write(filename, 'xyz')
    remove_last_line(filename)


def scale_vector(starting_point, vector, length):
    """ Scales a vector with a given length

    :param starting_point:
    :param vector:
    :param length:
    :return: scaled vector
    """
    vector = vector/np.linalg.norm(vector)
    return starting_point + vector*length


def convert_list_of_string_to_np_array(array_string):
    """Pandas is importing np arrays as strings, use this converter to convert the list of 1 string to a np.array
    https://stackoverflow.com/questions/42755214/how-to-keep-numpy-array-when-saving-pandas-dataframe-to-csv
    :param array_string: example: ['[-0.33332174004836124 0.9428131403470853 0.0]'] dtype=list
    :return: np.array(interpretation_of_array_string)
    """
    array_string = str(array_string).replace('[', '').replace(']', '').replace(' ', ', ').replace("'", '').split(', ')
    return np.array([float(x) for x in list(array_string)])


def generate_random_rotation_matrix():
    """https://math.stackexchange.com/questions/442418/random-generation-of-rotation-matrices

    :return: uniformly random rotation matrix
    """
    q, r = np.linalg.qr(np.random.rand(3, 3))
    return q


def read_connectivity_from_mol_file(source_file, n_atoms):
    # ToDO: delimiter is not always a space if index1=56 and index2=103 there is no space
    # https://chem.libretexts.org/Courses/University_of_Arkansas_Little_Rock/ChemInformatics_(2017)%3A_Chem_4399%2F%2F5399/2.2%3A_Chemical_Representations_on_Computer%3A_Part_II/2.2.2%3A_Anatomy_of_a_MOL_file
    skip_rows = n_atoms + 4  # title line, whiteline, comment line, n_atoms & n_bonds line = 4 lines to skip
    connectivity = pd.read_table(source_file, skiprows=skip_rows, delim_whitespace=True, header=None)

    # print(connectivity.loc[connectivity[0] == 56103]) #= int(connectivity[0])[3:]
    # drop last 'M END' line
    connectivity = connectivity.drop([len(connectivity) - 1])

    connectivity = connectivity.fillna(0)
    connectivity = connectivity.astype(int)

    # if there is no space between idx1 and idx2 the numbers still need to be separated
    # this happens if idx2 > 99, for example: 56 103 --> 56103 or if idx1 and idx2 are > 99 for ex: 100 100 --> 100100
    for i in range(len(connectivity)):
        current_row = connectivity.loc[i, [0]]
        if current_row[0] > 1100:  # 1 100 is the lowest set of integers for which this problem will occur
            row_string = str(current_row[0]).strip()
            if len(row_string) == 5:
                connectivity.loc[i, [0, 1, 2]] = int(row_string[:2]), int(row_string[2:]), int(connectivity.loc[i, [2]].
                                                                                               values)
            elif len(row_string) == 6:
                connectivity.loc[i, [0, 1, 2]] = int(row_string[:3]), int(row_string[3:]), int(connectivity.loc[i, [2]].
                                                                                               values)
    return connectivity


def convert_xyz_2_mol_file(source_file):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "mol")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, source_file)
    target_filename = source_file[:-4] + '.mol'
    obConversion.WriteFile(mol, target_filename)


def convert_mol_2_xyz_file(source_file):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "xyz")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, source_file)
    target_filename = source_file[:-4] + '.xyz'
    obConversion.WriteFile(mol, target_filename)


def print_mol_counts_block(old_string, n_atoms, n_bonds):
    n_atoms = str(n_atoms)
    n_bonds = str(n_bonds)
    static_part = old_string[6:]

    # the first 6 places are reserved for n_atoms and n_bonds, these need to be formattted correctly
    if len(n_atoms) == 1:
        new_string = '  ' + n_atoms
    elif len(n_atoms) == 2:
        new_string = ' ' + n_atoms
    else:
        new_string = n_atoms

    if len(n_bonds) == 1:
        new_string += '  ' + n_bonds
    elif len(n_bonds) == 2:
        new_string += ' ' + n_bonds
    else:
        new_string += n_bonds
    # the rest of the string stays the same and can be added back
    new_string += static_part

    return new_string


def print_correct_connectivity_line(line):
    line_list = line.split('  ')
    to_return = ['  ' for i in range(16)]  # initialize empty list with 2 spaces as separator

    # 6 cases, either idx 1 is 1 10 100 or idx is 1 10 100
    atom_index_1 = line_list[0]
    if len(atom_index_1) == 3:
        to_return[0] = atom_index_1[0]
        to_return[1] = atom_index_1[1]
        to_return[2] = atom_index_1[2]
    elif len(atom_index_1) == 2:
        to_return[0] = ' '
        to_return[1] = atom_index_1[0]
        to_return[2] = atom_index_1[1]
    elif len(atom_index_1) == 1:
        to_return[0] = ' '
        to_return[1] = ' '
        to_return[2] = atom_index_1

    atom_index_2 = line_list[1]
    if len(atom_index_2) == 3:
        to_return[3] = atom_index_2[0]
        to_return[4] = atom_index_2[1]
        to_return[5] = atom_index_2[2]
    elif len(atom_index_2) == 2:
        to_return[3] = ' '
        to_return[4] = atom_index_2[0]
        to_return[5] = atom_index_2[1]
    elif len(atom_index_2) == 1:
        to_return[3] = ' '
        to_return[4] = ' '
        to_return[5] = atom_index_2

    # the rest of the line is static with 2 spaces as separator
    to_return[7] = line_list[2]
    to_return[9] = line_list[3]
    to_return[11] = line_list[4]
    to_return[13] = line_list[5]
    to_return[15] = line_list[6]

    to_return = ''.join([str(elem) for elem in to_return])  # turn list into string
    return to_return


if __name__ == '__main__':
    # molec = 'H2O'
    # create_molecule_and_write_xyz('H2O', 'substituents_xyz/automatically_generated/' + molec + '.xyz')
    # visualize_xyz_file('substituents_xyz/automatically_generated/something.xyz', True, False)
    # visualize_xyz_file('skeletons/RuPNP_iPr_skl.xyz', False, False)
    # print(convert_list_of_string_to_np_array(['[-0.33332174004836124 0.9428131403470853 0.0]']))
    # print(read_central_atom_index('substituents_xyz/automatically_generated/CH4.xyz'))
    # print(find_distance('substituents_xyz/automatically_generated/CH4.xyz', 2, 3)==1.7473026804689453)
    # print(read_connectivity_from_mol_file('random.mol', 98))
    # convert_mol_2_xyz_file('random.mol')
    # convert_xyz_2_mol_file('substituents_xyz/automatically_generated/something_2.xyz')
    # print(print_mol_counts_block(15, 15, 0))
    print_correct_connectivity_line('120  113  1  0  0  0  0')
