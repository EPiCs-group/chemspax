# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
import ase.io as io
from ase.visualize import view
import ase.build
import numpy as np


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


if __name__ == '__main__':
    # molec = 'H2O'
    # create_molecule_and_write_xyz('H2O', 'substituents_xyz/automatically_generated/' + molec + '.xyz')
    # visualize_xyz_file('substituents_xyz/automatically_generated/something.xyz', True, False)
    visualize_xyz_file('skeletons/RuPNP_iPr_skl.xyz', False, False)
    # print(convert_list_of_string_to_np_array(['[-0.33332174004836124 0.9428131403470853 0.0]']))
    # print(read_central_atom_index('substituents_xyz/automatically_generated/CH4.xyz'))
    # print(find_distance('substituents_xyz/automatically_generated/CH4.xyz', 2, 3)==1.7473026804689453)
