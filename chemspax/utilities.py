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
from openbabel import pybel

"""This file contains utilities used in attach_substituent.py
"""

# reduce warnings from openbabel
ob_log_handler = openbabel.OBMessageHandler()
ob_log_handler.SetOutputLevel(0)
openbabel.obErrorLog.SetOutputLevel(0)


def get_mol(source_mol_file):
    """Setting file format, preparing molecule from source file

    :param source_mol_file:
    """
    # initalize openbabel classes
    obconversion = openbabel.OBConversion()
    # both xyz and mol can be used as input but mol contains an accurate graph representation
    if source_mol_file[-4:] == ".mol":
        obconversion.SetInFormat("mol")
    elif source_mol_file[-4:] == ".xyz":
        obconversion.SetInFormat("xyz")
    else:
        raise Exception(
            "file type is incorrect, .mol and .xyz are supported, not",
            source_mol_file[-4:],
        )
    mol = openbabel.OBMol()
    obconversion.ReadFile(mol, source_mol_file)

    return mol, obconversion


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
            io.write(
                "substituents_xyz/visualizations/" + filename[35:-4] + ".png",
                molecule,
                rotation="45x,45y,0z",
            )
        else:
            io.write(
                "substituents_xyz/visualizations/" + filename[40:-4] + ".png",
                molecule,
                rotation="45x,45y,0z",
            )
    else:
        return


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
        lines[last] = lines[last].replace("\r", "").replace("\n", "")
    with open(filename, "w") as wr:
        wr.writelines(lines)


def create_molecule_and_write_xyz(input_molecule, filename):
    """
    https://wiki.fysik.dtu.dk/ase/ase/build/build.html?highlight=ase%20build%20molecule#ase.build.molecule
    :param input_molecule:
    :param filename:
    :return:
    """
    molecule = ase.build.molecule(input_molecule)
    molecule.write(filename, "xyz")
    remove_last_line(filename)


def scale_vector(starting_point, vector, length):
    """Scales a vector with a given length

    :param starting_point:
    :param vector:
    :param length:
    :return: scaled vector
    """
    vector = vector / np.linalg.norm(vector)
    return starting_point + vector * length


def get_bonded_atoms(source_mol_file, atom_index):
    """Use openbabel's methods to find the coordinates of all atoms that are bonded to a given atom

    :param source_mol_file:
    :param atom_index:
    :return: numpy array of atoms bonded to a given atom
    """
    # make molecule from source file with openbabel class
    mol, _ = get_mol(source_mol_file)

    # make atom object for atom we want
    atom = mol.GetAtom(
        atom_index + 1
    )  # for obmol get functions indexing starts from 1

    coordinates_list = []
    for neighbour_atom in openbabel.OBAtomAtomIter(atom):
        # there's also a .GetCoordinate() function but it returns a SWIG object...
        coordinates = [
            neighbour_atom.x(),
            neighbour_atom.y(),
            neighbour_atom.z(),
        ]
        coordinates_list.append(coordinates)

    coordinates_array = np.array(
        [np.array(coords) for coords in coordinates_list]
    )

    return coordinates_array


def get_neighbour_bond_distance(
    source_mol_file, atom_index, search_this_atomic_num
):
    """Use Openbabel's meethods to find the distance of a specific atom type to a given atom, the atoms should be bonded
    and the specific atom type that is searched for should be searched by using it's atomic number.
    This function is useful if the given atom has a known index in the molecule
    if you don't know the index of the central atom use get_neighbour_distance_search.py or
    get_neighbour_bond_distance_search()

    :param source_mol_file:
    :param atom_index:
    :param search_this_atomic_num:
    :return: list of distance between indicated central atom and atoms that are bonded to it with given atomic number
    """
    atom_index = int(atom_index)
    search_this_atomic_num = int(search_this_atomic_num)

    # make molecule from source file with openbabel class
    mol, _ = get_mol(source_mol_file)

    # make atom object for atom of whose bonds we want to look at
    atom = mol.GetAtom(
        atom_index + 1
    )  # for obmol get functions indexing starts from 1

    distance_list = []
    for neighbour_atom in openbabel.OBAtomAtomIter(atom):
        atomic_number = neighbour_atom.GetAtomicNum()
        # change this atomic number to change the atom you're looking for
        if atomic_number == search_this_atomic_num:
            # if this is the atom we're looking for, append the distance
            distance = neighbour_atom.GetDistance(atom)
            distance_list.append(distance)

    return distance_list


def get_neighbour_bond_distance_search(
    source_mol_file, central_atom_atomic_num, search_this_atomic_num
):
    """Use Openbabel's meethods to find the distance of a specific atom type to a given atom type, the atoms should be
    bonded and the specific atom type that is searched for should be searched by using it's atomic number.
    This function is useful if the given atom does not have a known index in the molecule,
    if you know the index of the central atom use get_neighbour_distance.py or get_neighbour_bond_distance()

    :param source_mol_file:
    :param central_atom_atomic_num:
    :param search_this_atomic_num:
    :return: list of distance between a central atom with unknown index and atoms that are bonded to it
    with given atomic number
    """
    central_atom_atomic_num = int(central_atom_atomic_num)
    search_this_atomic_num = int(search_this_atomic_num)

    # make molecule from source file with openbabel class
    mol = get_mol(source_mol_file)

    distance_list = []
    # iterate over atoms in molecule to search the central atom
    for atom in openbabel.OBMolAtomIter(mol):
        atomic_number = atom.GetAtomicNum()
        if atomic_number == central_atom_atomic_num:
            for neighbour_atom in openbabel.OBAtomAtomIter(atom):
                atomic_number = neighbour_atom.GetAtomicNum()
                if atomic_number == search_this_atomic_num:
                    distance = neighbour_atom.GetDistance(atom)
                    distance_list.append(distance)

    return distance_list


def set_neighbour_bond_distance_search(
    source_mol_file,
    central_atom_atomic_num,
    search_this_atomic_num,
    new_length,
    new_atomic_num,
):
    """Use Openbabel's meethods to scale the distance of a specific atom type to a given atom type, the atoms should be
    bonded and the specific atom type that is searched for should be searched by using it's atomic number.
    For example, if you have a Mn pincer complex where Br is bonded to Mn and you want to replace Br with H
    set_neighbour_bond_distance_search(filename, 25, 35, 1.68, 1)

    :param source_mol_file:
    :param central_atom_atomic_num:
    :param search_this_atomic_num:
    :param new_length
    :param new_atomic_num
    :return: Same xyz or mol file with atom replaced
    """
    central_atom_atomic_num = int(central_atom_atomic_num)
    search_this_atomic_num = int(search_this_atomic_num)

    # make molecule from source file with openbabel class
    mol, obconversion = get_mol(source_mol_file)

    # iterate over atoms in molecule to search the central atom
    for atom in openbabel.OBMolAtomIter(mol):
        atomic_number = atom.GetAtomicNum()
        if atomic_number == central_atom_atomic_num:  # central atom is found
            central_atom_vector = np.array([atom.x(), atom.y(), atom.z()])
            for neighbour_atom in openbabel.OBAtomAtomIter(atom):
                atomic_number = neighbour_atom.GetAtomicNum()
                if (
                    atomic_number == search_this_atomic_num
                ):  # atom that you were searching for was found
                    search_this_atom_vector = np.array(
                        [
                            neighbour_atom.x(),
                            neighbour_atom.y(),
                            neighbour_atom.z(),
                        ]
                    )
                    d = neighbour_atom.GetDistance(atom)
                    distance_vector = (
                        search_this_atom_vector - central_atom_vector
                    ) / d
                    distance_vector = distance_vector / np.linalg.norm(
                        distance_vector
                    )
                    new_vector = (
                        central_atom_vector + distance_vector * new_length
                    )
                    neighbour_atom.SetAtomicNum(new_atomic_num)
                    neighbour_atom.SetVector(
                        new_vector[0], new_vector[1], new_vector[2]
                    )
                    obconversion.WriteFile(mol, source_mol_file)


def convert_list_of_string_to_np_array(array_string):
    """Pandas is importing np arrays as strings, use this converter to convert the list of 1 string to a np.array
    https://stackoverflow.com/questions/42755214/how-to-keep-numpy-array-when-saving-pandas-dataframe-to-csv
    :param array_string: example: ['[-0.33332174004836124 0.9428131403470853 0.0]'] dtype=list
    :return: np.array(interpretation_of_array_string)
    """
    array_string = (
        str(array_string)
        .replace("[", "")
        .replace("]", "")
        .replace(" ", ", ")
        .replace("'", "")
        .split(", ")
    )
    return np.array([float(x) for x in list(array_string)])


def generate_random_rotation_matrix():
    """https://math.stackexchange.com/questions/442418/random-generation-of-rotation-matrices

    :return: uniformly random rotation matrix
    """
    q, r = np.linalg.qr(np.random.rand(3, 3))
    return q


def read_connectivity_from_mol_file(source_file, n_atoms):
    """Reads connectivity from a .mol file, each number in .mol file has 3 allocated spaces and the file looks like:
    idx1 idx2 bond_type bond_stereochemistry 0 0 0
    :param source_file:
    :param n_atoms:
    :return: connectivity block from .mol file except 'M  END' line
    """
    # https://chem.libretexts.org/Courses/University_of_Arkansas_Little_Rock/ChemInformatics_(2017)%3A_Chem_4399%2F%2F5399/2.2%3A_Chemical_Representations_on_Computer%3A_Part_II/2.2.2%3A_Anatomy_of_a_MOL_file
    skip_rows = (
        n_atoms + 4
    )  # title line, whiteline, comment line, n_atoms & n_bonds line = 4 lines to skip
    connectivity_data = open(source_file).readlines()[skip_rows:-1]

    connectivity_list = []
    # each item has 3 spaces, based on this the item can be loaded in a dataframe
    for line in connectivity_data:
        idx1 = int(line[:3].replace(" ", ""))
        idx2 = int(line[3:6].replace(" ", ""))
        bond = int(line[6:9].replace(" ", ""))
        stereochem = int(line[9:12].replace(" ", ""))
        other_info = int(line[12:15].replace(" ", ""))
        other_info2 = int(line[15:18].replace(" ", ""))
        other_info3 = int(line[18:21].replace(" ", ""))
        # print([idx1, idx2, bond, stereochem, other_info, other_info2, other_info3])
        connectivity_list.append(
            [
                idx1,
                idx2,
                bond,
                stereochem,
                other_info,
                other_info2,
                other_info3,
            ]
        )
    connectivity = pd.DataFrame(
        connectivity_list, columns=[0, 1, 2, 3, 4, 5, 6]
    )

    # old approach
    # connectivity = pd.read_table(source_file, skiprows=skip_rows, delim_whitespace=True, header=None)
    #
    # # drop last 'M END' line
    # connectivity = connectivity.drop([len(connectivity) - 1])
    #
    # connectivity = connectivity.fillna(0)
    # connectivity = connectivity.astype(int)
    # #
    # # if there is no space between idx1 and idx2 the numbers still need to be separated
    # # this happens if idx2 > 99, for example: 56 103 --> 56103 or if idx1 and idx2 are > 99 for ex: 100 100 --> 100100
    # for i in range(len(connectivity)):
    #     current_row = connectivity.loc[i, [0]]
    #     if current_row[0] > 1100:  # 1 100 is the lowest set of integers for which this problem will occur
    #         row_string = str(current_row[0]).strip()
    #         if len(row_string) == 5:
    #             connectivity.loc[i, [0, 1, 2]] = int(row_string[:2]), int(row_string[2:]), int(connectivity.loc[i, [1]].
    #                                                                                            values)
    #         elif len(row_string) == 6:
    #             connectivity.loc[i, [0, 1, 2]] = int(row_string[:3]), int(row_string[3:]), int(connectivity.loc[i, [1]].
    #                                                                                            values)
    return connectivity


def convert_xyz_2_mol_file(source_file):
    """Converts .xyz to .mol file with the same filename

    :param source_file:
    :return: .mol file with same filename as .xyz file
    """
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "mol")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, source_file)
    target_filename = source_file[:-4] + ".mol"
    obConversion.WriteFile(mol, target_filename)


def convert_mol_2_xyz_file(source_file):
    """Converts .mol to .xyz file with the same filename

    :param source_file:
    :return: .xyz file with same filename as .mol file
    """
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "xyz")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, source_file)
    target_filename = source_file[:-4] + ".xyz"
    obConversion.WriteFile(mol, target_filename)


def print_mol_counts_block(old_string, n_atoms, n_bonds):
    """The counts block in a .mol file is corrected based on bonds and atoms given

    :param old_string:
    :param n_atoms:
    :param n_bonds:
    :return: corrected counts block of .mol file
    """
    n_atoms = int(n_atoms)
    n_bonds = int(n_bonds)

    if (n_bonds or n_atoms) > 999:
        raise ValueError("n_atoms or n_bonds can not have more than 3 digits")

    n_atoms = str(n_atoms)
    n_bonds = str(n_bonds)
    static_part = old_string[6:]

    # the first 6 places are reserved for n_atoms and n_bonds, these need to be formattted correctly
    if len(n_atoms) == 1:
        new_string = "  " + n_atoms
    elif len(n_atoms) == 2:
        new_string = " " + n_atoms
    else:
        new_string = n_atoms

    if len(n_bonds) == 1:
        new_string += "  " + n_bonds
    elif len(n_bonds) == 2:
        new_string += " " + n_bonds
    else:
        new_string += n_bonds
    # the rest of the string stays the same and can be added back
    new_string += static_part

    return new_string


def print_correct_connectivity_line(line):
    """A line delimited by 2 spaces is read and correctly formatted to correspond with the official .mol format.
    Each number has 3 allocated spaces, as explained in the docstring of read_connectivity_from_mol_file there are 7
    numbers. This means that the line is formatted as:
    '...''...''...''...''...''...''...' the most right space is allocated to the first number and if the tens are
    reached the middle space is allocated as well, the left space is allocated ones the hundreds are reached.
    So the possibilities are (where \w == a whitespace):
    '\w\w1' '\w10' '100'
    Notice that this function currently only supports indexes that are < 1000.
    :param line:
    :return: correctly formatted line in connectivity block of .mol file
    """
    line_list = line.split("  ")
    to_return = [
        "  " for i in range(16)
    ]  # initialize empty list with 2 spaces as separator

    # 6 cases, either idx 1 is 1 10 100 or idx is 1 10 100
    atom_index_1 = line_list[0]
    if len(atom_index_1) == 3:
        to_return[0] = atom_index_1[0]
        to_return[1] = atom_index_1[1]
        to_return[2] = atom_index_1[2]
    elif len(atom_index_1) == 2:
        to_return[0] = " "
        to_return[1] = atom_index_1[0]
        to_return[2] = atom_index_1[1]
    elif len(atom_index_1) == 1:
        to_return[0] = " "
        to_return[1] = " "
        to_return[2] = atom_index_1

    atom_index_2 = line_list[1]
    if len(atom_index_2) == 3:
        to_return[3] = atom_index_2[0]
        to_return[4] = atom_index_2[1]
        to_return[5] = atom_index_2[2]
    elif len(atom_index_2) == 2:
        to_return[3] = " "
        to_return[4] = atom_index_2[0]
        to_return[5] = atom_index_2[1]
    elif len(atom_index_2) == 1:
        to_return[3] = " "
        to_return[4] = " "
        to_return[5] = atom_index_2

    # the rest of the line is static with 2 spaces as separator
    to_return[7] = line_list[2]
    to_return[9] = line_list[3]
    to_return[11] = line_list[4]
    to_return[13] = line_list[5]
    to_return[15] = line_list[6]

    to_return = "".join(
        [str(elem) for elem in to_return]
    )  # turn list into string
    return to_return


def ff_optimize(source_file, ff_method="uff", list_of_indices_to_freeze=None):
    """Uses openbabenl's ff optimization to locally optimize a molecule and write it back to the same file

    :param source_file:
    :param ff_method:
    :param list_of_indices_to_freeze:
    :return: optimized .mol file with the same filename
    """

    # make molecule from source file with openbabel class
    mol, obconversion = get_mol(source_file)

    # old fashioned method: setup forcefield and do optimization. More customizable
    constr = openbabel.OBFFConstraints()
    # freeze skeleton atoms
    if list_of_indices_to_freeze is not None:
        for idx in list_of_indices_to_freeze:
            # indexing of atoms starts from 1 in openbabel
            constr.AddAtomConstraint(idx + 1)

    forcefield = openbabel.OBForceField.FindForceField(ff_method)
    s = forcefield.Setup(mol, constr)
    if s is not True:
        print("forcefield setup failed.")
        exit()
    else:
        forcefield.SteepestDescent(1000)
        forcefield.GetCoordinates(mol)
    obconversion.WriteFile(mol, source_file)

    # modern method: use pybel's Molecule class and let it do the ff opt, it is already implemented.
    # mol = pybel.Molecule(mol)
    # mol.localopt(ff_method, 1000)
    # mol.write('mol', filename=source_file, overwrite=True)


def xyz_2_smiles(file_name: str) -> str:
    # https://www.kaggle.com/roccomeli/easy-xyz-to-smiles-conversion
    mol = next(pybel.readfile("xyz", file_name))
    smi = mol.write(format="smi")

    return smi.split()[0].strip()


def check_overlap(dataframe_xyz):
    """Given a pandas dataframe, iterates over all the atoms and checks if they're overlapping

    :param dataframe_xyz:
    :return:
    """
    overlap = []
    n_atoms = len(dataframe_xyz)
    overlap_bool = False

    for i in range(0, n_atoms):
        current_atom = dataframe_xyz.iloc[i, :]
        for j in range(i + 1, n_atoms):
            if distance(current_atom, dataframe_xyz.iloc[j, :]) <= 0.8:
                overlap.append([i + 1, j + 1])
                overlap_bool = True
    return overlap_bool, overlap


def remove_hydrogens_and_write(source_mol_file):
    """Use Openbabel's meethods to delete all hydrogens of a molecule and return the same file type with the same
    filename

    :param source_mol_file:
    :return: output file of molecule without hydrogens, same file type as input
    """
    # make molecule from source file with openbabel class
    mol, obconversion = get_mol(source_mol_file)

    # remove all hydrogens and write file
    mol.DeleteHydrogens()
    obconversion.WriteFile(mol, source_mol_file)


def copy_functionalization_list_xyz_2_mol(source_xyz_file, target_mol_file):
    """Copy the functionalization list from a xyz file to a mol file

    :param source_xyz_file:
    :param target_mol_file:
    :return:
    """
    # get the functionalization list from the xyz file
    functionalization_list = open(source_xyz_file).readlines()[1]
    # write the functionalization list to the mol file
    molfile_handle = open(target_mol_file, "r").readlines()
    molfile_handle[0] = functionalization_list
    open(target_mol_file, "w").writelines(molfile_handle)


if __name__ == "__main__":
    # molec = 'H2O'
    # create_molecule_and_write_xyz('H2O', 'substituents_xyz/automatically_generated/' + molec + '.xyz')
    # visualize_xyz_file('substituents_xyz/automatically_generated/something.xyz', True, False)
    # visualize_xyz_file('skeletons/RuPNP_iPr_skl.xyz', False, False)
    # print(convert_list_of_string_to_np_array(['[-0.33332174004836124 0.9428131403470853 0.0]']))
    # print(read_central_atom_index('substituents_xyz/automatically_generated/CH4.xyz'))
    # print(find_distance('substituents_xyz/automatically_generated/CH4.xyz', 2, 3)==1.7473026804689453)
    # print(read_connectivity_from_mol_file('substituents_xyz/manually_generated/F.mol', 1))
    # convert_mol_2_xyz_file('random.mol')
    # convert_xyz_2_mol_file('substituents_xyz/automatically_generated/something_2.xyz')
    # print(print_mol_counts_block(15, 15, 0))
    # print_correct_connectivity_line('120  113  1  0  0  0  0')
    # ff_optimize('Ru(trop2dad)_sigma_func_14.mol', 'uff', None)
    # print(xyz_2_smiles('skeletons/RuPNP_aromatic_tBu.xyz'))
    print(get_bonded_atoms("substituents_xyz/manually_generated/C6H6.mol", 0))
    # xyz_coords_from_mol = pd.read_table('skeletons_temp/PCP-cy.mol', skiprows=4, delim_whitespace=True,
    #                                     names=['x', 'y', 'z', 'atom', 0,1,2,3,4,5,6,7,8,9,10,11,12])
    # print(xyz_coords_from_mol.loc[:, ['x', 'y', 'z']])
