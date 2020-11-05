# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
"""This script was used on the Cartesius cluster to extract bond distances for given atoms, this functionality was
also incorporated as a function in utilities.py
"""
from openbabel import openbabel
import sys

# reduce warnings from openbabel
ob_log_handler = openbabel.OBMessageHandler()
ob_log_handler.SetOutputLevel(0)
openbabel.obErrorLog.SetOutputLevel(0)


def get_neighbour_bond_distance_search(source_mol_file, central_atom_atomic_num, search_this_atomic_num):
    """Use Openbabel's meethods to find the distance of a specific atom type to a given atom type, the atoms should be
    bonded and the specific atom type that is searched for should be searched by using it's atomic number.
    This function is useful if the given atom does not have a known index in the molecule,
    if you know the index of the central atom use get_neighbour_distance.py

    :param source_mol_file:
    :param central_atom_atomic_num:
    :param search_this_atomic_num:
    :return:
    """
    central_atom_atomic_num = int(central_atom_atomic_num)
    search_this_atomic_num = int(search_this_atomic_num)

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

    # iterate over atoms in molecule to search the central atom
    for atom in openbabel.OBMolAtomIter(mol):
        atomic_number = atom.GetAtomicNum()
        if atomic_number == central_atom_atomic_num:
            for neighbour_atom in openbabel.OBAtomAtomIter(atom):
                atomic_number = neighbour_atom.GetAtomicNum()
                if atomic_number == search_this_atomic_num:
                    distance = neighbour_atom.GetDistance(atom)
                    print(distance)


if __name__ == '__main__':
    source_file = sys.argv[1]
    central_atom_atomic_number = sys.argv[2]
    atomic_number_to_look_for = sys.argv[3]
    get_neighbour_bond_distance_search(source_file, central_atom_atomic_number, atomic_number_to_look_for)