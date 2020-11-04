# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
"""This script was used on the Cartesius cluster to extract bond distances for given atoms, this functionality was
also incorporated as a function in utitilities.py
"""
from openbabel import openbabel
import sys

# reduce warnings from openbabel
ob_log_handler = openbabel.OBMessageHandler()
ob_log_handler.SetOutputLevel(0)
openbabel.obErrorLog.SetOutputLevel(0)


def get_neighbour_bond_distance(source_mol_file, atom_index, search_this_atomic_num):
    """Use Openbabel's meethods to find the distance of a specific atom to a given atom, the atoms should be bonded
    and the atom should be searched by using it's atomic number

    :param source_mol_file:
    :param atom_index:
    :param search_this_atomic_num:
    :return: list of distance between indicated central atom and atoms that are bonded to it with given atomic number
    """
    atom_index = int(atom_index)
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

    # make atom object for atom of whose bonds we want to look at
    atom = mol.GetAtom(atom_index + 1)  # for obmol get functions indexing starts from 1

    for neighbour_atom in openbabel.OBAtomAtomIter(atom):
        atomic_number = neighbour_atom.GetAtomicNum()
        # change this atomic number to change the atom you're looking for
        if atomic_number == search_this_atomic_num:
            # if this is the atom we're looking for, append the distance
            distance = neighbour_atom.GetDistance(atom)
            print(distance)


if __name__ == '__main__':
    source_file = sys.argv[1]
    central_atom_index = sys.argv[2]
    atomic_number_to_look_for = sys.argv[3]
    get_neighbour_bond_distance(source_file, central_atom_index, atomic_number_to_look_for)