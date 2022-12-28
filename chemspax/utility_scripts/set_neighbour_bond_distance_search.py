# -*- coding: utf-8 -*-
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
"""This script was used on the Cartesius cluster to scale bond distances for given atoms, this functionality was
also incorporated as a function in utilities.py
"""
from openbabel import openbabel
import numpy as np
import sys

from chemspax.utilities import get_mol

# reduce warnings from openbabel
ob_log_handler = openbabel.OBMessageHandler()
ob_log_handler.SetOutputLevel(0)
openbabel.obErrorLog.SetOutputLevel(0)


def scale(starting_point, vector, length):
    vector = vector / np.linalg.norm(vector)
    return starting_point + vector * length  # length in angstrom


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
                    new_vector = scale(
                        central_atom_vector, distance_vector, new_length
                    )
                    neighbour_atom.SetAtomicNum(new_atomic_num)
                    neighbour_atom.SetVector(
                        new_vector[0], new_vector[1], new_vector[2]
                    )
                    obconversion.WriteFile(mol, source_mol_file)


if __name__ == "__main__":
    source_file = sys.argv[1]
    central_atom_atomic_number = sys.argv[2]
    atomic_number_to_look_for = sys.argv[3]
    new_bond_length = sys.argv[4]
    new_atomic_number = sys.argv[5]
    set_neighbour_bond_distance_search(
        source_file,
        central_atom_atomic_number,
        atomic_number_to_look_for,
        new_bond_length,
        new_atomic_number,
    )
