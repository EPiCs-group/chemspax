# -*- coding: utf-8 -*-
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
import ast
import glob
import os

"""Script to automatically replace H with X in mol file of a skeleton in current path for a given functionalization 
strategy or to  just replace all hydrogens of the skeleton with X. This mol file can then be used for 
visualization purposes.

The script assumes that the functionalization_list is on the first line of the mol file
"""


def find_n_atoms(filename):
    atoms_and_bonds_block = open(filename).readlines()[3]
    n_atoms = int(atoms_and_bonds_block[:3].replace(" ", ""))
    return n_atoms


def replace_func_strategy_with_x(source_filename):
    # functionalization list on first line of .mol file
    functionalization_list = ast.literal_eval(open(source_filename).readline())
    atoms_to_be_functionalized = []

    for item in functionalization_list:
        # add atom_to_be_functionalized to a list (first item of nested list)
        atoms_to_be_functionalized.append(item[0])

    lines = open(source_filename).readlines()
    comment_n_atom_lines = lines[:4]

    atoms_and_connectivity = lines[4:]
    for idx in atoms_to_be_functionalized:
        # if atom has a single letter (which is the case for H) it is located on the 31st line
        # str does not support item assignment
        lines_list = list(atoms_and_connectivity[idx])
        lines_list[31] = "X"
        lines_str = "".join(lines_list)
        atoms_and_connectivity[idx] = lines_str

    total_lines = comment_n_atom_lines + atoms_and_connectivity
    target_filename = source_filename[:-4] + "_func_strat_x.mol"
    with open(target_filename, "w") as wr:
        wr.writelines(total_lines)


def replace_all_h_with_x(source_filename):
    n_atoms = find_n_atoms(source_filename)
    n_comment_lines = 4

    # read everything, correcting for shift in indices from comment lines
    lines = open(source_filename).readlines()
    comment_n_atom_lines = lines[:n_comment_lines]
    atoms_block = lines[n_comment_lines : n_comment_lines + n_atoms]
    connectivity_block = lines[n_comment_lines + n_atoms :]

    for idx, line in enumerate(atoms_block):
        line_list = list(line)
        # print('line_list[31] =', line_list[31])
        if line_list[31] == "H":
            line_list[31] = "X"
            atoms_block[idx] = "".join(line_list)

    total_lines = comment_n_atom_lines + atoms_block + connectivity_block
    target_filename = source_filename[:-4] + "_all_h_x.mol"
    with open(target_filename, "w") as wr:
        wr.writelines(total_lines)


if __name__ == "__main__":
    for file in glob.glob("*x.mol"):
        os.remove(file)

    for file in glob.glob("*mol"):
        replace_func_strategy_with_x(file)
        replace_all_h_with_x(file)
