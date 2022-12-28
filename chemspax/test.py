# -*- coding: utf-8 -*-
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
import pytest
import os
from utilities import *

"""Some testing of the used utilities
"""

test_folder = "tests/"


def create_test_ch3_file():
    if not os.path.exists(test_folder + "CH3.xyz"):
        lines = [
            "4\n",
            "\n",
            " C                 -2.06633184    1.26687443    0.00000000\n",
            " H                 -1.70965897    1.77127288    0.87365134\n",
            " H                 -1.70965897    1.77127288   -0.87365134\n",
            " H                 -3.13633184    1.26688758    0.00000000",
        ]

        wr = open(test_folder + "CH3.xyz", "w")
        wr.writelines(lines)


def create_random_file(filename):
    if not os.path.exists(test_folder + filename):
        lines = ["1", "\n", "2", "\n", "3", "\n"]
        wr = open(test_folder + filename, "w")
        wr.writelines(lines)


def return_last_line(filename):
    f = open(filename, "r")
    return f.readlines()[-1]


# test functions from utilities.py
def test_find_distance():
    create_test_ch3_file()
    assert find_distance(test_folder + "CH3.xyz", 0, 1) == 1.0699999983365585
    os.remove(test_folder + "CH3.xyz")


def test_remove_last_line():
    filename = "random_file.txt"
    create_random_file(filename)
    remove_last_line(test_folder + filename)
    new_last_line = return_last_line(test_folder + filename)
    assert new_last_line == "3"
    os.remove(test_folder + filename)


def test_create_molecule_and_write_xyz():
    filename = "CH4.xyz"
    create_molecule_and_write_xyz("CH4", test_folder + filename)
    lines = open(test_folder + filename).readlines()
    ch4_molecule_list = [
        "5\n",
        'Properties=species:S:1:pos:R:3 pbc="F F F"\n',
        "C        0.00000000       0.00000000       0.00000000\n",
        "H        0.62911800       0.62911800       0.62911800\n",
        "H       -0.62911800      -0.62911800       0.62911800\n",
        "H        0.62911800      -0.62911800      -0.62911800\n",
        "H       -0.62911800       0.62911800      -0.62911800",
    ]
    assert lines == ch4_molecule_list
    os.remove(test_folder + filename)


def test_scale_vector():
    starting_point = np.array([0, 0, 0])
    vector = np.around(
        np.array([-0.62911800, 0.62911800, -0.62911800]), decimals=8
    )
    length = 2.0
    scaled_vector = np.around(
        scale_vector(starting_point, vector, length), decimals=8
    )
    # print(scaled_vector)
    comparison = scaled_vector == np.array(
        [-1.15470054, 1.15470054, -1.15470054]
    )
    assert comparison.all()


def test_convert_list_of_strings_to_np_array():
    converted_array = convert_list_of_string_to_np_array(
        ["[-0.33332174004836124 0.9428131403470853 0.0]"]
    )
    comparison = converted_array == np.array(
        [-0.33332174004836124, 0.9428131403470853, 0.0]
    )
    assert comparison.all()


def test_convert_xyz_2_mol():
    filename = "CH3"
    create_test_ch3_file()
    extension = ".xyz"
    convert_xyz_2_mol_file(test_folder + filename + extension)
    extension = ".mol"
    lines = open(test_folder + filename + extension).readlines()
    lines[1] = "\n"
    ch3_molfile_list = [
        "tests/CH3.xyz\n",
        "\n",
        "\n",
        "  4  3  0  0  0  0  0  0  0  0999 V2000\n",
        "   -2.0663    1.2669    0.0000 C   0  0  0  0  0  3  0  0  0  0  0  0\n",
        "   -1.7097    1.7713    0.8737 H   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "   -1.7097    1.7713   -0.8737 H   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "   -3.1363    1.2669    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "  1  4  1  0  0  0  0\n",
        "  1  2  1  0  0  0  0\n",
        "  3  1  1  0  0  0  0\n",
        "M  END\n",
    ]
    assert lines == ch3_molfile_list
    os.remove(test_folder + filename + extension)
    extension = ".xyz"
    os.remove(test_folder + filename + extension)


def test_print_mol_counts_block():
    old_mol_counts_block = " 59 67  0  0  1  0  0  0  0  0999 V2000"
    mol_counts_block_single_digits = print_mol_counts_block(
        old_mol_counts_block, 1, 1
    )
    assert (
        mol_counts_block_single_digits
        == "  1  1  0  0  1  0  0  0  0  0999 V2000"
    )
    mol_counts_block_double_digits = print_mol_counts_block(
        old_mol_counts_block, 10, 10
    )
    assert (
        mol_counts_block_double_digits
        == " 10 10  0  0  1  0  0  0  0  0999 V2000"
    )
    mol_counts_block_triple_digits = print_mol_counts_block(
        old_mol_counts_block, 100, 100
    )
    assert (
        mol_counts_block_triple_digits
        == "100100  0  0  1  0  0  0  0  0999 V2000"
    )
    with pytest.raises(ValueError):
        print_mol_counts_block(old_mol_counts_block, 1000, 1000)


def test_print_correct_connectivity_line():
    old_connectivity_line_single_digits = "1  1  1  0  0  0  0"
    new_connectivity_line_single_digits = print_correct_connectivity_line(
        old_connectivity_line_single_digits
    )
    assert new_connectivity_line_single_digits == "  1  1  1  0  0  0  0"
    old_connectivity_line_double_digits = "10  1  1  0  0  0  0"
    new_connectivity_line_double_digits = print_correct_connectivity_line(
        old_connectivity_line_double_digits
    )
    assert new_connectivity_line_double_digits == " 10  1  1  0  0  0  0"
