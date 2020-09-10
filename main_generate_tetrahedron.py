# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
import sys
from generate_tetrahedron_folder.generate_tetrahedron import Complex  # conflicting folder name gave error


def main():
    source_file = sys.argv[1]
    target_molecule = sys.argv[2]
    bonded_atom_sub_0_length = sys.argv[3]
    substituent_0 = sys.argv[4]
    sub_0_sub_1_length = sys.argv[5]
    substituent_1_atom = sys.argv[6]
    sub_0_sub_2_length = sys.argv[7]
    substituent_2_atom = sys.argv[8]
    sub_0_sub_3_length = sys.argv[9]
    substituent_3_atom = sys.argv[10]

    substituent = Complex(source_file)
    substituent.generate_substituent_and_write_xyz(target_molecule, bonded_atom_sub_0_length, substituent_0,
                                                   sub_0_sub_1_length, substituent_1_atom,
                                                   sub_0_sub_2_length, substituent_2_atom,
                                                   sub_0_sub_3_length, substituent_3_atom)
    # if len(sys.argv) == 10:
    #     source_file = sys.argv[1]
    #     target_molecule = sys.argv[2]
    #     functionalization_list = None if sys.argv[3].lower() == 'none' else sys.argv[3]
    #     recursive_or_intial = sys.argv[4]
    #     substituent_0 = sys.argv[5]  # will replace atom_to_be_functionalized with this atom
    #     substituent_1_atom = sys.argv[6]
    #     substituent_2_atom = sys.argv[7]
    #     substituent_3_atom = sys.argv[8]
    #     view_xyz_file = True if sys.argv[9].lower() == 'true' else False
    #     save_xyz_picture = True if sys.argv[10].lower() == 'true' else False
    #
    #     substituent = Complex(source_file, functionalization_list, recursive_or_intial)
    #     substituent.generate_substituent_and_write_xyz(target_molecule, substituent_0, substituent_1_atom, substituent_2_atom,
    #                                                    substituent_3_atom, view_xyz_file, save_xyz_picture)
    # else:
    #     source_file = sys.argv[1]
    #     target_molecule = sys.argv[2]
    #     substituent_0 = sys.argv[3]
    #     substituent_1_atom = sys.argv[4]
    #     substituent_2_atom = sys.argv[5]
    #     substituent_3_atom = sys.argv[6]
    #
    #     substituent = Complex(source_file)
    #     substituent.generate_substituent_and_write_xyz(target_molecule, substituent_0, substituent_1_atom, substituent_2_atom,
    #                                                    substituent_3_atom)


if __name__ == "__main__":
    main()
