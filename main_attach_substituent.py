# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
from attach_substituent_folder.attach_substituent import Complex
import sys


def main():
    source_data = sys.argv[1]
    target_file = sys.argv[2]
    substituent_molecule = sys.argv[3]
    path_to_database = sys.argv[4]
    length_skeleton_bonded_atom_substituent_central_atom = sys.argv[5]
    use_xtb_script_after = True if sys.argv[6].lower() == 'true' else False

    some_complex = Complex(source_data, substituent_molecule, path_to_database)
    some_complex.generate_substituent_and_write_xyz(target_file, length_skeleton_bonded_atom_substituent_central_atom,
                                                    use_xtb_script_after)


if __name__ == '__main__':
    main()
