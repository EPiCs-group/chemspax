# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
from attach_substituent_folder.attach_substituent import Complex
import sys
"""Main class that takes command line arguments and passes them to attach_substituent.py's Complex class
"""


def main():
    original_skeleton_name = sys.argv[1]
    source_data = sys.argv[2]
    target_file = sys.argv[3]
    substituent_molecule = sys.argv[4]
    path_to_database = sys.argv[5]
    length_skeleton_bonded_atom_substituent_central_atom = sys.argv[6]
    functionalize_all_h = True if sys.argv[7].lower() == 'true' else False
    use_xtb_script_after = True if sys.argv[8].lower() == 'true' else False

    some_complex = Complex(original_skeleton_name, source_data, substituent_molecule, path_to_database)
    some_complex.generate_substituent_and_write_xyz(target_file, length_skeleton_bonded_atom_substituent_central_atom,
                                                    functionalize_all_h, use_xtb_script_after)


if __name__ == '__main__':
    main()
