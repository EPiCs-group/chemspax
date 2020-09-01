# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
from find_centroid import Substituent
from generate_tetrahedron import Complex
import sys


def main():
    source_file = sys.argv[1]
    target_molecule = sys.argv[2]
    functionalization_list = sys.argv[3]
    recursive_or_intial = sys.argv[4]
    substituent_1 = sys.argv[5]
    substituent_2 = sys.argv[6]
    substituent_3 = sys.argv[7]
    view_xyz_file = True if sys.argv[8].lower() == 'true' else False
    save_xyz_picture = True if sys.argv[9].lower() == 'true' else False
    
    substituent = Complex(source_file, functionalization_list, recursive_or_intial)
    substituent.generate_substituent_and_write_xyz(target_molecule, substituent_1, substituent_2, substituent_3,
                                                   view_xyz_file, save_xyz_picture)


if __name__ == "__main__":
    main()