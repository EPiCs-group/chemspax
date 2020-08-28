# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
from find_centroid import Substituent
from generate_tetrahedron import Complex


def main():
    # make methyl class to find bond between C-H
    methyl = Complex('substituents_xyz/manually_generated/CH3.xyz')
    # Write to .xyz file and visualize, H is in center and C is pointing outwards and will be attached to catalyst
    methyl.generate_substituent_and_write_xyz('CH4', 'H', 'H', 'H', True, False)
    # find centroid for this newly generated substituent
    substituent = Substituent('substituents_xyz/automatically_generated/CH4.xyz', 1.7)
    centroid, indices = substituent.first_coordination()
    print(centroid)


if __name__ == "__main__":
    main()