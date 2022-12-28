# -*- coding: utf-8 -*-
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
"""This script was used on the Cartesius cluster to remove hydrogens from a molecule to calculate the
HRMSD, this functionality was also incorporated as a function in utilities.py (note: this script writes the OBMol object
to a different filename!)
"""
from openbabel import openbabel
import sys

from chemspax.utilities import get_mol

# reduce warnings from openbabel
ob_log_handler = openbabel.OBMessageHandler()
ob_log_handler.SetOutputLevel(0)
openbabel.obErrorLog.SetOutputLevel(0)


def remove_hydrogens_and_write(source_mol_file):
    """Use Openbabel's meethods to delete all hydrogens of a molecule and return the same file type with a different
    filename

    :param source_mol_file:
    :return: output file of molecule without hydrogens, same file type as input
    """
    # make molecule from source file with openbabel class
    mol, obconversion = get_mol(source_mol_file)

    # remove all hydrogens and write file
    mol.DeleteHydrogens()
    obconversion.WriteFile(mol, source_mol_file[:-4] + "_no_H" + extension)


if __name__ == "__main__":
    source_file = sys.argv[1]
    remove_hydrogens_and_write(source_file)
