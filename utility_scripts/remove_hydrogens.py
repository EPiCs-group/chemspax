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
    # initalize openbabel classes
    obconversion = openbabel.OBConversion()
    # both xyz and mol can be used as input but mol contains an accurate graph representation
    if source_mol_file[-4:] == '.mol':
        extension = '.mol'
        obconversion.SetInFormat('mol')
        obconversion.SetOutFormat('mol')
    elif source_mol_file[-4:] == '.xyz':
        extension = '.xyz'
        obconversion.SetInFormat('xyz')
        obconversion.SetOutFormat('xyz')
    else:
        raise Exception('file type is incorrect, .mol and .xyz are supported, not', source_mol_file[-4:])
    # intialize class and read file
    mol = openbabel.OBMol()
    obconversion.ReadFile(mol, source_mol_file)
    # remove all hydrogens and write file
    mol.DeleteHydrogens()
    obconversion.WriteFile(mol, source_mol_file[:-4] + '_no_H' + extension)


if __name__ == '__main__':
    source_file = sys.argv[1]
    remove_hydrogens_and_write(source_file)