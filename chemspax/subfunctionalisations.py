#!/bin/python3

# Modules
from rdkit import Chem
import numpy as np
import os
import subprocess

# ChemSpaX setup
import glob
from chemspax.main import main

# Own functions
from pick_subfunctionalisations import *

"""Subfunctionalisations per structure"""
sub_sites = {'hydrazone': [H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl],
             'semicarbazone':[H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl],
             'aminal':[H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl],
             '12-diol':[H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl],
             '12-aminoalcohol':[H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl],
             'hydroxylamine':[H_alkyl_aryl,H_alkyl_aryl],
             'iminohetarene':[H_alkyl_aryl],
             'thiourea': [H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl],
             'isothiourea': [H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl],
             'guanidine': [H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl],
             'alpha-aminoacid': [H_alkyl_aryl],
             'lactone': [H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl],
             'thiolactone': [],
             'thioketone': [alkyl_aryl],
             'prim-alcohol':[],
             'carbonicaciddiester':[alkyl_aryl],
             'thiocarbonicacidmonoester':[],
             'imidothioester':[alkyl_aryl,H_alkyl_aryl],
             'thiocarboxylicacidderiv':[heteroatom],
             'alkyliodide':[],
             'aryliodide': [H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl],
             'sulfonicacidderiv':[heteroatom],
             'phosphoricacidester':[heteroatom,heteroatom],
             'thiophosphoricacidamide':[H_alkyl_aryl,heteroatom,heteroatom],
             'phosphonicacidderiv':[heteroatom,heteroatom]}

def initialise(folders: list):
    "This function empties all folders to prevent confusion from data from a previously failed calculation"
    # Function must be here to prevent circular import
    def delete_items_in_folder(folder_path):
        items = os.listdir(folder_path)
        for item in items:
            item_path = os.path.join(folder_path, item)
            if os.path.isfile(item_path):
                os.remove(item_path)
            else:
                delete_items_in_folder(item_path)
    
    for folder in folders:
        folder_path = os.path.join(os.getcwd(), folder)
        if os.path.isfile(folder_path):
            os.remove(folder_path)
        else:
            delete_items_in_folder(folder_path)

"""Information about sub sites dictionary for main program"""
def len_sub_sites():
    return len(sub_sites)

def sub_sites_name(n:int):
    return str(list(sub_sites.keys())[n])

def sub_sites_value_length(key:str):
    return len(sub_sites[key])


# Remove dummy atoms for final functionalisation
def remove_dummy_atom_MOL(working_directory,input_path,subsub,subfunc,dummy_atom="Li"):
    "This function removes the dummy atom from the MDL Molfile"
    if subfunc==True:
        os.chdir(os.path.join(working_directory, "tmp"))   
        subprocess.run("xtb "+input_path+ " --opt > /dev/null 2>&1",shell=True)
        os.chdir(working_directory)
        mol = Chem.MolFromMolFile(os.path.join(working_directory,"tmp","xtbopt.mol"), removeHs=False)
    else:
        mol = Chem.MolFromMolFile(input_path, removeHs=False)
    
    mol = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts('['+dummy_atom+']'))
    Chem.MolToMolFile(mol, os.path.join(working_directory,'functionalisations','sub'+str(subsub)+'.mol'),kekulize=False)

def remove_dummy_atom_XYZ(working_directory,input_path,subsub,dummy_atom="Li"):
    "This function removes the dummy atom from the XYZ file"
    f = open(input_path)
    xyz_file=f.readlines(0)
    f.close()
    f = open(os.path.join(working_directory,'functionalisations','sub'+str(subsub)+'.xyz'),'w')
    f.write(str(int(xyz_file[0])-1)+'\n')
    for line in xyz_file[1:]:
        if dummy_atom not in line:
            f.write(line)
    f.close()



def subfunctionalisation(chosen_sub,chosen_sub_name,working_directory,path_to_output):
    # Remove old files from functionalisations
    initialise(["tmp"])

    # Select sub substitutions
    substituent_list=[]
    for subsub in range(sub_sites_value_length(chosen_sub_name)):
        substituent_list.append(sub_sites[chosen_sub_name][subsub](np.random.randint(1,7)))
    print("Subfunctionalisations: ", substituent_list)
    # Setup ChemSpaX variables
    path_to_substituents = os.path.join(working_directory, "subfunctionalisations")
    path_to_database = os.path.join(path_to_substituents, "central_atom_centroid_database.csv")
    path_to_skeletons = os.path.join(working_directory, "subskeletons")
    skeleton_list = glob.glob(os.path.join(working_directory, "subskeletons", chosen_sub_name+'.xyz'))
    
    # Run ChemSpaX
    main(skeleton_list, substituent_list, path_to_database, path_to_substituents, path_to_skeletons, working_directory, path_to_output,True)

    # Save structures without Li dummy atom
    remove_dummy_atom_MOL(working_directory,os.path.join(working_directory, "tmp",chosen_sub_name+"_func_"+str(len(sub_sites[chosen_sub_name]))+".mol"),chosen_sub,True)
    remove_dummy_atom_XYZ(working_directory,os.path.join(working_directory, "tmp",chosen_sub_name+"_func_"+str(len(sub_sites[chosen_sub_name]))+".xyz"),chosen_sub)