# -*- coding: utf-8 -*- 
#                                                                                   #
#  __authors__ = Mark Heezen                                                        #
#  __institution__ = Vrije Universiteit Brussel & Universidad Autónoma de Madrid    #
#  

import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

import glob
import numpy as np
import shutil

from main import main
from data_preparation import prepare_data
from subfunctionalisations import initialise, sub_sites_value_length, remove_dummy_atom_MOL, remove_dummy_atom_XYZ
from pick_subfunctionalisations import H_alkyl_aryl
from run_structures import check_skeleton_sites

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

def final_functionalisation(skeleton):
    substituent_list=['sub0'] # Initialise functionalisation list
    for position in range(1,int(func_sites[list(func_sites.keys())[skeleton]])): # Add the right number of substituents
        substituent_list.append('sub'+str(position))

    # ChemSpaX setup
    path_to_substituents = os.path.join(working_directory, "functionalisations/")
    path_to_database = os.path.join(path_to_substituents, "central_atom_centroid_database.csv")
    path_to_skeletons = os.path.join(working_directory, "skeletons")
    skeleton_list = glob.glob(os.path.join(working_directory, "skeletons", str(list(func_sites.keys())[skeleton])+'.xyz'))
    prepare_data(path_to_substituents, path_to_skeletons, path_to_database)
    
    # Run ChemSpaX
    main(skeleton_list, substituent_list, path_to_database, path_to_substituents, path_to_skeletons, working_directory, path_to_output,False)
    
    # Copy files to result folder
    shutil.copy(os.path.join(os.getcwd(),'tmp',str(list(func_sites.keys())[skeleton])+'_func_'+str(func_sites[list(func_sites.keys())[skeleton]])+'.mol'),os.path.join(os.getcwd(),'results','struc_'+str(struccounter)+'.mol'))
    shutil.copy(os.path.join(os.getcwd(),'tmp',str(list(func_sites.keys())[skeleton])+'_func_'+str(func_sites[list(func_sites.keys())[skeleton]])+'.xyz'),os.path.join(os.getcwd(),'results','struc_'+str(struccounter)+'.xyz'))

if __name__ == "__main__":
    """The test performs necessary subfunctionalisations on the first substituent in sub_sites below and attaches it to the alphabetically first substituent in the substituents folder."""
    # Default program settings
    working_directory = os.getcwd()
    path_to_output = os.path.join(working_directory, "tmp")
    path_to_substituents = os.path.join(working_directory, "subfunctionalisations/")
    path_to_database = os.path.join(path_to_substituents, "central_atom_centroid_database.csv")
    path_to_skeletons = os.path.join(working_directory, "subskeletons")

    # Example settings
    sub_sites = {'hydrazone': [H_alkyl_aryl,H_alkyl_aryl,H_alkyl_aryl]}
    func=1
    skeleton=0
    struccounter=0

    # Example program
    initialise(["tmp","functionalisations","results"])
    prepare_data(path_to_substituents, path_to_skeletons, path_to_database)
    func_sites=check_skeleton_sites('skeletons')
    
    

    while True: # This makes sure the script tries again if a functionalisation fails
        try:
            print("START TEST")
            for chosen_sub in range(func_sites[list(func_sites.keys())[skeleton]]): # Over the number of functionalisation sites
                random_index = np.random.randint(0, len(sub_sites)) # Pick a substituent
                chosen_sub_name=str(list(sub_sites.keys())[random_index])

                if sub_sites_value_length(chosen_sub_name) > 0: # Check if subfunctionalisation is necessary
                    print("Chosen substituent: ",chosen_sub,chosen_sub_name,'subfunctionalisation: yes')
                    subfunctionalisation(chosen_sub,chosen_sub_name,working_directory,path_to_output) # Perform subfunctionalisation
                else:
                    print("Chosen substituent: ",chosen_sub,chosen_sub_name,'subfunctionalisation: no')
                    remove_dummy_atom_MOL(working_directory,os.path.join(working_directory, "subskeletons",chosen_sub_name+".mol"),chosen_sub,False)
                    remove_dummy_atom_XYZ(working_directory,os.path.join(working_directory, "subskeletons",chosen_sub_name+".xyz"),chosen_sub)
                
            # Final functionalisation of the just created substituents to the initial skeleton
            initialise(["tmp"])
            final_functionalisation(skeleton)
            initialise(["tmp","functionalisations"])
            break
        except:
            print("Functionalisation failed, trying again")
            initialise(["tmp","functionalisations"])
                    
    print("END TEST. The functionalised structure can be found in the results folder.")