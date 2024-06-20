# -*- coding: utf-8 -*- 
#                                                                                   #
#  __authors__ = Mark Heezen                                                        #
#  __institution__ = Vrije Universiteit Brussel & Universidad Autónoma de Madrid    #
#  

"""Program setup"""
# Modules import
import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

import numpy as np
import argparse
import shutil
import datetime
import resource
from subfunctionalisations import *

# ChemSpaX setup
import glob
from main import main
from data_preparation import prepare_data
from rdkit import Chem



"""Initialisation"""
def check_skeleton_sites(folder_name,dummy_atom='Br',test=False):
    """This function counts the number of functionalisation sites (by default the number of Br atoms) per skeleton"""
    func_sites={}
    for file in sorted(os.listdir(os.path.join(os.getcwd(),folder_name))):
        filename = os.fsdecode(file)
        if filename.endswith(".xyz"): 
            try:
                f = open(os.path.join(os.getcwd(),folder_name,filename))
                xyz=f.readlines()
                counter=0
                for line in xyz:
                    if dummy_atom in line:
                        counter+=1
                func_sites[filename[:-4]]=counter
            except:
                FileExistsError("Something wrong with ", filename)
    return func_sites


def parse_command_line():
    """This function parses the command line and returns an array with the number of structures to generate per skeleton"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--number', help='Number of molecules to generate', required=False, action='store',
                            default=10)
    total_strucs = parser.parse_args().number
    strucs=np.zeros(len(func_sites))

    strucs[:]=int(int(total_strucs)/(len(func_sites)))
    left_over=int(total_strucs)-np.sum(strucs)
    for i in range(int(left_over)):
        strucs[i]+=1
    return strucs

"""(Sub)functionalisations"""

def final_functionalisation(skeleton):
    print("Starting on final functionalisation for struc",str(struccounter))
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
    print("END: struc",str(struccounter),"\n \n")


def close_all_files():
    """This function closes all open files python opened to prevent a 'too many files open' error"""
    # Get a list of all file descriptors currently open
    open_fds = set()
    for fd in range(os.sysconf("SC_OPEN_MAX")):
        try:
            open_fds.add(os.path.realpath(os.readlink(f"/proc/self/fd/{fd}")))
        except Exception:
            pass

    # Close all open files
    for fd in open_fds:
        try:
            os.close(fd)
        except Exception:
            pass


if __name__ == "__main__":
    """Initialisation"""
    start=datetime.datetime.now()

    working_directory = os.getcwd()
    path_to_output = os.path.join(working_directory, "tmp")
    path_to_substituents = os.path.join(working_directory, "subfunctionalisations/")
    path_to_database = os.path.join(path_to_substituents, "central_atom_centroid_database.csv")
    path_to_skeletons = os.path.join(working_directory, "subskeletons")
    prepare_data(path_to_substituents, path_to_skeletons, path_to_database)

    print('''      -----------------------------------------------------------      
        |                   =====================                   |     
        |                    Recursive ChemSpaX                     |     
        |                   =====================                   |     
        |                         M. Heezen                         |     
        |Master in Theoretical Chemistry and Computational Modelling|     
        |                 Master Thesis, April 2024                 |     
        ----------------------------------------------------------- \n''')

    print("Program started at ",str(start))

    print("\n START: initialisation")
    initialise(["tmp","functionalisations","results"])
    print("Folders cleared")
    func_sites=check_skeleton_sites('skeletons')
    print("Functionalisation sites found")
    strucs=parse_command_line()
    print("Command line parsed")
    print("END: initialisation \n \n")


    """Functionalisations"""
    print("START: FUNCTIONALISATIONS")
    struccounter=0

    for skeleton in range(len(func_sites)): # Over all skeletons
        for func in range(int(strucs[skeleton])): # Over the number of times you want this skeleton
            while True: # This makes sure the script tries again if a functionalisation fails
                try:
                    print("START: struc",str(struccounter),'skeleton: ',str(list(func_sites.keys())[skeleton]))
                    for chosen_sub in range(func_sites[list(func_sites.keys())[skeleton]]): # Over the number of functionalisation sites
                        random_index = np.random.randint(0, len_sub_sites()) # Pick a substituent
                        chosen_sub_name=sub_sites_name(random_index)

                        if sub_sites_value_length(chosen_sub_name) > 0: # Check if subfunctionalisation is necessary
                            print("Chosen substituent: ",chosen_sub,chosen_sub_name,'subfunctionalisation: yes')
                            subfunctionalisation(chosen_sub,chosen_sub_name,working_directory,path_to_output)
                        else:
                            print("Chosen substituent: ",chosen_sub,chosen_sub_name,'subfunctionalisation: no')
                            remove_dummy_atom_MOL(working_directory,os.path.join(working_directory, "subskeletons",chosen_sub_name+".mol"),chosen_sub,False)
                            remove_dummy_atom_XYZ(working_directory,os.path.join(working_directory, "subskeletons",chosen_sub_name+".xyz"),chosen_sub)
                        
                    # Final functionalisation
                    initialise(["tmp"])
                    final_functionalisation(skeleton)
                    struccounter+=1
                    initialise(["tmp","functionalisations"])
                    break
                except:
                    print("Functionalisation failed, trying again")
                    initialise(["tmp","functionalisations"])
                    
        close_all_files()
    print("END: FUNCTIONALISATIONS")

    """Termination"""
    finish=datetime.datetime.now()
    print("Successful  termination at ",finish, "after",str(finish-start))