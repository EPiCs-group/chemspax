# auto_func

auto_func should be used to correctly functionalize TM complexes.
Two approaches were considered (currently for tetrahedral substituents only): 
  1) From a library of existing substituents (with a lone pair of electrons) construct the 
  centroid vector and align it with the site to be functionalized on the skeleton and then attach it.
  See attach_substituent.py.
  2) The user inputs an atom to be functionalized on a skeleton complex, 
  then a tetrahedron is formed around this atom and this atom can be replaced. 
  See generate_tetrahedron.py.

## Use case
This code should be able to attach substituents given a functionalization site. In the end this would allow 
for automatic functionalization of complexes.

## Contents
  **skeletons/**
  - contains .xyz files with skeletons of complexes to be functionalized. 
  Currently a nested list (without spaces, otherwise sys.argv() will split the list into 2 arguments)
  with functionalizations is needed on the comment line of the .xyz file.
  In this nested list the first element should be the atom_to_be_functionalized and 
  the second element should be the bonded_atom. 
  
  **substituents_xyz/**  
  - contains .xyz files for substituents that will be tested  
      * automatically_generated/: output of functionalized skeletons by generate_tetrahedron.py and attach_substituent.py.
      * manually_generated/: manually generated substituents that will be used for approach 1. These substituents should
      be an .xyz file with a free bonding site on the central atom of the substituent. 
      So if you want to attach a methyl group; make a CH4 .xyz file, remove one H and put the .xyz file in this folder.
      * old/: contains substituents that became either obsolete or will be dealt with later.
      * visualizations/: contains .png files of the functionalizations made in generate_tetrahedron.py
  
  **attach_substituent.py, run.sh and main_attach_substituent.py** 
  - Class **Substituent** 
  - Load .xyz files and loads substituent in class object 
  - Create first coordination shell to find centroid vector of the whole group
  - Writes group name, central atom and centroid vector to .csv file
  - Class **Complex**
  - Load .xyz file to load substituent and skeleton in class 
  - Generate a matrix of correctly rotated and translated substituent vectors
  - Using the bond length between skeleton_bonded_atom and substituent_central_atom, 
  generate these vectors and write functionalized skeleton to .xyz
  - Takes system arguments in main.py and usage is shown in **run.sh**
   
  
  **generate_tetrahedron.py, run.sh and main_generate_tetrahedron.py** 
  - Load .xyz files and loads complex in class object
  - Find centroid
  - Create and correctly rotate equilateral triangle
  - Write final tetrahedron to .xyz file with 'initial' flag or 
  creates tetrahedron around atom_to_be_functionalized with 'recursive flag' 
  and appends new atoms to source .xyz file
  - Takes system argument in main.py and usage is shown in **run.sh** 
  
  **utilities.py**
  - Uses the ase module for certain utilities
  - Find distances between atoms
  - Visualize molecules
  - Build molecules and write .xyz file from the ase g2 database
  - Read comment line of .xyz files (obsolete, central atom indices used to be written to 
  comment line of .xyz file but are now contained in a .csv file)
  - Remove last line of files 
  
  **exceptions.py**
  - Contains custom made exceptions that can be used for bugfixing
  
  **gjf_to_xyz.py and gjf_to_xyz.sh**
  - Run gjf_to_xyz.py to convert all .gjf files in current path to .xyz files
  
## ToDo  
  **generate_tetrahedron.py**
  - Use ASE to visualize created substituents (done, but make pictures better?)
  - convert .xyz to .gjf (done)  
  - make generate_tetrahedron usable with system arguments (done)
  - add bonded_atom and atom_to_be_functionalized to comment line of .xyz file (done)
  - Test test test
  
  **attach_substituent.py**
  - find centroid vector of substituent group (done)
  - write central atom and centroid vector to .csv (done)
  - attach substituent to skeleton (done)
  - make bash script that is able to use python file and deliver intermediates 
  and final version of functionalized & optimized complex (done)
  - Solve problem with steric hindrance when placing substituents (done)
  - Test test test

## Installation
Python 3.6.0 or higher is required. 

```
git clone https://github.com/EPiCs-group/auto_func

cd auto_func

pip install -r requirements.txt  #to install all required packages.
``` 

It is recommended to use a virtual environment. 
  
  A virtual environment can be created using 
  ```
pip install virtualenv
virtualenv auto_func  # create environment named auto_func
source auto_func/bin/activate  # to activate environment
pip install -r requirements.txt
  ```
  If a built-in virtual environment is used it is necessary to compile openbabel with python bindings from source:  
  The procedure is explained in OpenBabel's 
   [installation guide](https://openbabel.org/wiki/Install_(source_code)#Installing_locally_without_root_access).  
   Perform this compilation in a separate folder to prevent errors.  
   **Note:  In order to install the python bindings the user needs to run:**  
cmake ../ob-src -DCMAKE_INSTALL_PREFIX=/home/noel/tree -DPYTHON_BINDINGS=ON 2>&1 | tee cmake.out
 
   **in step B3.**  
  
   Alternatively, Anaconda can be used to manage virtual environments, **this doesn't require compilation from source of openbabel**:  
     Download the latest installer from 
     [Anaconda's website](https://www.anaconda.com/products/individual).    
     For example: 
   ``` 
   wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh
   ```
   
   Then run the installation script with: 
   ```
     bash Anaconda3-2020.07-Linux-x86_64.sh  
   ```
   The installation is pretty self explanatory, afterwards create a virtual environment and activate it.  
The environment will be named 'auto_func', the conda_env.yml file can be changed if a different name is required. 
        
  ```
  conda env create -f conda_env.yml
  
  conda activate auto_func
  ```
  Check if the environment is installed correctly by running 
  ```
  conda env list
  ```

~~REDUNDANT~~   
~~After OpenBabel is compiled succesfully the packages can be downloaded from the conda-forge channel
using the requirements.txt file, in the auto_func folder run:
conda install -c conda-forge --file requirements.txt~~
   
  More information on virtual environments can be found at the
  [venv homepage](https://docs.python.org/3/library/venv.html) or 
  [this Anaconda cheatsheet](https://uoa-eresearch.github.io/eresearch-cookbook/recipe/2014/11/20/conda/)
## Instructions
An example .xyz file of a skeleton in skeletons/ looks like this: 
![example](images/example_functionalization_list.jpg)
  
Make sure to put no spaces in this functionalization_list since 
sys.argv() is used in the main.py files, each space is thus interpreted as a new argument
to the function. Also note that index 0 is the **first atom in the list**.
In the example figure above, the index of Ru == 0.  
  
In the functionalization_list each list is a functionalization. 
The first item of each list is the atom_to_be_functionalized (the atom that will be removed).
![example](images/example_atom_to_be_functionalized.jpg)
  
  The second item of each list is the bonded_atom, this is the 
  atom that is bonded to the atom_to_be_functionalized and thus the 
  skeleton atom that will be bonded to the substituent. 
  
  ![example](images/example_bonded_atom.jpg)
  
  
  With these instructions the user should be able to 
  write the functionalization_list of the skeletons correctly 
  in skeletons/ (**note: no newline at the end of the file!**).
  For the 2 different approaches of functionalization an 
  explanation should be given on how to proceed.
  
  **generate_tetrahedron.py**
  
  **Note: Make sure to check the relative path on line 124 
  is correctly set, relatively to your working directory.**
  
To functionalize a skeleton using this approach, now only the bond lengths
between bonded_atom and central_atom of substituent and 
surrounding atoms of the substituent with central_atom of substituent is needed.
The main_generate_tetrahedron.py takes input values as follows:
  
  1) path_to_source_file
  2) target_name (**Note: name only, the file will always be placed in 
  substituents_xyz/automatically_generated with extension .xyz**)      
  3) bond length of new central atom of substituent with bonded_atom
  4) actual new central atom of substituent 
  5) bond length of new central atom of substituent with first surrounding
  substituent atom
  6) actual first surrounding substituent atom
  7) same as 5) for second surrounding substituent atom
  8) same as 6) for second surrounding substituent atom
  9) same as 5) for third surrounding substituent atom
  10) same as 6) for third surrounding substituent atom
  
  Where each input value is given as a space delimited system argument.
  The user can modify the example of generate_tetrahedron_folder/run.sh
  to their needs or use generate_tetrahedron.py directly.
  
  **attach_substituent.py**
  
  **Note: Make sure to check the relative path on line 175, 88, 66 
  and 21 are correctly set, relatively to your working directory.**
  
  The user can upload manually generated substituents 
  to substituents_xyz/manually_generated/ or use the pre-made substituents
  contained in that folder. Then the data_preparation.py script can be run 
  to generate a .csv database
  of the central atom and centroid vector per substituent as these will
  be used to align the substituent with the skeleton. This also generates .mol files
  for the skeleton and substituents since these are used for their connectivity data.
  **data_preparation.py assumes that the central atom of the substituent is the first atom in the 
  .xyz file of the substituent! Change this if it's necessary.**
  An example of 
  the .csv database is shown below:
  ![example](images/example_csv_database.jpg)
  
  
  Like in generate_tetrahedron.py the Complex.generate_substituent_and_write_xyz()
  function is used to functionalize a skeleton. The main_attach_substituent.py takes the following
  input:
  1) path_to_source_file
  2) target_name (**Note: name only, the file will always be placed in 
  substituents_xyz/automatically_generated with extension .xyz**)
  3) name of substituent group (same as key in .csv database)
  4) relative path to .csv database file
  5) bond length between central atom of the substituent and bonded_atom of skeleton 
  6) whether the user wants to use the python script with the xtb bash script, in the 
  xtb bash script the conversion of the optimized .mol file to .xyz file doesn't happen in python.
  If this is set to false only ff optimization will be done and the python script will handle file conversions.
  
  Where each input value is given as a space delimited system argument.
  The user can modify the example of attach_substituent_folder/run.sh
  to their needs or use attach_substituent.py directly.
  
## Example usage
  - ```bash generate_tetrahedron_folder/run.sh``` to functionalize a RuPNP skeleton where 4 H sites 
  will be substituted for methyl groups using generate_tetrahedron.py
  - ```bash attach_substituent_folder/run.sh``` to functionalize a RuPNP skeleton where
  6 sites will be functionalized with substituents
## Authors
Adarsh Kalikadien & Vivek Sinha