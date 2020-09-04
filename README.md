# correct_rotation_substituent

correct_rotation_substituent should be used to correctly functionalize TM complexes.
Two approaches were considered (currently for tetrahedral substituents only): 
  1) From a library of existing substituents (with 1 lone pair of electrons) construct the 
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
  - contains .xyz files with skeletons of complexes to be functionalized 
  
  **substituents_xyz/**  
  - contains .xyz files for substituents that will be tested  
      * automatically_generated/: output of functionalized skeletons by generate_tetrahedron.py
      * manually_generated/: manually generated substituents that will be used for approach 1.
      * old/: contains substituents that became either obsolete or will be dealt with later.
      * visualizations/: contains .png files of the functionalizations made in generate_tetrahedron.py
  
  **attach_substituent**  
  - Load .xyz files and loads substituent in class object 
  - Create first coordination shell to find centroid vector of the whole group
  - Writes group name, central atom and centroid vector to .csv file
   
  
  **generate_tetrahedron.py, run.sh and main.py** 
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
  
 ### Explanation of generating a tetrahedron 
  Given a vector (v), say a H-C bond length where C is to be functionalized and converted into a tetrahedral group.
  This code generates the tetrahedral functionalization.

  Explanation: If one edge of a tetrahedral is determined already, one only needs to place an equilateral triangle such that
  the centroid of the triangle is at a distance of b/3 where b is the bond length from C to be functionalized.
  Another challenge is to oreate the triangle such that the normal vector n of the triangle has a 0 (or 180 degrees) angle
  with the edge of the tetrahedral which is given. That is n.v = |n||v|.
  This achieved via a rotation matrix which rotates the equilateral triangle.
  The equilateral triangle is pre-defined with origin as center and all edges are in the XY plane (z = 0) the normal vector
  n is set to be 0,0,1.
  **The final distances seem to have an error 0.1 probably from the errors in arccos and arctan calculations**
  **or there is some silly mistake somewhere if the edge of equilateral triangle is a then a = 2sqrt(2/3)b **

## ToDo  
  **generate_tetrahedron.py**
  - Use ASE to visualize created substituents (done, but make pictures better?)
  - convert .xyz to .gjf (done)  
  - make generate_tetrahedron usable with system arguments (done)
  - add bonded_atom and atom_to_be_functionalized to comment line of .xyz file (done)
  - make bash script that is able to use python file and deliver intermediates 
  and final version of functionalized & optimized complex 
  using generate_tetrahedron.py (partly done, optimization can be added) 
  - Test test test
  
  **attach_substituent.py**
  - find centroid vector of substituent group (done)
  - write central atom and centroid vector to .csv (done)
  - attach substituent to skeleton 
  - Test test test
## Example usage
  - ```./run.sh``` to functionalize a RuPNP skeleton where 4 H sites 
  will be substituted for methyl groups using generate_tetrahedron.py
## Authors
Adarsh Kalikadien & Vivek Sinha