# correct_rotation_substituent

correct_rotation_substituent should be used to correctly functionalize TM complexes.

## Use case
The code should be able to generate a library of substituents and be 
able to find the centroid of a substituent. In the end this would allow 
for automatic functionalization of complexes.

## Contents
  **substituents_xyz/**  
  - contains .xyz files for substituents that will be tested  

  **correct_rotation.py**  
  - load .xyz files and loads substituent in class object 
  - Create first coordination shell
  - check distance between molecules  
  
  **generate_tetrahedron.py** 
  - load .xyz files and loads complex in class object
  - find centroid
  - create and correctly rotate equilateral triangle
  - write final tetrahedron to .xyz file  
   
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
  - Make *generate_substituent_vectors()* look nicer (done)
  - Add check if c in *generate_substituent_vectors()* doesn't become 0 (done)
  - In *generate_tetrahedron.py* add central atom to *write_xyz()*  (done manually, think about automating this)
  to allow automation of *find_centroid*
  - Use central atom as input in *find_centroid* (done)
  - Use ASE to combine *generate_tetrahedron.py* to get more substituents (done)
  - Use ASE to visualize created substituents (done, but make pictures better?)
  - Make code usable to create molecules like ethyl and methoxide, recursively (Done)
  - Use ASE to find distances between molecules (done)
  - convert .xyz to .gjf (optional, skipped for now)  
  - make python file usable with system arguments  
  - add bonded_atom and atom_to_be_functionalized to comment line of .xyz file  
  - make bash script that is able to use python file and deliver intermediates and final version of functionalized & optimized complex
  - Test test test
## Example usage
  - For the initial creation of a CH4 molecule from the CLI run: ```python main.py substituents_xyz/manually_generated/CH3.xyz CH4 [[1,0],[4,0]] initial H H C False False```
  - To create a C2H6 from this generated CH4 and view the molecule run: ```python main.py substituents_xyz/automatically_generated/CH4.xyz CH4 None recursive H H H True False```
## Authors
Adarsh Kalikadien & Vivek Sinha