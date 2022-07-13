# ChemSpaX in an automated mode
To use ChemSpaX in an automated mode please follow the following steps:  
1. Install ChemSpaX in the same way as described in the regular README.  
2. Place all the skeletons with the functionalisation list in the skeletons folder.  
3. Make sure your wished substituents are placed in the substituents\_xyz folder in both XYZ and MDL Molfile format with the atom to be connected to as the first atom in the XYZ file.  
4. Run data\_preperation.py.  
5. Generate the right list of substitutes through the substitutes\_combined.py in the substitutes\_list folder. In case of many substituents a slurmfile is accompanied to generate the list of substituents on a supercomputer.  
6. Use the gen\_chemspax.sh bash script in the command line with the arguments:  
	* skeleton\_name (this name is used to keep track of the logfiles and output files)  
	* number of substituent sites (this way the right substituents file is taken, make sure you generated this file in step 3)

Have fun with your structures!