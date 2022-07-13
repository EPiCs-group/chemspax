# -*- coding: utf-8 -*-
#                                                     #
#  __authors__ = Mark Heezen                          #
#  __institution__ = TU Delft                         #
#                                                     #
"""This file is used to make lists of substituents combinatorically. 

This way many functionalisations can be made using gen_chemspax.sh script.
"""


import itertools
import sys

n = sys.argv[1] #Number of substitution sites

# Define all substitutes that you want to functionalise in a list with the names in the substitutes folder from chemspax
substitutes = ["H", "CH3", "CH2CH3", "CHCH3CH3", "CCH3CH3CH3", "OCCH3CH3CH3", "OCHCH3CH3", "C6H6", "C6H6-CH3-ortho-1-2-para", "C6H6-CH3-meta-1-2", "C6H6-CH3-para",  "C6H6-iPr-ortho-1-2-para", "C6H6-iPr-ortho-1-2-CH3-para"]
temp = itertools.product(substitutes, repeat=n-1) #n-1 since the first substitute has to be a seperate argument for ChemSpaX

# writing all the substitute combinations to a txt file
name = "substitutes_combined_%s.txt" %(n)

options = []
for i in temp:
    options.append(i)


with open(name, 'w') as f:
    for j in substitutes:
        for i in options:
            k = j + " " + '\"' + ' '.join(i) + '\"\n'
            f.write(k)
