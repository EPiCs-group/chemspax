# -*- coding: utf-8 -*- 
#                                                                                   #
#  __authors__ = Mark Heezen                                                        #
#  __institution__ = Vrije Universiteit Brussel & Universidad Autónoma de Madrid    #
#  

"""This file contains the functions to choose between the possible subfunctionalisations"""

def H_alkyl_aryl(n):
    if n in [1,2,3]:
        return 'H'
    elif n in [4,5]:
        return 'CH3'
    elif n==6:
        return 'CH2CH3'
    else:
        ValueError('Wrong input in function H_alkyl_aryl')

def alkyl_aryl(n):
    if n in [1,2,3]:
        return 'CH3'
    elif n in [4,5]:
        return 'CH2CH3'
    elif n==6:
        return 'C6H6'
    else:
        ValueError('Wrong input in function alkyl_aryl')

def alkyl(n):
    if n in [1,2,3]:
        return 'CH3'
    elif n==[4,5]:
        return 'CH2CH3'
    elif n==6:
        return 'CH2CH2CH3'
    else:
        ValueError('Wrong input in function alkyl')

def aryl(n):
    return 'C6H6'

def heteroatom(n):
    if n in [1,2,3]:
        return 'NH2'
    elif n in [4,5]:
        return 'OCH3'
    elif n==6:
        return 'OCH2CH3'
    else:
        ValueError('Wrong input in function heteroatom')