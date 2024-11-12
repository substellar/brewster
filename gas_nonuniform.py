#!/usr/bin/env python

""" Module of bits for non-uniform gas profile"""
import numpy as np

__author__ = "Fei Wang"
__copyright__ = "Copyright 2024 - Fei Wang"
__credits__ = ["Fei Wang", "Ben Burningham"]
__license__ = "GPL"
__version__ = "0.2"  
__maintainer__ = ""
__email__ = ""
__status__ = "Development"


def non_uniform_gas(press,logPt,logft,alpha):
    gas_f = np.zeros_like(press)
    
    Pt = 10.**logPt
    ft = 10.**logft
    
    for i in range(0,press.size):
        if (press[i] < Pt):
            gas_f[i]=logft+(np.log10(press[i])-logPt) / alpha
        else:
            gas_f[i]=float(logft)
    
    return gas_f


