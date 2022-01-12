#!/usr/bin/env python

""" Module of bits for setting up T-P profile for Brewster"""
import math
import gc
import numpy as np
import scipy as sp
from scipy import interpolate
from astropy.convolution import convolve, Box1DKernel,Gaussian1DKernel


__author__ = "Ben Burningham"
__copyright__ = "Copyright 2016 - Ben Burningham"
__credits__ = ["Ben Burningham"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ben Burningham"
__email__ = "burninghamster@gmail.com"
__status__ = "Development"



# This code takes in parameters for T-P profile and spit out the profile
# Type 1 is the knots for a spline
# Type 2 is a Madhusudhan & Seager 2009 parameterised profile, no inversion
# i.e. T1,P1 == T2,P2
# Type 3 is Madhusudhan & Seager 2009 with an inversion

def set_prof(proftype, coarsePress,press,intemp):

    if (proftype == 1 or proftype == 9):
        # interp temp onto finer grid coarsePress => press
        # spline fit with max smoothing
        tfit = sp.interpolate.splrep(np.log10(coarsePress),intemp,s=0)
        # we'll take the abs value to avoid negative T
        temp = np.asfortranarray(np.abs(sp.interpolate.splev(np.log10(press),tfit,der=0)),dtype='d')
        
        #tfit = sp.interpolate.UnivariateSpline(np.log10(coarsePress),intemp,k=2,s=0)
        #temp = tfit(np.log10(press))
                    
    elif (proftype == 2):
        # unpact the parameters
        beta = 0.5
        P0 = press[0]
        #  We can either have T0 or T3 as a parameter. We're taking T3
        a1, a2, logP1,logP3,T3 = intemp[0:5]

        P1 = 10.**logP1
        P3 = 10.**logP3


        # set T1 from T3
        T1 = T3 - (np.log(P3/P1) / a2)**(1/beta)
    
        # Set T0 from T1
        T0 = T1 - ((np.log(P1/P0) / a1)**(1/beta))

         
        temp = np.zeros_like(press)
        
        for i in range(0,press.size):
            if (press[i] < P1):
                temp[i] = T0 + (np.log(press[i] / P0) / a1)**(1/beta)
            elif (press[i] >= P1 and press[i] < P3):
                temp[i] = T1 + (np.log(press[i] / P1) / a2)**(1/beta)
            elif (press[i] >= P3):
                temp[i] = T3

        # then smooth with 5 layer box car
        temp1 = convolve(temp,Gaussian1DKernel(5),boundary='extend')
        temp = temp1

    elif (proftype == 3):
        # unpact the parameters
        beta = 0.5
        P0 = press[0]
        #  We can either have T0 or T3 as a parameter. We're taking T3
        a1, a2, logP1,logP2,logP3,T3 = intemp[0:6]

        P1 = 10.**logP1
        P2 = 10.**logP2
        P3 = 10.**logP3
        
        # set up Ts at boundaries from continuity

        T2 = T3 - (np.log(P3/P2) / a2)**(1/beta)
        T1 = T2 + (np.log(P1/P2) / a2)**(1/beta)
        T0 = T1 - (np.log(P1/P0) / a1)**(1/beta)

        
        temp = np.zeros_like(press)
        
        for i in range(0,press.size):
            if (press[i] < P1):
                temp[i] = T0 + (np.log(press[i] / P0) / a1)**(1/beta)
            elif (press[i] >= P1 and press[i] < P3):
                temp[i] = T2 + (np.log(press[i] / P2) / a2)**(1/beta)
            elif (press[i] >= P3):
                temp[i] = T3

        # then smooth with 5 layer box car 
        temp1 = convolve(temp,Gaussian1DKernel(5),boundary='extend')
        temp = temp1

        
    return temp
