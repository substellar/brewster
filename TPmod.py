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
# Type 7 is Molliere+2020,2021 w/o inversion, and dry adiabat (instead of moist)


def set_prof(proftype, coarsePress,press,intemp):
    temp = np.zeros_like(press)

    if (proftype == 1 or proftype == 9 or proftype == 6):
        # interp temp onto finer grid coarsePress => press
        # spline fit with max smoothing
        tfit = sp.interpolate.splrep(np.log10(coarsePress),intemp,s=0)
        # we'll take the abs value to avoid negative T
        temp1 = np.asfortranarray(np.abs(sp.interpolate.splev(np.log10(press),tfit,der=0)),dtype='d')
        
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

         
        
        for i in range(0,press.size):
            if (press[i] < P1):
                temp[i] = T0 + (np.log(press[i] / P0) / a1)**(1/beta)
            elif (press[i] >= P1 and press[i] < P3):
                temp[i] = T1 + (np.log(press[i] / P1) / a2)**(1/beta)
            elif (press[i] >= P3):
                temp[i] = T3

        # then smooth with 5 layer box car
        temp1 = convolve(temp,Gaussian1DKernel(5),boundary='extend')
 
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

        for i in range(0,press.size):
            if (press[i] < P1):
                temp[i] = T0 + (np.log(press[i] / P0) / a1)**(1/beta)
            elif (press[i] >= P1 and press[i] < P3):
                temp[i] = T2 + (np.log(press[i] / P2) / a2)**(1/beta)
            elif (press[i] >= P3):
                temp[i] = T3

        # then smooth with 5 layer box car 
        temp1 = convolve(temp,Gaussian1DKernel(5),boundary='extend')
 
    elif (proftype == 7):
        # this is Molliere's hybrid profile, hacked by Michelle Colantoni from
        # petitRadTran. But, using dry adiabat for H2/He atmosphere
        # a few variable names changed by BB, and a bit of restructuring

        # unpack the intemp parameters
        Tint,alpha,delta,T1,T2,T3 = intemp[0:6]
        
        n=len(press)

        # convert pressure to tau using free parameters alpha, delta
        tau=delta*(press)**alpha

        # Get the ratio of specific heats for the dry adiabat
        cp = 0.84*14.32 + 0.16*5.19
        cv = 0.84*10.16 + 0.16*3.12
        gamma=cp/cv

        # first the get T for top layers (stratosphere) where tau < 0.1
        Tconnect = (((3/4) * Tint**4) * ((2/3) + (0.1)))**(1/4)
        t_support=np.array([T1, T2, T3, Tconnect])
        Pstrat = ((0.1/delta)**(1/alpha))
        # prevent cubic spline fart due to decreasing values in P 
        if (Pstrat < press[0]):
            Pstrat = press[1]
        support_points=np.linspace(np.log10(press[0]), np.log10(Pstrat),4) 
        values=interpolate.CubicSpline(support_points, t_support) # fit it in log10(P)
        toplayers = np.where(press <= Pstrat)

        temp[toplayers] = values(np.log10(press[toplayers]))
        # Now get the temperature in the radiative zone, following eddington

        T_edd=(((3/4)*Tint**4)*((2/3)+(tau)))**(1/4)

        # set all temps deeper that Pstrat to T_edd, before replacing with
        # adiabat below RCbound
        Tedd_top = toplayers[-1][-1]
        temp[Tedd_top:] = T_edd[Tedd_top:]


        # Now get the boundary where the radiative gradient exceeds the adiabat

        nabla_ad=(gamma-1)/gamma
        nabla_rad = np.diff(np.log(T_edd))/np.diff(np.log(press))
        convtest = np.any(np.where(nabla_rad >= nabla_ad))

        # Now get temperatures on the adiabat from RC boundary downwards
        if convtest:            
            RCbound = np.where(nabla_rad >= nabla_ad)[0][0]
            temp[RCbound:] = temp[RCbound]*(press[RCbound:]/press[RCbound])**((gamma-1)/gamma)
        # weed out any temperatures that are above our opacity tables 
        temp2 = np.where(temp < 6000., temp, 6000.)
                        
        # finally smooth with a 5 layer boxcar
        temp1 = convolve(temp2,Gaussian1DKernel(5),boundary='extend')
 
    # weed out any negative temperatures
    temp = np.where(temp1 > 10., temp1, 10.)

    return temp
