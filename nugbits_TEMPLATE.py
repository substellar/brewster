#!/usr/bin/env python


""" Bits for McNuggets: the post-processing tool for brewster"""
from __future__ import print_function

from builtins import str
from builtins import range
import numpy as np
import scipy as sp
import testkit
import ciamod
import TPmod
import settings
import os
import gc
import sys
import pickle
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from mpi4py import MPI
  
def teffRM(theta,sigDist,sigPhot):



    dist = settings.runargs[2]
    chemeq = settings.runargs[1]
    fwhm = settings.runargs[15]
    gasnum = settings.runargs[5]
    gnostics = 0
    shiftspec, photspec,tauspec,cfunc = testkit.modelspec(theta,settings.runargs,gnostics)  
    wave = np.array(shiftspec[0,::-1])
    outflux = np.array(shiftspec[1,::-1])
    
    # now calculate Fbol by summing the spectrum across its wave bins
    fbol = 0.0

    if chemeq == 0:
        if (gasnum[gasnum.size-1] == 21):
            ng = gasnum.size - 1
        elif (gasnum[gasnum.size-1] == 23):
            ng = gasnum.size -2
        else:
            ng = gasnum.size
    else:
        ng = 2


    logg = theta[ng]
    if (fwhm < 0.0):
        if (fwhm == -1 or fwhm == -3):
            R2D2 = theta[ng+1:ng+4]
        elif (fwhm == -2):
            R2D2 = theta[ng+1:ng+3]
    else:
        r2d2 = theta[ng+1]
        
    # get r2d2 sorted for multi-instruments
    # For now we'll just use the r2d2 for the NIR for the postprocessing
    if (fwhm < 0.0):
        if (fwhm == -1 or fwhm == -3):
            r2d2 = R2D2[0]
            scale1 = R2D2[1]
            scale2 = R2D2[2]
        elif (fwhm == -2):
            r2d2 = R2D2[0]
            scale1 = R2D2[1]            
    else:
        r2d2 = theta[ng+1]

    if (fwhm < 0.0):
        # This is for multi-instrument cases
        # -1: spex + akari + IRS
        # -2: spex + IRS
        # -3: spex + Lband + IRS
        if (fwhm == -1 or fwhm == -3):

            # Spex
            mr1 = np.where(wave < 2.5)
            spec1 = outflux[mr1]

            # AKARI IRC or L' band
            mr2 = np.where(np.logical_and(wave > 2.5,wave < 5.0))
            spec2 = scale1 * outflux[mr2]

            # Spitzer IRS
            mr3 = np.where(wave > 5.0)
            spec3 = scale2 * outflux[mr3]

            flux = np.concatenate((spec1,spec2,spec3))
        elif (fwhm == -2):
            # Spex
            mr1 = np.where(wave < 2.5)
            spec1 = outflux[mr1]

            # Spitzer IRS
            mr2 = np.where(wave > 5.0)
            spec2 = scale1 * outflux[mr2]

            flux = np.concatenate((spec1,spec2))
    else:
        # For SpeX only data
        flux = outflux
            
    for j in range(1, (wave.size - 1)):
        sbin = ((wave[j] - wave[j-1]) + (wave[j+1] - wave[j])) / 2. 
        fbol = (sbin * flux[j]) + fbol

    # Get Lbol in log(L / Lsol)
    lsun = 3.828e26
    l_bol = np.log10((fbol * 4.*np.pi*(dist * 3.086e16)**2) / lsun)

    # now get T_eff
    t_ff = ((fbol/(r2d2 * 5.670367e-8))**(1./4.))

    # and Radius
    sigR2D2 = sigPhot * r2d2 * (-1./2.5)* np.log(10.)

    sigD = sigDist * 3.086e16
    D = dist * 3.086e16

    R = np.sqrt(((np.random.randn() * sigR2D2)+ r2d2)) \
        * ((np.random.randn()* sigD) + D)

    g = (10.**logg)/100.

    # and mass

    M = (R**2 * g/(6.67E-11))/1.898E27
    R = R / 71492e3

    # Now lets get the C/O and M/H ratios...
    # read molecules (elements!) to be used for M/H and C/O from theta
    # Think about these choices, and maybe experiment
    # Are the elements depleted by missed gases, or condensation
    # e.g. Fe/H, from retrieved FeH abundance will look subsolar
    # due to condensation of Fe.
    # Similarly, N2 is not observable, so N/H from NH3 may look subsolar...

    # Example here is for T dwarf gaslist:['h2o','ch4','co','co2','nh3','h2s','k','na']

    # We will base e
    h2o = theta[0]
    ch4 = theta[1]
    co = theta[2]
    co2 = theta[3]
    nak = theta[7]

    # first get the C/O

    O = 10**h2o + 10**co + 2.*10**(co2) 
    C = 10**(co) + 10**(co2) + 10**(ch4)

    CO_ratio = C/O

    # rest of the elements
    NaK = 10**nak

    # Determine "fraction" of H2 in the L dwarf
    gas_sum = 10**h2o + 10**co +10**co2 + 10**ch4 + 10**nak
    fH2 = (1-gas_sum)* 0.84  # fH2/(fH2+FHe) = 0.84
    fH = 2.*fH2

    
    # Determine linear solar abundance sum of elements in our L dwarf
    # abundances taken from Asplund+ 2009
    solar_H = 12.00
    solar_O = 10**(8.69-solar_H)
    solar_C = 10**(8.43-solar_H) 
    #solar_Ti = 10**(4.95-solar_H)
    #solar_V = 10**(3.93-solar_H)
    #solar_Cr = 10**(5.64-solar_H)
    #solar_Fe = 10**(7.50-solar_H)
    solar_NaK = 10**(6.24-solar_H) + 10**(5.03-solar_H)
    
    # Calculate the metallicity fraction in the star and the same for the sun and then make the ratio
    metallicity_target = (O/fH) + (C/fH) + (NaK/fH)
    metallicity_sun = solar_O + solar_C + solar_NaK

    MH = np.log10(metallicity_target / metallicity_sun)

    result = np.concatenate((theta,np.array([l_bol, t_ff, R, M, MH, CO_ratio])),axis=0)
    
    return result

def get_endchain(runname,fin):
    if (fin == 1):
        pic = runname+".pk1"
        with open(pic, 'rb') as input:
            sampler = pickle.load(input) 
        nwalkers = sampler.chain.shape[0]
        niter = sampler.chain.shape[1]
        ndim = sampler.chain.shape[2]
        flatprobs = sampler.lnprobability[:,:].reshape((-1))
        flatendchain = sampler.chain[:,niter-2000:,:].reshape((-1,ndim))
        flatendprobs = sampler.lnprobability[niter-2000:,:].reshape((-1))
 

    elif(fin ==0):
        pic = runname+"_snapshot.pic"
        with open(pic, 'rb') as input:
            chain,probs = pickle.load(input) 
        nwalkers = chain.shape[0]
        ntot = chain.shape[1]
        ndim = chain.shape[2]
        niter = int(np.count_nonzero(chain) / (nwalkers*ndim))
        flatprobs = probs[:,:].reshape((-1))
        flatendchain = chain[:,(niter-2000):niter,:].reshape((-1,ndim))
        flatendprobs = probs[(niter-2000):niter,:].reshape((-1))
    else:
        print("File extension not recognised")

        
    return flatendchain, flatendprobs,ndim

def getargs(runname):
    
    pic = runname+"_runargs.pic"
    with open(pic, 'rb') as input:
        gases_myP,chemeq,dist, cloudtype,do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge, prof,do_bff,bff_raw,ceTgrid,metscale,coscale = pickle.load(input) 

    args =  gases_myP,chemeq,dist, cloudtype,do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge, prof,do_bff,bff_raw,ceTgrid,metscale,coscale
    return args
