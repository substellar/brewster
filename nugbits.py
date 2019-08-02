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
            
    for j in range(1, (wave.size - 1)):
        sbin = ((wave[j] - wave[j-1]) + (wave[j+1] - wave[j])) / 2. 
        fbol = (sbin * flux[j]) + fbol
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

    result = np.concatenate((theta,np.array([t_ff, R, M])),axis=0)
    
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
        stop
        
    return flatendchain, flatendprobs,ndim

def getargs(runname):
    
    pic = runname+"_runargs.pic"
    with open(pic, 'rb') as input:
        gases_myP,chemeq,dist, cloudtype,do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge, prof,do_bff,bff_raw,ceTgrid,metscale,coscale = pickle.load(input) 

    args =  gases_myP,chemeq,dist, cloudtype,do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge, prof,do_bff,bff_raw,ceTgrid,metscale,coscale
    return args
