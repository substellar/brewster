#!/usr/bin/env python

""" Module of bits to plug into Brewster """
import math
import gc
import numpy as np
import scipy as sp
import forwardmodel
from scipy import interpolate
from astropy.convolution import convolve, convolve_fft
from astropy.convolution import Gaussian1DKernel
#from pysynphot import observation
#from pysynphot import spectrum

__author__ = "Ben Burningham"
__copyright__ = "Copyright 2015 - Ben Burningham"
__credits__ = ["Ben Burningham","The EMCEE DOCS"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ben Burningham"
__email__ = "burninghamster@gmail.com"
__status__ = "Development"

# This bit is for flux conservation rebin of spectrum
#def rebinspec(wave, specin, wavenew,):
#    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
#    f = np.ones(len(wave))
#    filt = spectrum.ArraySpectralElement(wave, f, waveunits='microns')
#    obs = observation.Observation(spec, filt, binset=wavenew, force='taper')
 
#    return obs.binflux



def lnlike(intemp, invmr, pcover, cloudparams, r2d2, logg, dlam, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,logf):
    # get the ngas
    ngas = invmr.shape[0] + 1
    # Hard code nlayers
    nlayers = press.shape[0]
    # interp temp onto finer grid coarsePress => press
    # spline fit with max smoothing
    tfit = sp.interpolate.splrep(np.log10(coarsePress),intemp,s=0)
    temp = np.asfortranarray(sp.interpolate.splev(np.log10(press),tfit,der=0),dtype='d')
    # now loop through gases and get VMR for model
    # check if its a fixed VMR or a profile
    # VMR is log10(VMR) !!!
    logVMR = np.empty((ngas,nlayers),dtype='d')
    alkratio = 16.2 #  from Asplund et al (2009)
    if invmr.size > invmr.shape[0]:
        # now sort Na and K
        tmpvmr = np.empty((ngas,nlayers),dtype='d')
        tmpvmr[0:(ngas-2),:] = invmr[0:(ngas-2),:]
        tmpvmr[ngas-2,:] = np.log10(10.**invmr[ngas-2,:] / (alkratio+1.))
        tmpvmr[ngas-1,:] = np.log10(10.**invmr[ngas-2,:] * (alkratio / (alkratio+1.)))                                
        for i in range(0,ngas):
            vfit = sp.interpolate.splrep(np.log10(coarsepress),tmpvmr[i,:],s=0)
            logVMR[i,:] = sp.interpolate.splev(np.log10(press),vfit,der=0)
    else:
        # now sort Na and K
        tmpvmr = np.empty(ngas,dtype='d')
        tmpvmr[0:(ngas-2)] = invmr[0:(ngas-2)]
        tmpvmr[ngas-2] = np.log10(10.**invmr[ngas-2] / (alkratio+1.))
        tmpvmr[ngas-1] = np.log10(10.**invmr[ngas-2] * (alkratio / (alkratio+1.)))
        for i in range(0,ngas):                              
            logVMR[i,:] = tmpvmr[i]

    # now need to translate cloudparams in to cloud profile even
    # if do_clouds is zero..
    # 5 entries for cloudparams for simple slab model are:
    # 0) log10(number density / gas number density)
    # 1) top layer id (or pressure)
    # 2) base ID (these are both in 64 layers)
    # 3) rg
    # 4) rsig
    # in the case of a simple mixto cloud (i.e. cloudnum = 99) we have:
    # 0) ndens = dtau
    # 1) top layer ID
    # 2) bottom later ID
    # 3) rg = albedo
    # 4) rsig = asymmetry
    if (do_clouds == 1):
        npatch = cloudparams.shape[0]
        ncloud = cloudparams.shape[1]
        cloudrad = np.empty((npatch,nlayers,ncloud),dtype='d')
        cloudsig = np.empty_like(cloudrad)
        cloudprof = np.zeros_like(cloudrad)
        ndens= np.reshape(cloudparams['f0'],(npatch,ncloud))
        c1 = np.reshape(cloudparams['f1'],(npatch,ncloud))
        c2 = np.reshape(cloudparams['f2'],(npatch,ncloud))
        rad = np.reshape(cloudparams['f3'],(npatch,ncloud))
        sig = np.reshape(cloudparams['f4'],(npatch,ncloud))
        for i in range(0, npatch):
            for j in range(0, ncloud):
                b1 = c1[i,j] - 1
                b2 = c2[i,j] -1 
                cloudprof[i,b1:b2+1,j] = ndens[i,j]
                cloudrad[i,:,j] = rad[i,j]
                cloudsig[i,:,j] = sig[i,j]        
    else:
        npatch = 1
        ncloud = 1
        cloudrad = np.ones((npatch,nlayers,ncloud),dtype='d')
        cloudsig = np.ones_like(cloudrad)
        cloudprof = np.ones_like(cloudrad)

    # now we can call the forward model
    outspec = forwardmodel.marv(temp,logg,r2d2,gasnum,logVMR,pcover,do_clouds,cloudnum,cloudrad,cloudsig,cloudprof,inlinetemps,press,inwavenum,linelist,cia,ciatemps,use_disort)
    # Trim to length where it is defined.
    nwave = inwavenum.size
    trimspec = np.zeros((2,nwave),dtype='d')
    trimspec[:,:] = outspec[:,:nwave]
 
    # now shift wavelen by delta_lambda
    shiftspec = np.empty_like(trimspec)
    shiftspec[0,:] =  trimspec[0,:] + dlam
    shiftspec[1,:] =  trimspec[1,:]
 
    # length and interval for later
    wlen = shiftspec.shape[1]
    wint =  shiftspec[0,0] - shiftspec[0,wlen-1]

    # convolve with instrumental profile
    # start by setting up kernel
    # First step is finding the array index length of the FWHM
    disp = wint / wlen
    gwidth = int((((fwhm / disp) // 2) * 2) +1)

    # needs to be odd
    # now get the kernel and convolve
    gauss = Gaussian1DKernel(gwidth)
    cspec = convolve(shiftspec[1,:],gauss,boundary='extend')
    spec = np.array([shiftspec[0,::-1],cspec[::-1]])
    
    # rebin to observed dispersion
    wfit = sp.interpolate.splrep(spec[0,:],spec[1,:],s=0)
    modspec = sp.interpolate.splev(obspec[0,:],wfit,der=0)
    
    # Below is method for rebinning using conserve flux method
    #    oblen = obspec.shape[1]
    #    modspec = np.empty((2,oblen),dtype='d')
    #    modspec[1,:] =  rebinspec(spec[0,:], spec[1,:], obspec[0,:])

    # get log-likelihood
    # We've lifted this from Mike's code, below is original from emcee docs
    # Just taking every 3rd point to keep independence
    s2=obspec[2,::3]**2 #+ 10.**logf
    lnLik=-0.5*np.sum((((obspec[1,::3] - modspec[::3])**2) / s2) + np.log(2.*np.pi*s2))


    return lnLik
    #chi2 log likelihood--can modify this
    #invsigma2 = 1.0/((obspec[2,::3])**2 + modspec[1,::3]**2 * np.exp(2*lnf))
    #return -0.5*(np.sum((obspec[1,::3] - modspec[1,::3])**2 * invsigma2 - np.log(invsigma2)))
    
    
def lnprob(theta,pcover, cloudparams, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec):
    invmr = theta[0:7]
    logf = 0.0 #theta[5]
    #logbeta = theta[6]
    logg = theta[7]
    r2d2 = theta[8]
    dlam = theta[9]
    gam = theta[10]
    intemp = theta[11:]
    # now check against the priors, if not beyond them, run the likelihood
    lp = lnprior(theta,obspec)
    if not np.isfinite(lp):
        return -np.inf
    # else run the likelihood
    lnlike_value = lnlike(intemp, invmr,pcover, cloudparams,r2d2, logg, dlam, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,logf)

    lnprb = lp+lnlike_value
    if np.isnan(lnprb):
        lnprb = -np.inf
    return lnprb


def lnprior(theta,obspec):
    # set up the priors here
    invmr = theta[0:7]
    logf = 0.0 #theta[5]
    #logbeta = theta[6]
    logg = theta[7]
    r2d2 = theta[8]
    dlam = theta[9]
    gam = theta[10]
    T = theta[11:]

    diff=np.roll(T,-1)-2.*T+np.roll(T,1)
    pp=len(T)

    #for mass prior
    D = 3.086e+16 * 5.84
    R = np.sqrt(r2d2) * D
    g = (10.**logg)/100.
    M = (R**2 * g/(6.67E-11))/1.898E27
    
    #        ((0.001*np.max(obspec[2,:]**2)) < 10.**logf < (100.*np.max(obspec[2qq,:]**2))) and  and (-5. < logbeta < 0))
    if (all(invmr[0:7] > -12.0) and (np.sum(10.**(invmr[0:7])) < 1.0) 
        and  0.0 < logg < 6.0 
        and 1. < M < 80. 
        and  0. < r2d2 < 1. 
        and -0.01 < dlam < 0.01 
        and (min(T) > 1.0) and (max(T) < 5000.) 
        and (gam > 0.)): 
        logbeta = -5.0
    	beta=10.**logbeta
    	alpha=1.0
    	x=gam
    	invgamma=((beta**alpha)/math.gamma(alpha)) * (x**(-alpha-1)) * np.exp(-beta/x)
        prprob = (-0.5/gam)*np.sum(diff[1:-1]**2) - 0.5*pp*np.log(gam) + np.log(invgamma)
#        print -0.5*np.sum(diff[1:]**2/gam), -0.5*np.sum(np.log(2.*np.pi*gam)), np.log(invgamma)
        return prprob 
    return -np.inf


