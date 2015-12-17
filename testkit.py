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



def lnlike(w1,w2,intemp, invmr, pcover, cloudparams, r2d2, logg, dlam, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,logf):
    # get the ngas
    ngas = invmr.shape[0]
    # Hard code nlayers
    nlayers = press.shape[0]
    # interp temp onto finer grid coarsePress => press
    # spline fit with max smoothing
    tfit = sp.interpolate.splrep(np.log10(coarsePress),np.log10(intemp),s=10)
    temp = 10.**(np.asfortranarray(sp.interpolate.splev(np.log10(press),tfit,der=0),dtype='d'))
    # now loop through gases and get VMR for model
    # check if its a fixed VMR or a profile
    # VMR is log10(VMR) !!!
    logVMR = np.empty((ngas,nlayers),dtype='d')
    if invmr.size > invmr.shape[0]:
        for i in range(0,ngas):
            vfit = sp.interpolate.splrep(np.log10(coarsepress),invmr[i,:],s=0)
            logVMR[i,:] = sp.interpolate.splev(np.log10(press),vfit,der=0)
    else:
        for i in range(0,ngas):
            logVMR[i,:] = invmr[i]

    # now need to translate cloudparams in to cloud profile even
    # if do_clouds is zero..
    # 5 entries for cloudparams for simple slab model are:
    # 0) log10(number density)
    # 1) top layer id (or pressure)
    # 2) base ID (these are both in 61 layers)
    # 3) rg
    # 4) rsig
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
    outspec = forwardmodel.marv(w1,w2,temp,logg,r2d2,gasnum,logVMR,pcover,do_clouds,cloudnum,cloudrad,cloudsig,cloudprof,inlinetemps,press,inwavenum,linelist,cia,ciatemps,use_disort)

    # Trim to length where it is defined.
    trimspec =  outspec[:,np.logical_not(np.logical_or(outspec[0,:] > w2, outspec[0,:] < w1))] 

    # now shift wavelen by delta_lambda
    shiftspec = np.array([trimspec[0,:]+dlam,trimspec[1,:]])

 
    # length and interval for later
    wlen = shiftspec.shape[1]
    wint =  shiftspec[0,0] - shiftspec[0,wlen-1]
    
    # convolve with instrumental profile
    # start by setting up kernel
    # First step is finding the array index length of the FWHM
    disp = wint / wlen
    gwidth = int(round(fwhm / disp))
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
    s2=obspec[2,5::3]**2 #+ 10.**logf
    lnLik=-0.5*np.sum((((obspec[1,5::3] - modspec[5::3])**2) / s2) + np.log(2.*np.pi*s2))


    return lnLik
    #chi2 log likelihood--can modify this
    #invsigma2 = 1.0/((obspec[2,::3])**2 + modspec[1,::3]**2 * np.exp(2*lnf))
    #return -0.5*(np.sum((obspec[1,::3] - modspec[1,::3])**2 * invsigma2 - np.log(invsigma2)))
    
    
def lnprob(theta,w1,w2,pcover, cloudparams, r2d2,logg, dlam, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec):
    invmr = theta[0:5]
    logf = 0.0 #theta[5]
    gam = theta[5]
#    logbeta = theta[6]    
    intemp = theta[6:]
    
    # now check against the priors, if not beyond them, run the likelihood
    lp = lnprior(theta,obspec)
    if not np.isfinite(lp):
        return -np.inf
    # else run the likelihood
    lnlike_value = lnlike(w1,w2,intemp, invmr,pcover, cloudparams, r2d2, logg, dlam, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,logf)
    lnprb = lp+lnlike_value
    if np.isnan(lnprb):
        lnprb = -np.inf
    return lnprb


def lnprior(theta,obspec):
    # set up the priors here
    invmr = theta[0:5]
#    logf = theta[5]
    gam = theta[5]
#    logbeta = theta[6]
    T = theta[6:]
    diff=np.roll(T,-1)-2.*T+np.roll(T,1)
    pp=len(T)

    if (all(invmr[0:5] > -12.0) and (np.sum(10.**(invmr[0:5])) < 1.0) and (min(T) > 10.0) and (max(T) < 4000) and (gam > 0)):
#        ((0.001*np.max(obspec[2,:]**2)) < 10.**logf < (100.*np.max(obspec[2,:]**2))) and  and (-5. < logbeta < 0))
        logbeta = -5.0
    	beta=10.**logbeta
    	alpha=1.0
    	x=gam
    	invgamma=((beta**alpha)/math.gamma(alpha)) * (x**(-alpha-1)) * np.exp(-beta/x)
        prprob = (-0.5/gam)*np.sum(diff[1:-1]**2) - 0.5*pp*np.log(gam) + np.log(invgamma)
#        print -0.5*np.sum(diff[1:]**2/gam), -0.5*np.sum(np.log(2.*np.pi*gam)), np.log(invgamma)
        return prprob 
    return -np.inf


