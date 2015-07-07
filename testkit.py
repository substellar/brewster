#!/usr/bin/env python

""" Module of bits to plug into Brewster """

import numpy as np
import scipy as sp
import forwardmodel
from scipy import interpolate
from astropy.convolution import convolve, convolve_fft
from astropy.convolution import Gaussian1DKernel
from pysynphot import observation
from pysynphot import spectrum

__author__ = "Ben Burningham"
__copyright__ = "Copyright 2015 - Ben Burningham"
__credits__ = ["Ben Burningham","The EMCEE DOCS"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ben Burningham"
__email__ = "burninghamster@gmail.com"
__status__ = "Development"


def rebinspec(wave, specin, wavenew,):
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='microns')
    obs = observation.Observation(spec, filt, binset=wavenew, force='taper')
 
    return obs.binflux



def lnlike(w1,w2,intemp, invmr, pcover, cloudparams, r2d2, logg, dlam, do_clouds,gasnum,cloudnum,fwhm,obspec):
    # get the ngas
    ngas = invmr.shape[0]
    # interp temp onto finer grid (16 => 61)
    inlayer = np.arange(0,15.25,1)
    layer = np.arange(0,15.25,0.25)
    # Hard code nlayers
    nlayers = 61
    # spline fit with no smoothing)
    tfit = sp.interpolate.splrep(inlayer,intemp,s=0)
    temp = np.asfortranarray(sp.interpolate.splev(layer,tfit, der=0),dtype='f')

    # now loop through gases and get VMR for model
    # check if its a fixed VMR or a profile
    # VMR is log10(VMR) !!!
    logVMR = np.empty((ngas,nlayers),dtype='d')
    if invmr.size > invmr.shape[0]:
        for i in range(0,ngas):
            vfit = sp.interpolate.splrep(inlayer,invmr[i,:],s=0)
            logVMR[i,:] = sp.interpolate.splev(layer,vfit,der=0)
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

    outspec = forwardmodel.marv(w1,w2,temp,logg,r2d2,gasnum,logVMR,pcover,do_clouds,cloudnum,cloudrad,cloudsig,cloudprof)

    # Trim to length where it is defined.
    trimspec =  outspec[:,np.logical_not(np.logical_or(outspec[0,:] > w2, outspec[0,:] < w1))] 

    # now shift wavelen by delta_lambda
    shiftspec = np.array([trimspec[0,:],trimspec[0,:]+dlam])

 
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
    spec = np.array([shiftspec[0,:],cspec])
    
    # rebin to observed dispersion

    oblen = obspec.shape[1]
    modspec = np.empty((2,oblen),dtype='d')
    modspec[1,:] =  rebinspec(spec[0,:], spec[1,:], obspec[0,:])

    # get log-likelihood

    invsigma2 = 1/(obspec[2,:])**2
    return -0.5*(np.sum((obspec[1,:] - modspec[1,:])**2 * invsigma2 - np.log(2*np.pi*invsigma2)))
    
    
def lnprob(theta,fixvmrs,w1,w2,intemp, pcover, cloudparams, r2d2, logg, dlam, do_clouds,gasnum,cloudnum,fwhm,obspec):

    invmr = np.array([theta[0],theta[1], fixvmrs]).reshape(3,)
    # run the likelihood
    lnlike_value = lnlike(w1,w2,intemp, invmr,pcover, cloudparams, r2d2, logg, dlam, do_clouds,gasnum,cloudnum,fwhm,obspec)
    
    # now check against the priors
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike_value


def lnprior(theta):
    # set up the priors here
    invmr = theta
    if -4.5 < invmr[0] < -2.5 and -4.5 < invmr[1] < -2.5: 
        return 0.0
    return -np.inf


