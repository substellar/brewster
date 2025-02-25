#!/usr/bin/env python

""" Module of bits to process mags to fluxes and spectra to band fluxes """
from __future__ import print_function

import scipy as sp
import numpy as np
from scipy import interpolate
from scipy.interpolate import interp1d




__author__ = "Ben Burningham"
__copyright__ = "Copyright 2015 - Ben Burningham"
__credits__ = ["Ben Burningham", "The EMCEE DOCS"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ben Burningham"
__email__ = "burninghamster@gmail.com"
__status__ = "Development"

""" This is all about taking a magnitude and getting total flux through a filter
or taking a spectrum and getting the same. 

All in Vega until further notice. 

Vega fluxes are hardwired for speed. 

Units are W/m2/um

Where not calculated WISE details are taken from Jarrett et al (2011)

"""


def getfilt(filtname):
    
    if filtname == "Jmko":
        rawfilt = np.loadtxt("data/UKIRT-UKIDSS.J.dat", unpack=True, skiprows=0)
        tempfilt = rawfilt[1, :] * rawfilt[0, :]
        rawfilt[1, :] = tempfilt / np.amax(tempfilt)
    elif filtname == "Hmko":
        rawfilt = np.loadtxt("data/UKIRT-UKIDSS.H.dat", unpack=True, skiprows=0)
        tempfilt = rawfilt[1, :] * rawfilt[0, :]
        rawfilt[1, :] = tempfilt / np.amax(tempfilt)
    elif filtname == "nirc_Lp":
        rawfilt = np.loadtxt("data/nirc_Lp.txt", unpack=True, skiprows=1)
        bw = 0.700
        isow = 3.776
        tempfilt = rawfilt[1, :] * rawfilt[0, :]
        rawfilt[1, :] = tempfilt / np.amax(tempfilt)
    elif filtname == "w1":
        rawfilt = np.loadtxt("data/RSR-W1.EE.txt", unpack=True, skiprows=0)
        bw = 6.6256e-01
        isow = 3.35
    elif filtname == "w2":
        rawfilt = np.loadtxt("data/RSR-W2.EE.txt", unpack=True, skiprows=0)
        bw = 1.0423
        isow = 4.60
    elif filtname == "w3":
        rawfilt = np.loadtxt("data/RSR-W3.EE.txt", unpack=True, skiprows=0)
        bw = 5.5069
        isow = 11.56
    elif filtname == "w4":
        rawfilt = np.loadtxt("data/RSR-W4.EE.txt", unpack=True, skiprows=0)
        bw = 4.1013
        isow = 22.08
    else:
        print("Filter ", filtname, " not recognised")
        return np.nan

    return rawfilt, bw, isow


def mag2flux(mag, magerr, filtname, iso=False):
 
    rawfilt, bw, isow = getfilt(filtname)
    
    # First trim vega to match filter, and put filter on same wave grid as vega

    rawvega = np.loadtxt("data/STSci_Vega.txt", unpack=True, skiprows=0)

    rawvega[0, :] = rawvega[0, :] / 10000.
    rawvega[1, :] = rawvega[1, :] * 10.  # erg/cm/s/A to W/m2/um

    w1 = rawfilt[0, 0]
    w2 = rawfilt[0, rawfilt.shape[1] - 1]

    vega = rawvega[:, np.logical_not(np.logical_or(rawvega[0, :] > w2, rawvega[0, :] < w1))]

    filt = np.zeros_like(vega)
    wfit = sp.interpolate.splrep(rawfilt[0, :], rawfilt[1, :], s=0)
    filt[0, :] = vega[0, :]
    filt[1, :] = sp.interpolate.splev(filt[0, :], wfit, der=0)

    # Now we'll sum flux across the filter
    # need the bin width first
    sbin = np.zeros(filt.shape[1])
    for i in range(1, filt.shape[1]-1):
        sbin[i] = ((vega[0, i] - vega[0, i-1]) + (vega[0, i+1] - vega[0, i])) / 2.
    # Deal with the ends...
    sbin[0] = sbin[1]
    sbin[sbin.size-1] = sbin[sbin.size-2]

    bandvega = sum(sbin*filt[1, :]*vega[1, :])

    bandflux = bandvega * 10.**(-mag/2.5)

    banderr = abs(bandflux * (-1./2.5) * np.log(10) * magerr)

    if iso:
        isoflux = bandflux / bw
        isofluxerr = banderr / bw
        return isow, isoflux, isofluxerr

    return bandflux, banderr


def spec2flux(rawspec, filtname, iso=False):

    rawfilt, bw, isow = getfilt(filtname)

    # Now trim target spectrum to wavelength region of filter, and rebin filter onto new target grid
    w1 = rawfilt[0, 0]
    w2 = rawfilt[0, rawfilt.shape[1] - 1]
    trimspec = rawspec[:, np.logical_not(np.logical_or(rawspec[0, :] > w2, rawspec[0, :] < w1))]
    filt = np.zeros_like(trimspec)
    filt[0, :] = trimspec[0, :]
    wfit = sp.interpolate.splrep(rawfilt[0, :], rawfilt[1, :], s=0)
    filt[1, :] = sp.interpolate.splev(filt[0, :], wfit, der=0)

    # Now we want to sum over the filter
    sbin = np.zeros(filt.shape[1])
    for i in range(1, filt.shape[1]-1):
        sbin[i] = ((trimspec[0, i] - trimspec[0, i-1]) + (trimspec[0, i+1] - trimspec[0, i])) / 2.
    # Deal with the ends...
    sbin[0] = sbin[1]
    sbin[sbin.size-1] = sbin[sbin.size-2]
    
    bandflux = sum(sbin*filt[1, :]*trimspec[1, :])

    if iso:
        isoflux = bandflux / bw
        return isow, isoflux

    return bandflux
