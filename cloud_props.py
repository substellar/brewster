import math
import time
import gc
import numpy as np
import scipy as sp
import cloud
import TPmod
import cloudpost
import settings
import sys
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import interp1d
from astropy.convolution import convolve, convolve_fft
from astropy.convolution import Gaussian1DKernel
from bensconv import prism_non_uniform
from bensconv import conv_uniform_R
from bensconv import conv_uniform_FWHM


__author__ = "Ben Burningham"
__copyright__ = "Copyright 2021 - Ben Burningham"
__credits__ = ["Ben Burningham"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ben Burningham"
__email__ = "burninghamster@gmail.com"
__status__ = "Development"


def get(theta,chemeq,gasnum,fwhm,do_fudge,cloudtype,do_clouds,cloudnum,coarsePress,press,inwavenum):

    # This returns cloud column density in g / cm2 / layer
    # and cloud column number density in /cm2 / layer
    # and optical depth of each cloud by layer
    nlayers = press.size
    if chemeq == 0:
        if (gasnum[gasnum.size-1] == 21):
            ng = gasnum.size - 1
        elif (gasnum[gasnum.size-1] == 23):
            ng = gasnum.size -2
        else:
            ng = gasnum.size
    else:
        ng = 2

    if (fwhm < 0.0):
        if (fwhm == -1 or fwhm == -3 or fwhm == -4):
            if (do_fudge == 1):
                logf = theta[ng+5:ng+8]
                nb = 8
            else:
                nb = 5
        elif (fwhm == -2):
            if (do_fudge == 1):
                nb = 6
            else:
                nb = 4

    else:
        if (do_fudge == 1):
            nb = 4
        else:
            nb = 3

        
    npatches = do_clouds.size
    if (npatches > 1):
        prat = theta[ng+nb]
        pcover = np.array([prat,(1.-prat)])
        pc = ng + nb + 1
    else:
        pc = ng + nb 
        pcover = 1.0
        
    # use correct unpack method depending on situation
    
    if ((npatches > 1) and np.all(do_clouds != 0)):
        cloudparams, nc = cloud.unpack_patchy(theta,pc,cloudtype,cloudnum,do_clouds)
    else:
        cloudparams, nc = cloud.unpack_default(theta,pc,cloudtype,cloudnum,do_clouds)
    
    # now need to translate cloudparams in to cloud profile even
    # if do_clouds is zero..

    cloudprof,cloudrad,cloudsig = cloud.atlas(do_clouds,cloudnum,cloudtype,cloudparams,press)
    cloudprof = np.asfortranarray(cloudprof,dtype = 'float64')
    cloudrad = np.asfortranarray(cloudrad,dtype = 'float64')
    cloudsig = np.asfortranarray(cloudsig,dtype = 'float64')
    pcover = np.asfortranarray(pcover,dtype = 'float32')
    cloudnum = np.asfortranarray(cloudnum,dtype='i')
    do_clouds = np.asfortranarray(do_clouds,dtype = 'i')

    nclouds = cloudtype.size


    tau1_cloud,mass_cloud,num_cloud = cloudpost.properties(press,inwavenum,nclouds,do_clouds,cloudnum,cloudprof,cloudrad,cloudsig)
    
    
    
    #trim these to lengths where they have data
    nwave = inwavenum.size
    npatch = do_clouds.size
    nlayers = press.size

    taucloud1_press = tau1_cloud[0:npatch,0:nwave,0:nclouds]
    cloudmass = mass_cloud[0:npatch,0:nlayers,0:nclouds]
    cloudnumdens =  num_cloud[0:npatch,0:nlayers,0:nclouds]

    return taucloud1_press,cloudmass,cloudnumdens
