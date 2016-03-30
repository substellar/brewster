#!/usr/bin/env python

""" Module of processes to interpret cloud parameters from Brewster in testkit"""

import numpy as np
import scipy as sp
from scipy import interpolate
from astropy.convolution import convolve, convolve_fft
from astropy.convolution import Gaussian1DKernel


__author__ = "Ben Burningham"
__copyright__ = "Copyright 2016 - Ben Burningham"
__credits__ = ["Ben Burningham","The EMCEE DOCS"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ben Burningham"
__email__ = "burninghamster@gmail.com"
__status__ = "Development"


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



def atlas(do_clouds,cloudnum,cloudtype,cloudparams,press):

    # Cloud types
    # 1:  slab cloud
    # 2: deep thick cloud , we only see the top
    # In both cases the cloud properties are density, rg, rsig for real clouds
    # and dtau, w0, and gg or cloudnum = 99
    nlayers = press.size
    npatch = 1 #cloudparams.shape[0]
    ncloud = 1 # cloudparams.shape[1]
    cloudrad = np.zeros((npatch,nlayers,ncloud),dtype='d')
    cloudsig = np.zeros_like(cloudrad)
    cloudprof = np.zeros_like(cloudrad)
    
    for i in range(0, npatch):
        if (do_clouds[i] == 1 and cloudtype[i] == 1):
            # 5 entries for cloudparams are:
            # 0) density
            # 1) log top pressure
            # 2) pressure thickness in dex
            # 3) rg
            # 4) rsig
            


            ndens= cloudparams[0]
            p1 = 10.**cloudparams[1]
            p2 = p1 * 10.**cloudparams[2]
            rad = cloudparams[3]
            sig = cloudparams[4]
            pdiff = np.empty(nlayers,dtype='f')
 
            for j in range(0, ncloud):
                pdiff = abs(np.log(press) - np.log(p1))
                l1 = np.argmin(pdiff)
                # Whichever layer-mean-P the cloud top is closest
                # to is the +1 layer of the cloud
                # same for the -1 layer of the base
                pdiff = abs(np.log(press) - np.log(p2))
                l2 = np.argmin(pdiff)

                if (cloudnum != 99):
                    cloudprof[i,:,j] = -150.00
                    if (l1 == l2):
                        pl1, pl2 = atlev(l1,press)
                        cloudprof[i,l1,j] = np.log10((10.**ndens) * 
                                                     ((np.log10(p2) - np.log10(p1))
                                                      /
                                                      (np.log10(pl2) - np.log10(pl1))))
                    else:
                        cloudprof[i,l1+1:l2,j] = ndens
                        # now account for partial fill of bottom and top layers
                        # need levels for this
                        pl1, pl2 = atlev(l1,press)
                        # This is done in logP for now
                        cloudprof[i,l1,j] = np.log10((10.**ndens) * \
                                                 ((np.log10(pl2) - np.log10(p1))
                                                  /
                                                  (np.log10(pl2) - np.log10(pl1))))
                        # same for bottom
                        pl1, pl2 = atlev(l2,press)
                        cloudprof[i,l2,j] = np.log10(10.**ndens * \
                                                ((np.log10(p2) - np.log10(pl1))
                                                 /
                                                 (np.log10(pl2) - np.log10(pl1))))
                else:
                    # This is a slab cloud
                    # dtau/dP propto P
                    tau = ndens
                    if (l1 == l2):
                        cloudprof[i,l1,j] = tau
                    else:
                        const = tau / (p2**2 - p1**2)
                        # partial top fill
                        pl1, pl2 = atlev(l1,press)
                        cloudprof[i,l1,j] = const * (pl2**2 - p1**2)
                        # partial bottom fill
                        pl1, pl2 = atlev(l2,press)
                        cloudprof[i,l2,j] = const *  (p2**2 - pl1**2) 
                        for k in range (l1+1,l2):
                            l1,l2 = atlev(k,press)
                            cloudprof[j,k,j] = const * (l2**2 - l1**2)
 
                cloudrad[i,:,j] = rad
                cloudsig[i,:,j] = sig        

        if (do_clouds[i]==1 and cloudtype[i] == 2):

            # 5 entries for cloudparams are:
            # 0) reference density 
            # 1) top pressure
            # 2) scale height (in dex)
            # 3) rg
            # 4) rsig
            
            ndens= cloudparams[0]
            p0 = 10.**cloudparams[1]
            scale = ((p0 * 10.**cloudparams[2]) - p0)  / 10.**cloudparams[2]
            rad = cloudparams[3]
            sig = cloudparams[4]
            
            pdiff = np.empty(nlayers,dtype='f')
            
            for j in range(0, ncloud):
                
                if (cloudnum != 99):
                    cloudprof[i,:,j] = -150.00
                    pdiff = abs(np.log(press) - np.log(p0))
                    l0 = np.argmin(pdiff)
                    if (l0 <= nlayers-2):
                        cloudprof[i,l0+1:,j] = ndens
                    # now for top layer of cloud
                    # need levels for this
                    pl1, pl2 = atlev(l0,press)
                    # want geometric weighted mean of level above p0 and ndens
                    dl1 = (10.**ndens) * np.exp((pl1 - p0)/scale)
                    cloudprof[i,l0,j] = np.log10(np.sqrt(dl1 *  10.**ndens))
                    # now the rest of the layers
                    for k in range(0, l0):
                        cloudprof[i,k,j] = np.log10((10.**ndens) * np.exp((press[k] - p0)/scale))
                else:
                    #  cloud 99 case!!
                    # Here P0 is the pressure where tau = 1 for the cloud
                    # so dtau / dP = const * exp((P-P0) / scale)
                    # See notes for derivation of constant and integral
                    const = 1. / (1 - np.exp(-p0 / scale))
                    for k in range (0,nlayers):
                        pl1, pl2 = atlev(k,press)
                        # now get dtau for each layer, where tau = 1 at P0
                        term1 = (pl2 - p0) / scale
                        term2 = (pl1 - p0) / scale
                        if (term1 > 10 or term2 > 10):
                            cloudprof[i,k,j] = 100.00
                        else:
                            cloudprof[i,k,j] = const * (np.exp(term1) -
                                                        np.exp(term2))
                                                     
                cloudrad[i,:,j] = rad
                cloudsig[i,:,j] = sig       


                #print "dl1 = ",dl1
                    
        if (cloudtype[i] != 1 and cloudtype != 2):
            print "cloud layout not recognised. stopping" 

    
    return cloudprof,cloudrad,cloudsig

def atlev(l0,press):
    nlayers = press.size
    if (l0 <= nlayers-2):
        pl1 = np.exp(((1.5)*np.log(press[l0])) - ((0.5)*np.log(press[l0+1])))
        pl2 = np.exp((0.5)*(np.log(press[l0] * press[l0+1])))
    else:
        pl1 = np.exp((0.5 * np.log(press[l0-1] * press[l0])))
        pl2 = press[l0]**2 / pl1

    return pl1, pl2
