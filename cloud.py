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
            # 1) top pressure
            # 2) base pressure
            # 3) rg
            # 4) rsig
            


            ndens= cloudparams[0]
            p1 = cloudparams[1]
            p2 = cloudparams[2]
            rad = cloudparams[3]
            sig = cloudparams[4]
            pdiff = np.empty(nlayers,dtype='f')
 
            for j in range(0, ncloud):
                pdiff = abs(np.log(press) - np.log(p1))
                l1 = np.argmin(pdiff)
                # Whichever layer mean P the cloud top is closest
                # to is the +1 layer of the cloud
                # same for the -1 layer of the base
                pdiff = abs(np.log(press) - np.log(p2))
                l2 = np.argmin(pdiff)
                cloudprof[i,l1+1:l2-1,j] = ndens
                # now account for partial fill of bottom and top layers
                # need levels for this
                pl1 = np.exp(((1.5)*np.log(press[l1])) - ((0.5)*np.log(press[l1+1])))
                pl2 = np.exp((0.5)*(np.log(press[l1] * press[l1+1])))
                # This is done in logP for now
                if (cloudnum != 99):
                    cloudprof[i,l1,j] = np.log10(10.**ndens * \
                                                 ((np.log10(pl2) - np.log10(p1)) /
                                                  (np.log10(pl2) - np.log10(pl1))))
                else:
                     cloudprof[i,l1,j] = ndens * ((np.log10(pl2) - np.log10(p1)) /
                                                  (np.log10(pl2) - np.log10(pl1)))
                    
                # same for bottom
                pl1 = np.exp(((1.5)*np.log(press[l2])) - ((0.5)*np.log(press[l2+1])))
                pl2 = np.exp((0.5)*(np.log(press[l2] * press[l2+1])))
                # This is done in logP for now
                if (cloudnum != 99):
                    cloudprof[i,l2,j] = np.log10(10.**ndens * \
                                                 ((np.log10(p2) - np.log10(pl1)) /
                                                  (np.log10(pl2) - np.log10(pl1))))
                else:
                    cloudprof[i,l2,j] = ndens * ((np.log10(p2) - np.log10(pl1)) /
                                                  (np.log10(pl2) - np.log10(pl1)))
 
                cloudrad[i,:,j] = rad
                cloudsig[i,:,j] = sig        

        if (do_clouds[i]==1 and cloudtype[i] == 2):

            # 5 entries for cloudparams are:
            # 0) reference density 
            # 1) top pressure
            # 2) scale height
            # 3) rg
            # 4) rsig
            
            ndens= cloudparams[0]
            p0 =cloudparams[1]
            scale = cloudparams[2]
            rad = cloudparams[3]
            sig = cloudparams[4]
            
            pdiff = np.empty(nlayers,dtype='f')
            
            for j in range(0, ncloud):
                
                pdiff = abs(np.log(press) - np.log(p0))
                l0 = np.argmin(pdiff)
                cloudprof[i,l0+1:,j] = ndens
                # now for top layer of cloud
                # need levels for this
                pl1 = np.exp(((1.5)*np.log(press[l0])) - ((0.5)*np.log(press[l0+1])))
                pl2 = np.exp((0.5)*(np.log(press[l0] * press[l0+1])))
                print "l0 =",l0
                print "p0 =",p0
                print "pl1 =", pl1
                print "pl2 =",pl2
                # This is done in logP for now
                # want geometric weighted mean of top layer, and level above
                if (cloudnum != 99):
                    dl1 = (10.**ndens) * np.exp((pl1 - p0)/scale)
                    cloudprof[i,l0,j] = np.log10(np.sqrt((dl1 * ((p0 - pl1)/(pl2-pl1))) *
                                                ((10.**ndens)*((pl2 - p0)/(pl2-pl1)))))
                    # now the rest of the layers
                    for k in range(0, l0):
                        cloudprof[i,k,j] = np.log10((10.**ndens) * np.exp((press[k] - p0)/scale))
                else:
                    dl1 = ndens * np.exp((pl1 - p0)/scale)
                    cloudprof[i,l0,j] = np.sqrt((dl1 * ((p0 - pl1)/(pl2-pl1))) *
                                                ((ndens)*((pl2 - p0)/(pl2-pl1))))
                    # now the rest of the layers
                    for k in range(0, l0):
                        cloudprof[i,k,j] = ndens * np.exp((press[k] - p0)/scale)

                cloudrad[i,:,j] = rad
                cloudsig[i,:,j] = sig       


                print "dl1 = ",dl1
                    
        if (cloudtype[i] != 1 and cloudtype != 2):
            print "cloud layout not recognised. stopping" 

    
    return cloudprof,cloudrad,cloudsig
