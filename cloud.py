#!/usr/bin/env python

""" Module of processes to interpret cloud parameters from Brewster in testkit"""
from __future__ import print_function

import numpy as np
import scipy as sp
from scipy import interpolate



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
    # 0) dtau at 1um
    # 1) top layer id (or pressure)
    # 2) base ID (these are both in 64 layers)
    # 3) rg
    # 4) rsig
    # in the case of a simple mixto cloud (i.e. cloudnum = 99 or 89) we have:
    # 0) dtau (at 1um for non-grey) 
    # 1) top layer ID
    # 2) bottom later ID
    # 3) rg = albedo
    # 4) rsig = power for tau power law



def atlas(do_clouds,cloudnum,cloudtype,cloudparams,press):

    # Cloud types
    # 1:  slab cloud
    # 2: deep thick cloud , we only see the top
    # 3: slab with fixed thickness log dP = 0.005 (~1% height)
    # 4: deep thick cloud with fixed height log dP = 0.005
    # In both cases the cloud properties are density, rg, rsig for real clouds
    # and dtau, w0, and power law for cloudnum = 99 or 89
    nlayers = press.size
    npatch = do_clouds.size
    ncloud = 1
    if (cloudparams.size > 5):
        ncloud = cloudparams.shape[2]
    
    cloudrad = np.zeros((npatch,nlayers,ncloud),dtype='d')
    cloudsig = np.zeros_like(cloudrad)
    cloudprof = np.zeros_like(cloudrad)

    
    for i in range(0, npatch):
        if (do_clouds[i] != 0):
            for j in range(0,ncloud):

                if (cloudtype[i,j] == 1 or cloudtype[i,j] == 3):
                    # 5 entries for cloudparams are:
                    # 0) total tau for cloud at 1 micron
                    # 1) log top pressure
                    # 2) pressure thickness in dex
                    # 3) rg
                    # 4) rsig
            


                    tau = cloudparams[0,i,j]
                    p1 = 10.**cloudparams[1,i,j]
                    if (cloudtype[i,j] == 1):
                        dP = cloudparams[2,i,j]
                    else:
                        dP = 0.005
                    p2 = p1 * 10.**dP
                    rad = cloudparams[3,i,j]
                    sig = cloudparams[4,i,j]
                    pdiff = np.empty(nlayers,dtype='f')
 
                    pdiff = abs(np.log(press) - np.log(p1))
                    l1 = np.argmin(pdiff)
                    # Whichever layer-mean-P the cloud top is closest
                    # to is the +1 layer of the cloud
                    # same for the -1 layer of the base
                    pdiff = abs(np.log(press) - np.log(p2))
                    l2 = np.argmin(pdiff)

                    # This is a slab cloud
                    # dtau/dP propto P
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
                            cloudprof[i,k,j] = const * (l2**2 - l1**2)

                    # We're sampling particle radius in log space        
                    if (cloudnum[i,j] < 50.):
                        cloudrad[i,:,j] = 10.**rad
                    else:
                        cloudrad[i,:,j] = rad
                    cloudsig[i,:,j] = sig        

                if (cloudtype[i,j] == 2 or cloudtype[i,j] == 4):

                    # 5 entries for cloudparams are:
                    # 0) empty
                    # 1) top pressure
                    # 2) scale height (in dex)
                    # 3) rg
                    # 4) rsig
            
                    p0 = 10.**cloudparams[1,i,j]
                    if (cloudtype[i,j] == 2):
                        dP = cloudparams[2,i,j]
                    else:
                        dP = 0.005
                    scale = ((p0 * 10.**dP) - p0)  / 10.**dP
                    rad = cloudparams[3,i,j]
                    sig = cloudparams[4,i,j]
            
                    pdiff = np.empty(nlayers,dtype='f')
            
                
                    # In cloud 99/89 case rsig is power law for tau~lambda^alpha 
                    # Here P0 is the pressure where tau= 1 for the cloud
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

                    # We're sampling particle radius in log space        
                    if (cloudnum[i,j] < 50.):
                        cloudrad[i,:,j] = 10.**rad
                    else:
                        cloudrad[i,:,j] = rad
                    cloudsig[i,:,j] = sig       

                if (cloudtype[i,j] == 0):
                    cloudprof[i,:,j] = 0.0
                    cloudrad[i,:,j] = 0.0
                    cloudsig[i,:,j] = 0.0

                if (cloudtype[i,j] >  4):
                    print ("cloud layout not recognised. stopping") 

    
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


def unpack_default(theta,pc,cloudtype,cloudnum,do_clouds):
    nc =0
    if (cloudtype.size > cloudtype.shape[1]):
        nclouds = cloudtype.shape[1]
    else:
        nclouds = cloudtype.size


 
    npatches = do_clouds.size

    cloudparams = np.ones([5,npatches,nclouds],dtype='d')
    cloudparams[0,:,:] = 0.
    cloudparams[1,:,:] = 0.0
    cloudparams[2,:,:] = 0.1
    cloudparams[3,:,:] = 0.0
    cloudparams[4,:] = 0.5

    for i in range (0,npatches):
        if (do_clouds[i] != 0):
            for j in range (0, nclouds):
                if ((cloudtype[i,j] == 2) and (cloudnum[i,j] == 99)):
                    cloudparams[1:4,i,j] = theta[pc+nc:pc+3+nc]
                    cloudparams[4,i,j] = 0.0
                    nc = nc + 3
                elif ((cloudtype[i,j] == 1) and (cloudnum[i,j] == 99)):
                    cloudparams[0:4,i,j] = theta[pc+nc:pc+4+nc]
                    cloudparams[4,i,j] = 0.0
                    nc = nc + 4
                elif ((cloudtype[i,j] == 2) and (cloudnum[i,j] < 90)):
                    cloudparams[1:5,i,j] = theta[pc+nc:pc+4+nc]
                    nc = nc +4
                elif ((cloudtype[i,j] == 3) and (cloudnum[i,j] == 99)):
                    cloudparams[0:2,i,j] = theta[pc+nc:pc+nc+2]
                    cloudparams[3,i,j] =  theta[pc+nc+2]
                    nc = nc +3
                elif ((cloudtype[i,j] == 3) and (cloudnum[i,j] < 90)):
                    cloudparams[0:2,i,j] = theta[pc+nc:pc+nc+2]
                    cloudparams[3:5,i,j] =  theta[pc+nc+2:pc+nc+4]
                    nc = nc + 4
                elif ((cloudtype[i,j] == 4) and (cloudnum[i,j] == 99)):
                    cloudparams[1,i,j] = theta[pc+nc]
                    cloudparams[3,i,j] = theta[pc+nc+1]
                    nc = nc +2
                elif ((cloudtype[i,j] == 4) and (cloudnum[i,j] < 90)):
                    cloudparams[1,i,j] = theta[pc+nc]
                    cloudparams[3:5,i,j] = theta[pc+nc+1:pc+nc+3]
                    nc = nc +3
                elif (cloudtype[i,j] == 0):
                    cloudparams[:,i,j] = 0.0
                else:
                    cloudparams[:,i,j] = theta[pc+nc:pc+5+nc]
                    nc = nc + 5

    return(cloudparams,nc)


def unpack_patchy(theta,pc,cloudtype,cloudnum,do_clouds):
    # This unpacks a patchy cloud
    # This is achieved by equated patch 2 to patch 1,
    # but all but some clouds in patch 2 are set to cloudtype = 0
    nc =0 
    if (cloudtype.size > cloudtype.shape[1]):
        nclouds = cloudtype.shape[1]
    else:
        nclouds = cloudtype.size
        

    npatches = 2

    cloudparams = np.ones([5,npatches,nclouds],dtype='d')
    cloudparams[0,:,:] = 0.
    cloudparams[1,:,:] = 0.0
    cloudparams[2,:,:] = 0.1
    cloudparams[3,:,:] = 0.0
    cloudparams[4,:] = 0.0

    # First patch
    if (do_clouds[0] != 0):
        for j in range (0, nclouds):
            if ((cloudtype[0,j] == 2) and (cloudnum[0,j] == 99)):
                cloudparams[1:4,0,j] = theta[pc+nc:pc+3+nc]
                cloudparams[4,0,j] = 0.0
                nc = nc + 3
            elif ((cloudtype[0,j] == 1) and (cloudnum[0,j] == 99)):
                cloudparams[0:4,0,j] = theta[pc+nc:pc+4+nc]
                cloudparams[4,0,j] = 0.0
                nc = nc + 4
            elif ((cloudtype[0,j] == 2) and (cloudnum[0,j] < 90)):
                cloudparams[1:5,0,j] = theta[pc+nc:pc+4+nc]
                nc = nc +4
            elif ((cloudtype[0,j] == 3) and (cloudnum[0,j] == 99)):
                cloudparams[0:2,0,j] = theta[pc+nc:pc+nc+2]
                cloudparams[3,0,j] =  theta[pc+nc+2]
                nc = nc +3
            elif ((cloudtype[0,j] == 3) and (cloudnum[0,j] < 90)):
                cloudparams[0:2,0,j] = theta[pc+nc:pc+nc+2]
                cloudparams[3:5,0,j] =  theta[pc+nc+2:pc+nc+4]
                nc = nc + 4
            elif ((cloudtype[0,j] == 4) and (cloudnum[0,j] == 99)):
                cloudparams[1,0,j] = theta[pc+nc]
                cloudparams[3,0,j] = theta[pc+nc+1]
                nc = nc +2
            elif ((cloudtype[0,j] == 4) and (cloudnum[0,j] < 90)):
                cloudparams[1,0,j] = theta[pc+nc]
                cloudparams[3:5,0,j] = theta[pc+nc+1:pc+nc+3]
                nc = nc +3                    
            else:
                cloudparams[:,0,j] = theta[pc+nc:pc+5+nc]
                nc = nc + 5
                    
    # 2nd patch
    if (do_clouds[1] != 0):
        cloudparams[:,1,:] = cloudparams[:,0,:]
        
    return(cloudparams,nc)
