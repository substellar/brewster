#!/usr/bin/env python

""" Module of bits to plug into Brewster """
from __future__ import print_function
import math
import time
import gc
import numpy as np
import scipy as sp
import forwardmodel
import cloud
import TPmod
import settings
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import interp1d

from astropy.convolution import convolve, convolve_fft
from astropy.convolution import Gaussian1DKernel
from bensconv import spex_non_uniform
from bensconv import conv_uniform_R
from bensconv import conv_uniform_FWHM

__author__ = "Ben Burningham"
__copyright__ = "Copyright 2015 - Ben Burningham"
__credits__ = ["Ben Burningham","The EMCEE DOCS"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ben Burningham"
__email__ = "burninghamster@gmail.com"
__status__ = "Development"


    
def lnprob(theta):


    # now check against the priors, if not beyond them, run the likelihood
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    # run the likelihood
    lnlike_value = lnlike(theta)

    lnprb = lp+lnlike_value
    if np.isnan(lnprb):
        lnprb = -np.inf
    return lnprb


def lnprior(theta):

    gases_myP,chemeq,dist,cloudtype, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge,prof,do_bff,bff_raw,ceTgrid,metscale,coscale = settings.runargs

    # set up the priors here
    if (chemeq != 0):
        invmr = np.array([-3.,-3.])
        mh = theta[0]
        co = theta[1]
        ng = 2
    else:
        if (gasnum[gasnum.size-1] == 21):
            ng = gasnum.size - 1
            mh = 0.0
            co = 1.0
            invmr = theta[0:ng]
        elif (gasnum[gasnum.size-1] == 23):
            mh = 0.0
            co = 1.0
            ng = gasnum.size - 2
            invmr = theta[0:ng]
        else:
            ng = gasnum.size
            invmr = theta[0:ng]
            mh = 0.0
            co = 1.0
            
    logg = theta[ng]
    if (fwhm < 0.0):
        if (fwhm == -1 or fwhm == -3):
            s1  = np.where(obspec[0,:] < 2.5)
            s2  = np.where(np.logical_and(obspec[0,:] > 2.5,obspec[0,:] < 5.0))
            s3 =  np.where(obspec[0,:] > 5.)
            r2d2 = theta[ng+1]
            scale1 = theta[ng+2]
            scale2 = theta[ng+3]
            dlam = theta[ng+4]
            if (do_fudge == 1):
                logf1 = theta[ng+5]
                logf2 = theta[ng+6]
                logf3 = theta[ng+7]
                logf = np.log10(0.1*(max(obspec[2,10::3]))**2) 
                pc = ng+8
            else:
                # This is a place holder value so the code doesn't break
                logf = np.log10(0.1*(max(obspec[2,:]))**2)
                logf1 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf2 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf3 = np.log10(0.1*(max(obspec[2,:]))**2)
                pc = ng + 5
        elif (fwhm == -2):
            s1  = np.where(obspec[0,:] < 2.5)
            s2 = s1
            s3 =  np.where(obspec[0,:] > 5.)
            r2d2 = theta[ng+1]
            scale1 = 1.0 # dummy value
            scale2 = theta[ng+2]
            dlam = theta[ng+3]
            if (do_fudge == 1):
                logf1 = theta[ng+4]
                logf2 = np.log10(0.1*(max(obspec[2,10::3]))**2) # dummy
                logf3 = theta[ng+5]
                logf = np.log10(0.1*(max(obspec[2,10::3]))**2) 
                pc = ng+6
            else:
                # This is a place holder value so the code doesn't break
                logf = np.log10(0.1*(max(obspec[2,:]))**2)
                logf1 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf2 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf3 = np.log10(0.1*(max(obspec[2,:]))**2)
                pc = ng + 4

    else:
        # this just copes with normal, single instrument data
        s1 = np.where(obspec[0,:] > 0.0)
        s2 = s1
        s3 = s1
        
        r2d2 = theta[ng+1]
        dlam = theta[ng+2]
        if (do_fudge == 1):
            logf = theta[ng+3]
            logf1 = np.log10(0.1*(max(obspec[2,:]))**2)
            logf2 = np.log10(0.1*(max(obspec[2,:]))**2)
            logf3 = np.log10(0.1*(max(obspec[2,:]))**2)
            scale1 = 1.0
            scale2 = 1.0
            pc = ng + 4
        else:
            # This is a place holder value so the code doesn't break
            logf = np.log10(0.1*(max(obspec[2,10::3]))**2)
            logf1 = np.log10(0.1*(max(obspec[2,:]))**2)
            logf2 = np.log10(0.1*(max(obspec[2,:]))**2)
            logf3 = np.log10(0.1*(max(obspec[2,:]))**2)
            scale1 = 1.0
            scale2 = 1.0
            pc = ng + 3


    npatches = do_clouds.size
    if (npatches > 1):
        prat = theta[pc]
        pcover = np.array([prat,(1.-prat)])
        pc = pc + 1
    else:
        pcover = np.array([0.5,0.5])
        
    # use correct unpack method depending on situation
    
    if ((npatches > 1) and np.all(do_clouds != 0)):
        cloudparams, nc = cloud.unpack_patchy(theta,pc,cloudtype,cloudnum,do_clouds)
    else:
        cloudparams, nc = cloud.unpack_default(theta,pc,cloudtype,cloudnum,do_clouds)


    if (cloudtype.size > cloudtype.shape[0]):
        nclouds = cloudtype.shape[1]
    else:
        nclouds = cloudtype.size
    
    cloud_tau0 = np.empty([npatches,nclouds])
    cloud_top = np.empty_like(cloud_tau0)
    cloud_bot = np.empty_like(cloud_tau0)
    cloud_height  = np.empty_like(cloud_tau0)
    w0 = np.empty_like(cloud_tau0)
    taupow =  np.empty_like(cloud_tau0)
    loga = np.empty_like(cloud_tau0)
    b = np.empty_like(cloud_tau0)

    if (np.abs(sum(do_clouds)) >= 1):    
        for i in range(0,npatches):
            if (do_clouds[i] != 0):
                for j in range (0, nclouds):
                    if (cloudnum[i,j] == 99):
                        if (cloudtype[i,j] == 1):
                            cloud_tau0[i,j] = cloudparams[0,i,j]
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = cloudparams[2,i,j]
                            cloud_bot[i,j] = cloud_top[i,j] + cloud_height[i,j]
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = 0.0
                            loga[i,j] = 0.0
                            b[i,j] = 0.5
                        elif (cloudtype[i,j] == 2):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[-1])
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = cloudparams[2,i,j]
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = 0.0
                            loga[i,j] = 0.0
                            b[i,j] = 0.5
                        elif (cloudtype[i,j] == 3):
                            cloud_tau0[i,j] = cloudparams[0,i,j]
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = 0.005
                            cloud_bot[i,j] = cloud_top[i,j] + cloud_height[i,j]
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = 0.0
                            loga[i,j] = 0.0
                            b[i,j] = 0.5
                        elif (cloudtype[i,j] == 4):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[-1])
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = 0.005
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = 0.0
                            loga[i,j] = 0.0
                            b[i,j] = 0.5
                        elif (cloudtype[i,j] == 0):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[press.size-1])
                            cloud_top[i,j] = np.log10(press[0])
                            cloud_height[i,j] = 0.1
                            w0[i,j] = +0.5
                            taupow[i,j] = 0.0
                            loga[i,j] =  0.0
                            b[i,j] = 0.5
                    elif (cloudnum[i,j] == 89):
                        if (cloudtype[i,j] == 1):
                            cloud_tau0[i,j] = cloudparams[0,i,j]
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = cloudparams[2,i,j]
                            cloud_bot[i,j] = cloud_top[i,j] + cloud_height[i,j]
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = cloudparams[4,i,j]
                            loga[i,j] = 0.0
                            b[i,j] = 0.5
                        elif (cloudtype[i,j] == 2):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[-1])
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = cloudparams[2,i,j]
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = cloudparams[4,i,j]
                            loga[i,j] = 0.0
                            b[i,j] = 0.5
                        elif (cloudtype[i,j] == 3):
                            cloud_tau0[i,j] = cloudparams[0,i,j]
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = 0.005
                            cloud_bot[i,j] = cloud_top[i,j] + cloud_height[i,j]
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = cloudparams[4,i,j]
                            loga[i,j] = 0.0
                            b[i,j] = 0.5
                        elif (cloudtype[i,j] == 4):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[-1])
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = 0.005
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = cloudparams[4,i,j]
                            loga[i,j] = 0.0
                            b[i,j] = 0.5
                        elif (cloudtype[i,j] == 0):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[press.size-1])
                            cloud_top[i,j] = np.log10(press[0])
                            cloud_height[i,j] = 0.1
                            w0[i,j] = +0.5
                            taupow[i,j] = 0.0
                            loga[i,j] =  0.0
                            b[i,j] = 0.5
                    else:
                        if (cloudtype[i,j] == 1):
                            cloud_tau0[i,j] =  cloudparams[0,i,j]
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = cloudparams[2,i,j]
                            cloud_bot[i,j] = cloud_top[i,j] + cloud_height[i,j]
                            w0[i,j] = 0.5
                            taupow[i,j] = 0.0
                            loga[i,j] = cloudparams[3,i,j]
                            b[i,j] = cloudparams[4,i,j]
                        elif (cloudtype[i,j] == 2):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[press.size-1])
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = cloudparams[2,i,j]
                            w0[i,j] = +0.5
                            taupow[i,j] =0.0
                            loga[i,j] =  cloudparams[3,i,j]
                            b[i,j] =  cloudparams[4,i,j]
                        elif (cloudtype[i,j] == 3):
                            cloud_tau0[i,j] = cloudparams[0,i,j]
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = 0.005
                            cloud_bot[i,j] = cloud_top[i,j] + cloud_height[i,j]
                            w0[i,j] = 0.5
                            taupow[i,j] = 0.0
                            loga[i,j] = cloudparams[3,i,j]
                            b[i,j] = cloudparams[4,i,j]
                        elif (cloudtype[i,j] == 4):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[-1])
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = 0.005
                            w0[i,j] = +0.5
                            taupow[i,j] =0.0
                            loga[i,j] =  cloudparams[3,i,j]
                            b[i,j] =  cloudparams[4,i,j]
                        elif (cloudtype[i,j] == 0):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[-1])
                            cloud_top[i,j] = np.log10(press[0])
                            cloud_height[i,j] = 0.1
                            w0[i,j] = +0.5
                            taupow[i,j] = 0.0
                            loga[i,j] =  0.0
                            b[i,j] = 0.5
    else:
        cloud_tau0[:,:] = 1.0
        cloud_bot[:,:] = np.log10(press[-1])
        cloud_top[:,:] = np.log10(press[0])
        cloud_height[:,:] = 0.1
        w0[:,:] = +0.5
        taupow[:,:] = 0.0
        loga[:,:] =  0.0
        b[:,:] = 0.5
                
    junkP = np.ones([13])
    if (proftype == 1):
        gam = theta[pc+nc]
        T = theta[pc+nc+1:]
        diff=np.roll(T,-1)-2.*T+np.roll(T,1)
        pp=len(T)
    
        #for mass prior
        D = 3.086e+16 * dist
        R = -1.0
        if (r2d2 > 0.):
            R = np.sqrt(r2d2) * D 
        g = (10.**logg)/100.
        M = (R**2 * g/(6.67E-11))/1.898E27
        Rj = R / 69911.e3 
        #         and  and (-5. < logbeta < 0))
        if (all(invmr[0:ng] > -12.0) and all(invmr[0:ng] < 0.0) and (np.sum(10.**(invmr[0:ng])) < 1.0)
            and all(pcover > 0.) and (np.sum(pcover) == 1.0)
            and  metscale[0] <= mh <= metscale[-1]
            and  coscale[0] <= co <= coscale[-1]
            and  0.0 < logg < 6.0 
            and 1.0 < M < 80
            and  0. < r2d2 < 1.
            and 0.1 < scale1 < 10.0
            and 0.1 < scale2 < 10.0
            and  0.5 < Rj < 2.0
            and -0.01 < dlam < 0.01 
            and (min(T) > 1.0) and (max(T) < 6000.) 
            and (gam > 0.)
            and ((0.01*np.min(obspec[2,:]**2)) < 10.**logf
                 < (100.*np.max(obspec[2,:]**2)))
            and ((0.01*np.min(obspec[2,s1]**2)) < 10.**logf1
                 < (100.*np.max(obspec[2,s1]**2)))
            and ((0.01*np.min(obspec[2,s2]**2)) < 10.**logf2
                 < (100.*np.max(obspec[2,s2]**2)))
            and ((0.01*np.min(obspec[2,s3]**2)) < 10.**logf3
                 < (100.*np.max(obspec[2,s3]**2)))
            and (np.all(cloud_tau0 >= 0.0))
            and (np.all(cloud_tau0 <= 100.0))
            and np.all(cloud_top < cloud_bot)
            and np.all(cloud_bot <= np.log10(press[-1]))
            and np.all(np.log10(press[0]) <= cloud_top)
            and np.all(cloud_top < cloud_bot)
            and np.all(0. < cloud_height)
            and np.all(cloud_height < 7.0)
            and np.all(0.0 < w0)
            and np.all(w0 <= 1.0)
            and np.all(-10.0 < taupow)
            and np.all(taupow < +10.0)
            and np.all( -3.0 < loga)
            and np.all (loga < 3.0)
            and np.all(b < 1.0)
            and np.all(b > 0.0)):
                 
            logbeta = -5.0
            beta=10.**logbeta
            alpha=1.0
            x=gam
            invgamma=((beta**alpha)/math.gamma(alpha)) * (x**(-alpha-1)) * np.exp(-beta/x)
            prprob = (-0.5/gam)*np.sum(diff[1:-1]**2) - 0.5*pp*np.log(gam) + np.log(invgamma)
            return prprob 
        
        return -np.inf

    elif (proftype == 2):
        a1 = theta[pc+nc]
        a2 = theta[pc+1+nc]
        P1 = theta[pc+2+nc]
        P3 = theta[pc+3+nc]
        T3 = theta[pc+4+nc]
        T = np.empty([press.size])
        T[:] = -100.

        if (0. < a1 < 1. and 0. < a2 < 1.0
            and T3 > 0.0 and P3 >= P1 and P1 >= np.log10(press[0])
            and P3 < 5.):
            T = TPmod.set_prof(proftype,junkP,press,theta[pc+nc:])
        #for mass prior
        D = 3.086e+16 * dist
        R = -1.0
        if (r2d2 > 0.):
            R = np.sqrt(r2d2) * D 
        g = (10.**logg)/100.
        M = (R**2 * g/(6.67E-11))/1.898E27
        Rj = R / 69911.e3 
        #         and  and (-5. < logbeta < 0))
        if (all(invmr[0:ng] > -12.0) and all(invmr[0:ng] < 0.0) and (np.sum(10.**(invmr[0:ng])) < 1.0)
            and  all(pcover > 0.) and (np.sum(pcover) == 1.0)
            and  metscale[0] <= mh <= metscale[-1]
            and  coscale[0] <= co <= coscale[-1]
            and  0.0 < logg < 6.0 
            and 1.0 < M < 80. 
            and  0. < r2d2 < 1.
            and 0.1 < scale1 < 10.0
            and 0.1 < scale2 < 10.0
            and  0.5 < Rj < 2.0
            and -0.01 < dlam < 0.01 
            and ((0.01*np.min(obspec[2,:]**2)) < 10.**logf
                 < (100.*np.max(obspec[2,:]**2)))
            and ((0.01*np.min(obspec[2,s1]**2)) < 10.**logf1
                 < (100.*np.max(obspec[2,s1]**2)))
            and ((0.01*np.min(obspec[2,s2]**2)) < 10.**logf2
                 < (100.*np.max(obspec[2,s2]**2)))
            and ((0.01*np.min(obspec[2,s3]**2)) < 10.**logf3
                 < (100.*np.max(obspec[2,s3]**2)))
            and (np.all(cloud_tau0 >= 0.0))
            and (np.all(cloud_tau0 <= 100.0))
            and np.all(cloud_top < cloud_bot)
            and np.all(cloud_bot <= np.log10(press[-1]))
            and np.all(np.log10(press[0]) <= cloud_top)
            and np.all(cloud_top < cloud_bot)
            and np.all(0. < cloud_height)
            and np.all(cloud_height < 7.0)
            and np.all(0.0 < w0)
            and np.all(w0 <= 1.0)
            and np.all(-10.0 < taupow)
            and np.all(taupow < +10.0)
            and np.all( -3.0 < loga)
            and np.all (loga < 3.0)
            and np.all(b < 1.0)
            and np.all(b > 0.0)
            and  (min(T) > 1.0)
            and (max(T) < 6000.)):
            return 0.0
        return -np.inf
        
    elif (proftype == 3):
        a1 = theta[pc+nc]
        a2 = theta[pc+1+nc]
        P1 = theta[pc+2+nc]
        P2 = theta[pc+3+nc]
        P3 = theta[pc+4+nc]
        T3 = theta[pc+5+nc]

        T = np.empty([press.size])
        T[:] = -100.
        if  (0. < a1 < 1. and 0. < a2 < 1.0
            and T3 > 0.0 and P3 >= P2 and P2 >= P1 and P1 >= np.log10(press[0])
             and P3 < 5.):
            T = TPmod.set_prof(proftype,junkP,press,theta[pc+nc:])

            
        #for mass prior
        D = 3.086e+16 * dist
        R = -1.0
        if (r2d2 > 0.):
            R = np.sqrt(r2d2) * D 
        g = (10.**logg)/100.
        M = (R**2 * g/(6.67E-11))/1.898E27
        Rj = R / 69911.e3 
        #         and  and (-5. < logbeta < 0))
        if (all(invmr[0:ng] > -12.0) and all(invmr[0:ng] < 0.0) and (np.sum(10.**(invmr[0:ng])) < 1.0) 
            and  all(pcover > 0.) and (np.sum(pcover) == 1.0)
            and  metscale[0] <=  mh <= metscale[-1]
            and  coscale[0] <= co <= coscale[-1]
            and  0.0 < logg < 6.0 
            and 1.0 < M < 80. 
            and  0. < r2d2 < 1.
            and 0.1 < scale1 < 10.0
            and 0.1 < scale2 < 10.0
            and  0.5 < Rj < 2.0
            and -0.01 < dlam < 0.01 
            and ((0.01*np.min(obspec[2,:]**2)) < 10.**logf
                 < (100.*np.max(obspec[2,:]**2)))
            and ((0.01*np.min(obspec[2,s1]**2)) < 10.**logf1
                 < (100.*np.max(obspec[2,s1]**2)))
            and ((0.01*np.min(obspec[2,s2]**2)) < 10.**logf2
                 < (100.*np.max(obspec[2,s2]**2)))
            and ((0.01*np.min(obspec[2,s3]**2)) < 10.**logf3
                 < (100.*np.max(obspec[2,s3]**2)))
            and (np.all(cloud_tau0 >= 0.0))
            and (np.all(cloud_tau0 <= 100.0))
            and np.all(cloud_top < cloud_bot)
            and np.all(cloud_bot <= np.log10(press[-1]))
            and np.all(np.log10(press[0]) <= cloud_top)
            and np.all(cloud_top < cloud_bot)
            and np.all(0. < cloud_height)
            and np.all(cloud_height < 7.0)
            and np.all(0.0 < w0)
            and np.all(w0 <= 1.0)
            and np.all(-10.0 < taupow)
            and np.all(taupow < +10.0)
            and np.all( -3.0 < loga)
            and np.all (loga < 3.0)
            and np.all(b < 1.0)
            and np.all(b > 0.0)
            and  (min(T) > 1.0) and (max(T) < 6000.)):               
            return 0.0
        return -np.inf
    
    elif (proftype == 9):
        #for mass prior
        D = 3.086e+16 * dist
        R = -1.0
        if (r2d2 > 0.):
            R = np.sqrt(r2d2) * D 
        g = (10.**logg)/100.
        M = (R**2 * g/(6.67E-11))/1.898E27
        Rj = R / 69911.e3 
        #         and  and (-5. < logbeta < 0))
        print("Rj = ", Rj)
        print("M = ", M)
        print("[M/H] = ",mh)
        print("[C/O] = ", co)
        print("logg = ", logg)
        print("R2D2 = ", r2d2)
        print("scale1 = ",scale1)
        print("scale2 = ",scale1)
        print("Rj = ", Rj)
        print("dlam = ", dlam)
        print("logf = ", logf)
        print("logf1 = ", logf1)
        print("logf2 = ", logf2)
        print("logf3 = ", logf3)
        print(((0.01*np.min(obspec[2,:]**2))))
        print((100.*np.max(obspec[2,:]**2)))
        print("VMRs = ", invmr)
        print("ng = ", ng)
        print("sum VMRs = ", np.sum(10.**(invmr[0:ng])))
        print("cloud_top = ", cloud_top)
        print("cloud_bot = ", cloud_bot)
        print("cloud_height = ", cloud_height)
        print("loga = ", loga)
        print("b = ", b)
        print("cloud_tau0 = ", cloud_tau0)
        print("w0 = ", w0)
        print("taupow = ", taupow)
        print("pcover = ", pcover)
        print("sum(pcover) = ", np.sum(pcover))
        print("metscale = ", metscale)
        print("coscale = ", coscale)
        print("press[press.size-1]  = ", press[press.size-1])
        print("press[0] = ",press[0]) 
        if (all(invmr[0:ng] > -12.0)
            and all(invmr[0:ng] < 0.0)
            and (np.sum(10.**(invmr[0:ng])) < 1.0)
            and all(pcover > 0.) and (np.sum(pcover) == 1.0)
            and (metscale[0] <= mh <= metscale[-1])
            and (coscale[0] <= co <= coscale[-1])
            and (0.0 < logg < 6.0)
            and (1.0 < M < 80.) 
            and (0. < r2d2 < 1.)
            and (0.1 < scale1 < 10.0)
            and (0.1 < scale2 < 10.0)
            and (0.5 < Rj < 2.0)
            and (-0.01 < dlam < 0.01) 
            and ((0.01*np.min(obspec[2,:]**2)) < 10.**logf
                 < (100.*np.max(obspec[2,:]**2)))
            and ((0.01*np.min(obspec[2,s1]**2)) < 10.**logf1
                 < (100.*np.max(obspec[2,s1]**2)))
            and ((0.01*np.min(obspec[2,s2]**2)) < 10.**logf2
                 < (100.*np.max(obspec[2,s2]**2)))
            and ((0.01*np.min(obspec[2,s3]**2)) < 10.**logf3
                 < (100.*np.max(obspec[2,s3]**2)))
            and (np.all(cloud_tau0 >= 0.0))
            and (np.all(cloud_tau0 <= 100.0))
            and np.all(cloud_top < cloud_bot)
            and np.all(cloud_bot <= np.log10(press[press.size-1]))
            and np.all(np.log10(press[0]) <= cloud_top)
            and np.all(cloud_top < cloud_bot)
            and np.all(0. < cloud_height)
            and np.all(cloud_height < 7.0)
            and np.all(0.0 < w0)
            and np.all(w0 <= 1.0)
            and np.all(-10.0 < taupow)
            and np.all(taupow < +10.0)
            and np.all( -3.0 < loga)
            and np.all (loga < 3.0)
            and np.all(b < 1.0)
            and np.all(b > 0.0)):
            return 0.0

        return -np.inf

def lnlike(theta):

    gases_myP,chemeq,dist,cloudtype, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge,prof,do_bff,bff_raw,ceTgrid,metscale,coscale = settings.runargs

    #intemp, invmr, pcover, cloudtype, cloudparams, r2d2, logg, dlam, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,logf,proftype,do_fudge,do_bff,bff_raw,ceTgrid):


    # get the spectrum
    # for MCMC runs we don't want diagnostics
    gnostics = 0
    shiftspec, photspec,tauspec,cfunc = modelspec(theta,settings.runargs,gnostics)
    if chemeq == 0:
        if (gasnum[gasnum.size-1] == 21):
            ng = gasnum.size - 1
        elif (gasnum[gasnum.size-1] == 23):
            ng = gasnum.size -2
        else:
            ng = gasnum.size
        invmr = theta[0:ng]

    else:
        ng = 2

    if (fwhm < 0.0):
        if (fwhm == -1 or fwhm == -3):
            scale1 = theta[ng+2]
            scale2 = theta[ng+3]
            if (do_fudge == 1):
                logf = theta[ng+5:ng+8]
                nb = 8
            else:
                # This is a place holder value so the code doesn't break
                logf = np.log10(0.1*(max(obspec[2,10::3]))**2) 
                nb = 5
        elif (fwhm == -2):
            scale1 = theta[ng+2]
            if (do_fudge == 1):
                logf = theta[ng+4:ng+6]
                nb = 6
            else:
                # This is a place holder value so the code doesn't break
                logf = np.log10(0.1*(max(obspec[2,10::3]))**2) 
                nb = 4

    else:
        if (do_fudge == 1):
            logf = theta[ng+3]
            nb = 4
        else:
            # This is a place holder value so the code doesn't break
            logf = np.log10(0.1*(max(obspec[2,10::3]))**2) 
            nb = 3

    # If we've set a value for FWHM that we're using... 
    if (fwhm > 0.00 and fwhm < 1.00):
        # this is a uniform FWHM in microns
        
        modspec = conv_uniform_FWHM(shiftspec,obspec,fwhm)
   
        # Below is method for rebinning using conserve flux method
        #    oblen = obspec.shape[1]
        #    modspec = np.empty((2,oblen),dtype='d')
        #    modspec[1,:] =  rebinspec(spec[0,:], spec[1,:], obspec[0,:])
        # get log-likelihood
        # We've lifted this from Mike's code, below is original from emcee docs
        # Just taking every 3rd point to keep independence
    elif (fwhm == 0.0):
        # Use Mike's convolution for Spex
        wno = 1e4 / shiftspec[0,:]
        modspec = spex_non_uniform(obspec[0,:],wno,shiftspec[1,:])
    if (fwhm >= 0.0):
        if (do_fudge == 1):
            s2=obspec[2,::3]**2 + 10.**logf
        else:
            s2 = obspec[2,::3]**2

        lnLik=-0.5*np.sum((((obspec[1,::3] - modspec[::3])**2) / s2) + np.log(2.*np.pi*s2))
            
    elif (fwhm < 0.0):
        lnLik = 0.0
        # This is for multi-instrument cases
        # -1: spex + akari + IRS
        # -2: spex + IRS
        # -3: spex + Lband + IRS
        if (fwhm == -1):

            # Spex
            mr1 = np.where(shiftspec[0,:] < 2.5)
            or1  = np.where(obspec[0,:] < 2.5)
            wno = 1e4 / shiftspec[0,mr1]
            spec1 = spex_non_uniform(obspec[0,or1],wno,shiftspec[1,mr1])

            modspec = np.array([shiftspec[0,::-1],shiftspec[1,::-1]])
            # AKARI IRC
            # dispersion constant across order 0.0097um
            # R = 100 at 3.6um for emission lines
            # dL ~constant at 3.6 / 120
            dL = 0.03
            mr2 = np.where(np.logical_and(modspec[0,:] > 2.5,modspec[0,:] < 5.0))
            or2 = np.where(np.logical_and(obspec[0,:] > 2.5,obspec[0,:] < 5.0))
            spec2 = scale1 * conv_uniform_FWHM(modspec[:,mr2],obspec[:,or2],dL)

            # Spitzer IRS
            # R roughly constant within orders, and orders both appear to
            # have R ~ 100
            R = 100.0
            mr3 = np.where(modspec[0,:] > 5.0)
            or3 = np.where(obspec[0,:] > 5.0)
            spec3 = scale2 * conv_uniform_R(modspec[:,mr3],obspec[:,or3],R)

            if (do_fudge == 1):
                s1 = obspec[2,or1]**2 + 10.**logf[0]
                s2 = obspec[2,or2]**2 + 10.**logf[1]
                s3 = obspec[2,or3]**2 + 10.**logf[2]
            else:
                s1 = obspec[2,or1]**2
                s2 = obspec[2,or2]**2
                s3 = obspec[2,or3]**2


            lnLik1=-0.5*np.sum((((obspec[1,or1] - spec1)**2) / s1) + np.log(2.*np.pi*s1))
            lnLik2=-0.5*np.sum((((obspec[1,or2] - spec2)**2) / s2) + np.log(2.*np.pi*s2))
            lnLik3=-0.5*np.sum((((obspec[1,or3] - spec3)**2) / s3) + np.log(2.*np.pi*s3))
            lnLik = lnLik1 + lnLik2 + lnLik3

        elif (fwhm == -2):
            # This is just spex + IRS
            # Spex
            mr1 = np.where(shiftspec[0,:] < 2.5)
            or1  = np.where(obspec[0,:] < 2.5)
            wno = 1e4 / shiftspec[0,mr1]
            spec1 = spex_non_uniform(obspec[0,or1],wno,shiftspec[1,mr1])

            modspec = np.array([shiftspec[0,::-1],shiftspec[1,::-1]])

            # Spitzer IRS
            # R roughly constant within orders, and orders both appear to
            # have R ~ 100
            R = 100.0
            mr3 = np.where(modspec[0,:] > 5.0)
            or3 = np.where(obspec[0,:] > 5.0)
            spec3 = scale1 * conv_uniform_R(modspec[:,mr3],obspec[:,or3],R)

            if (do_fudge == 1):
                s1 = obspec[2,or1]**2 + 10.**logf[0]
                s3 = obspec[2,or3]**2 + 10.**logf[1]
            else:
                s1 = obspec[2,or1]**2
                s3 = obspec[2,or3]**2


            lnLik1=-0.5*np.sum((((obspec[1,or1] - spec1)**2) / s1) + np.log(2.*np.pi*s1))
            lnLik3=-0.5*np.sum((((obspec[1,or3] - spec3)**2) / s3) + np.log(2.*np.pi*s3))
            lnLik = lnLik1 + lnLik3
            
        elif (fwhm == -3):
            # This is spex + Mike Cushing's L band R = 425 + IRS
            # Spex
            mr1 = np.where(shiftspec[0,:] < 2.5)
            or1  = np.where(obspec[0,:] < 2.5)
            wno = 1e4 / shiftspec[0,mr1]
            spec1 = spex_non_uniform(obspec[0,or1],wno,shiftspec[1,mr1])

            modspec = np.array([shiftspec[0,::-1],shiftspec[1,::-1]])
            # Mike Cushing supplied L band R = 425 
            # dispersion constant across order 0.0097um
            # R = 425
            R = 425
            mr2 = np.where(np.logical_and(modspec[0,:] > 2.5,modspec[0,:] < 5.0))
            or2 = np.where(np.logical_and(obspec[0,:] > 2.5,obspec[0,:] < 5.0))
            spec2 = scale1 * conv_uniform_R(modspec[:,mr2],obspec[:,or2],R)

            # Spitzer IRS
            # R roughly constant within orders, and orders both appear to
            # have R ~ 100
            R = 100.0
            mr3 = np.where(modspec[0,:] > 5.0)
            or3 = np.where(obspec[0,:] > 5.0)
            spec3 = scale2 * conv_uniform_R(modspec[:,mr3],obspec[:,or3],R)

            if (do_fudge == 1):
                s1 = obspec[2,or1]**2 + 10.**logf[0]
                s2 = obspec[2,or2]**2 + 10.**logf[1]
                s3 = obspec[2,or3]**2 + 10.**logf[2]
            else:
                s1 = obspec[2,or1]**2
                s2 = obspec[2,or2]**2
                s3 = obspec[2,or3]**2


            lnLik1=-0.5*np.sum((((obspec[1,or1] - spec1)**2) / s1) + np.log(2.*np.pi*s1))
            lnLik2=-0.5*np.sum((((obspec[1,or2] - spec2)**2) / s2) + np.log(2.*np.pi*s2))
            lnLik3=-0.5*np.sum((((obspec[1,or3] - spec3)**2) / s3) + np.log(2.*np.pi*s3))
            lnLik = lnLik1 + lnLik2 + lnLik3
            
    return lnLik


def modelspec(theta, args,gnostics):

    gases_myP,chemeq,dist,cloudtype, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge,prof,do_bff,bff_raw,ceTgrid,metscale,coscale = args
    nlayers = press.size
    if chemeq == 0:
        if (gasnum[gasnum.size-1] == 21):
            ng = gasnum.size - 1
        elif (gasnum[gasnum.size-1] == 23):
            ng = gasnum.size -2
        else:
            ng = gasnum.size
        invmr = theta[0:ng]

    else:
        mh  = theta[0]
        co = theta[1]
        ng = 2
        mfit = interp1d(metscale,gases_myP,axis=0)
        gases_myM = mfit(mh)
        cfit = interp1d(coscale,gases_myM,axis=0)
        invmr = cfit(co)

    logg = theta[ng]
    if (fwhm < 0.0):
        if (fwhm == -1 or fwhm == -3):
            r2d2 = theta[ng+1:ng+4]
            dlam = theta[ng+4]
            if (do_fudge == 1):
                logf = theta[ng+5:ng+8]
                nb = 8
            else:
                # This is a place holder value so the code doesn't break
                logf = np.log10(0.1*(max(obspec[2,10::3]))**2) 
                nb = 5
        elif (fwhm == -2):
            r2d2 = theta[ng+1:ng+3]
            dlam = theta[ng+3]
            if (do_fudge == 1):
                logf = theta[ng+4:ng+6]
                nb = 6
            else:
                # This is a place holder value so the code doesn't break
                logf = np.log10(0.1*(max(obspec[2,10::3]))**2) 
                nb = 4

    else:
        r2d2 = theta[ng+1]
        dlam = theta[ng+2]
        if (do_fudge == 1):
            logf = theta[ng+3]
            nb = 4
        else:
            # This is a place holder value so the code doesn't break
            logf = np.log10(0.1*(max(obspec[2,10::3]))**2) 
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
    
    if (proftype == 1):
        gam = theta[pc+nc]
        intemp = theta[pc+1+nc:]
    elif (proftype == 2 or proftype ==3):
        intemp = theta[pc+nc:]
    elif (proftype == 9):
        intemp = prof
    else:
        raise ValueError("not valid profile type %proftype" % (char, string))

    # set the profile
    temp = TPmod.set_prof(proftype,coarsePress,press,intemp)

    ngas = gasnum.size
    bff = np.zeros([3,nlayers],dtype="float64")

    # check if its a fixed VMR or a profile from chem equilibrium 
    # VMR is log10(VMR) !!!
    if chemeq == 1:
        # this case is a profile
        ng = invmr.shape[2]
        ngas = ng - 3
        logVMR = np.zeros([ngas,nlayers],dtype='d')
        for p in range(0,nlayers):
            for g in range(0,ng):
                tfit = InterpolatedUnivariateSpline(ceTgrid,invmr[:,p,g])
                if (g < 3):
                    bff[g,p] = tfit(temp[p])
                else:
                    logVMR[g-3,p]= tfit(temp[p])
    else:
        # This case is fixed VMR
        # chemeq = 0
        logVMR = np.empty((ngas,nlayers),dtype='d')
        alkratio = 16.2 #  from Asplund et al (2009)

        # now sort Na and K
        # get the ngas for forward model (ngas, not ng
        if (gasnum[gasnum.size-1] == 21):
            ngas = invmr.shape[0] + 1
        elif (gasnum[gasnum.size-1] == 23):
            ngas = invmr.shape[0] + 2
        else:
            ngas = invmr.shape[0]

        tmpvmr = np.empty(ngas,dtype='d')
        if (gasnum[gasnum.size-1] == 21):
            tmpvmr[0:(ngas-2)] = invmr[0:(ngas-2)]
            tmpvmr[ngas-2] = np.log10(10.**invmr[ngas-2] / (alkratio+1.)) # K
            tmpvmr[ngas-1] = np.log10(10.**invmr[ngas-2] * (alkratio / (alkratio+1.))) # Na
        elif (gasnum[gasnum.size-1] == 23):
            #f values are ratios between Na and (K+Cs) and K and Cs respectively
            f1 = 1.348
            f2 = 8912.5
            tmpvmr[0:(ngas-3)] = invmr[0:(ngas-3)]
            tmpvmr[ngas-1] = np.log10(10.**invmr[ngas-3] / ((f1+1)*(f2+1))) # Cs
            tmpvmr[ngas-2] = np.log10(10.**invmr[ngas-3] * (f1 /(f1+1)) ) # Na
            tmpvmr[ngas-3] = np.log10(10.**invmr[ngas-3] - 10.**tmpvmr[ngas-2] - 10.**tmpvmr[ngas-1]) #K   
        else:
            tmpvmr[0:ngas] = invmr[0:ngas]
            
        for i in range(0,ngas):                              
            logVMR[i,:] = tmpvmr[i]

    # now need to translate cloudparams in to cloud profile even
    # if do_clouds is zero..

    cloudprof,cloudrad,cloudsig = cloud.atlas(do_clouds,cloudnum,cloudtype,cloudparams,press)
    cloudprof = np.asfortranarray(cloudprof,dtype = 'float64')
    cloudrad = np.asfortranarray(cloudrad,dtype = 'float64')
    cloudsig = np.asfortranarray(cloudsig,dtype = 'float64')
    pcover = np.asfortranarray(pcover,dtype = 'float32')
    cloudnum = np.asfortranarray(cloudnum,dtype='i')
    do_clouds = np.asfortranarray(do_clouds,dtype = 'i')

    # get r2d2 sorted for multi-instruments
    if (fwhm < 0.0):
        if (fwhm == -1 or fwhm == -3):
            R2D2 = r2d2[0]
            scale1 = r2d2[1]
            scale2 = r2d2[2]
        elif (fwhm == -2):
            R2D2 = r2d2[0]
            scale1 = r2d2[1]            
    else:
        R2D2 = r2d2
    

        
    # Now get the BFF stuff sorted
    if (chemeq == 0 and do_bff == 1):
        for gas in range(0,3):
            for i in range(0,nlayers):
                tfit = InterpolatedUnivariateSpline(ceTgrid,bff_raw[:,i,gas],k=1) 
                bff[gas,i] = tfit(temp[i])

    bff = np.asfortranarray(bff, dtype='float64')
    press = np.asfortranarray(press,dtype='float32')
    temp = np.asfortranarray(temp,dtype='float64')
    logVMR = np.asfortranarray(logVMR,dtype='float64')

    # Diagnostics below.
    # make_cf = get a contribution function
    # clphot = get pressure for cloud_tau = 1.0 as function of wavelength
    # ^^ i.e the cloud photosphere
    # ophot = get pressures for tau(not cloud) = 1.0 as function of wavelength]
    # ^^ i.e. the photosphere due to other (gas phase) opacities)
    
    # Set clphot,ophot and cfunc as we don't need these in the emcee run
    if (gnostics == 0):
        clphot = 0
        ophot = 0
        make_cf = 0
    else:
        clphot = 1
        ophot = 1
        make_cf = 1

    # now we can call the forward model
    outspec,tmpclphotspec,tmpophotspec,cf = forwardmodel.marv(temp,logg,R2D2,gasnum,logVMR,pcover,do_clouds,cloudnum,cloudrad,cloudsig,cloudprof,inlinetemps,press,inwavenum,linelist,cia,ciatemps,use_disort,clphot,ophot,make_cf,do_bff,bff)

    # Trim to length where it is defined.
    nwave = inwavenum.size
    trimspec = np.zeros([2,nwave],dtype='d')
    trimspec = outspec[:,:nwave]
    cloud_phot_press = tmpclphotspec[0:npatches,:nwave].reshape(npatches,nwave)
    other_phot_press = tmpophotspec[0:npatches,:nwave].reshape(npatches,nwave)
    cfunc = np.zeros([npatches,nwave,nlayers],dtype='d')
    cfunc = cf[:npatches,:nwave,:nlayers].reshape(npatches,nwave,nlayers)

    # now shift wavelen by delta_lambda
    shiftspec = np.empty_like(trimspec)
    shiftspec[0,:] =  trimspec[0,:] + dlam
    shiftspec[1,:] =  trimspec[1,:]


    return shiftspec, cloud_phot_press,other_phot_press,cfunc
