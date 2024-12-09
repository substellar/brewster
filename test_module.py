#!/usr/bin/env python

""" Module of bits to define prior,lnlike,and lnprob for MCMC/Multinest """
from __future__ import print_function
import time
import math
import numpy as np
import scipy as sp
import gc
import ciamod
import TPmod
import os
import sys
import pickle
import forwardmodel
import cloud_dic
from builtins import str
from builtins import range
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.convolution import convolve, convolve_fft
from astropy.convolution import Gaussian1DKernel
from bensconv import prism_non_uniform
from bensconv import conv_uniform_R
from bensconv import conv_uniform_FWHM
from collections import namedtuple
import utils
import settings
import gas_nonuniform


__author__ = "Fei Wang"
__copyright__ = "Copyright 2024 - Fei Wang"
__credits__ = ["Fei Wang", "Ben Burningham"]
__license__ = "GPL"
__version__ = "0.2"  
__maintainer__ = ""
__email__ = ""
__status__ = "Development"

def lnprob(theta,re_params):


    # now check against the priors, if not beyond them, run the likelihood
    lp = lnprior(theta,re_params)
    if not np.isfinite(lp):
        return -np.inf
    # run the likelihood
    lnlike_value = lnlike(theta,re_params)

    lnprb = lp+lnlike_value
    if np.isnan(lnprb):
        lnprb = -np.inf
    return lnprb



def lnprior(theta,re_params):


    all_params,all_params_values =utils.get_all_parametres(re_params.dictionary) 
    params_master = namedtuple('params',all_params)
    params_instance = params_master(*theta)

    args_instance=settings.runargs

# Unpack all necessary parameters into local variables
    (   gases_myP,
        chemeq,
        dist,
        cloudtype,
        do_clouds,
        gaslist,
        gasnames,
        gasmass,
        cloudnum,
        inlinetemps,
        coarsePress,
        press,
        inwavenum,
        linelist,
        cia,
        ciatemps,
        use_disort,
        fwhm,
        obspec,
        proftype,
        do_fudge,
        prof,
        do_bff,
        bff_raw,
        ceTgrid,
        metscale,
        coscale
    ) = (
        args_instance.gases_myP,
        args_instance.chemeq,
        args_instance.dist,
        args_instance.cloudtype,
        args_instance.do_clouds,
        args_instance.gaslist,
        args_instance.gasnames,
        args_instance.gasmass,
        args_instance.cloudnum,
        args_instance.inlinetemps,
        args_instance.coarsePress,
        args_instance.press,
        args_instance.inwavenum,
        args_instance.linelist,
        args_instance.cia,
        args_instance.ciatemps,
        args_instance.use_disort,
        args_instance.fwhm,
        args_instance.obspec,
        args_instance.proftype,
        args_instance.do_fudge,
        args_instance.prof,
        args_instance.do_bff,
        args_instance.bff_raw,
        args_instance.ceTgrid,
        args_instance.metscale,
        args_instance.coscale
    )

    
    knots = len(coarsePress)
    ndim = len(theta)

    # set up the priors here
    if (chemeq != 0):
        invmr = np.array([-3.,-3.])
        mh = params_instance.mh
        co = params_instance.co

    else:
        gas_keys = re_params.dictionary['gas'].keys()
        gas_keys=list(gas_keys)
        invmr=np.array([getattr(params_instance, key) for key in gas_keys])
        mh = 0.0
        co = 1.0
    
    logg = params_instance.logg

    if (fwhm < 0.0):
        if (fwhm == -1 or fwhm == -3 or fwhm == -4):
            s1  = np.where(obspec[0,:] < 2.5)
            s2  = np.where(np.logical_and(obspec[0,:] > 2.5,obspec[0,:] < 5.0))
            s3 =  np.where(obspec[0,:] > 5.)
            r2d2 = params_instance.r2d2
            scale1 = params_instance.scale1
            scale2 = params_instance.scale2
            dlam = params_instance.dlambda
            if (do_fudge == 1):
                logf1 = params_instance.tolerance_parameter_1
                logf2 = params_instance.tolerance_parameter_2
                logf3 = params_instance.tolerance_parameter_3
                logf = np.log10(0.1*(max(obspec[2,10::3]))**2)
            else:
                # This is a place holder value so the code doesn't break
                logf = np.log10(0.1*(max(obspec[2,:]))**2)
                logf1 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf2 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf3 = np.log10(0.1*(max(obspec[2,:]))**2)

        elif (fwhm == -2):
            s1  = np.where(obspec[0,:] < 2.5)
            s2 = s1
            s3 =  np.where(obspec[0,:] > 5.)
            r2d2 = params_instance.r2d2
            scale1 = 1.0 # dummy value
            scale2 = params_instance.scale2
            dlam = params_instance.dlambda
            if (do_fudge == 1):
                logf1 = params_instance.tolerance_parameter_1
                logf2 = np.log10(0.1*(max(obspec[2,10::3]))**2) # dummy
                logf3 = params_instance.tolerance_parameter_2
                logf = np.log10(0.1*(max(obspec[2,10::3]))**2)
            else:
                # This is a place holder value so the code doesn't break
                logf = np.log10(0.1*(max(obspec[2,:]))**2)
                logf1 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf2 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf3 = np.log10(0.1*(max(obspec[2,:]))**2)
        elif (fwhm == -5):
            # this is for JWST NIRSpec prism + MIRI MRS
            # we assume that this is all scaled well
            # as is currently just simulated data anyhow
            # so this unpack is the same as a single instrument
            s1 = np.where(obspec[0,:] > 0.0)
            s2 = s1
            s3 = s1
            r2d2 = params_instance.r2d2
            dlam = params_instance.dlambda
            if (do_fudge == 1):
                logf =  params_instance.tolerance_parameter_1
                logf1 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf2 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf3 = np.log10(0.1*(max(obspec[2,:]))**2)
                scale1 = 1.0
                scale2 = 1.0
            else:
                # This is a place holder value so the code doesn't break
                logf = np.log10(0.1*(max(obspec[2,10::3]))**2)
                logf1 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf2 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf3 = np.log10(0.1*(max(obspec[2,:]))**2)
                scale1 = 1.0
                scale2 = 1.0

        elif (fwhm == -6):  ### UKIRT first and second order (Geballe cuts)
            s1  = np.where(obspec[0,:] < 1.585) ### wavelength less than 1.585um (second order)
            s2 = s1
            s3 =  np.where(obspec[0,:] > 1.585) ### wavelength greater than 1.585um (first order)
            r2d2 = params_instance.r2d2
            dlam = params_instance.dlambda
            scale1 = 1.0
            scale2 = 1.0
            if (do_fudge == 1):
                logf =  params_instance.tolerance_parameter_1
                logf1 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf2 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf3 = np.log10(0.1*(max(obspec[2,:]))**2)
            else:
                # This is a place holder value so the code doesn't break
                logf = np.log10(0.1*(max(obspec[2,10::3]))**2)
                logf1 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf2 = np.log10(0.1*(max(obspec[2,:]))**2)
                logf3 = np.log10(0.1*(max(obspec[2,:]))**2)
    
    else:
        # this just copes with normal, single instrument data
        s1 = np.where(obspec[0,:] > 0.0)
        s2 = s1
        s3 = s1
        
        r2d2 = params_instance.r2d2
        dlam = params_instance.dlambda
        if (do_fudge == 1):
            logf =  params_instance.tolerance_parameter_1
            logf1 = np.log10(0.1*(max(obspec[2,:]))**2)
            logf2 = np.log10(0.1*(max(obspec[2,:]))**2)
            logf3 = np.log10(0.1*(max(obspec[2,:]))**2)
            scale1 = 1.0
            scale2 = 1.0

        else:
            # This is a place holder value so the code doesn't break
            logf = np.log10(0.1*(max(obspec[2,10::3]))**2)
            logf1 = np.log10(0.1*(max(obspec[2,:]))**2)
            logf2 = np.log10(0.1*(max(obspec[2,:]))**2)
            logf3 = np.log10(0.1*(max(obspec[2,:]))**2)
            scale1 = 1.0
            scale2 = 1.0


    npatches = do_clouds.size
    if (npatches > 1):
        prat = params_instance.fcld
        pcover = np.array([prat,(1.-prat)])
    else:
        pcover = np.array([0.5,0.5])

    # use correct unpack method depending on situation


    # attribute_names = list(params_instance._fields)
    # first_cloud_para=list(re_params.dictionary['cloud']['patch 1']['params'].keys())[0]
    # # Find the index of pc
    # pc = attribute_names.index(first_cloud_para)

    if ((npatches > 1) and np.all(do_clouds != 0)):
        cloudparams = cloud_dic.unpack_patchy(re_params,params_instance,cloudtype,cloudnum,do_clouds)
    else:
        cloudparams = cloud_dic.unpack_default(re_params,params_instance,cloudtype,cloudnum,do_clouds)


    if (cloudtype.size > cloudtype.shape[1]):
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
                cloud_tau0[i,:] = 1.0
                cloud_bot[i,:] = np.log10(press[-1])
                cloud_top[i,:] = np.log10(press[0])
                cloud_height[i,:] = 0.1
                w0[i,:] = +0.5
                taupow[i,:] = 0.0
                loga[i,:] =  0.0
                b[i,:] = 0.5

                
                    
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


    intemp_keys = list(re_params.dictionary['pt']['params'].keys())
    intemp = np.array([getattr(params_instance, key) for key in intemp_keys])

    gastype_values = [info['gastype'] for key, info in re_params.dictionary['gas'].items() if 'gastype' in info]
    count_N = gastype_values.count('N')

    if (proftype == 1):
        gam = params_instance.gamma
        T = intemp[1:]  #theta[:-2]  the pressure  deck  of  nonuniform gas
        diff=np.roll(T,-1)-2.*T+np.roll(T,1)
        pp=len(T)

        gas_profile=-1
        P_hgas=0
        if count_N>0:
            gas_profile = np.full((count_N, press.size), -1.0)
            gas_profile_index =0
            for i in range(len(gastype_values)):
                if  gastype_values[i]=="N":
                    P_gas= getattr(params_instance, "p_ref_%s"%gas_keys[i])
                    gas_alpha= getattr(params_instance, "alpha_%s"%gas_keys[i])
                    t_gas= getattr(params_instance, gas_keys[i])
                    if (0. < gas_alpha < 1. and -12.0 < t_gas < 0.0  and np.log10(press[0]) <= P_gas <= 2.4):
                        gas_profile[gas_profile_index,:]=gas_nonuniform.non_uniform_gas(press,P_gas,t_gas,gas_alpha)
                    else:
                        gas_profile[gas_profile_index,:]=-30
                    gas_profile_index+=1

                if  gastype_values[i]=="H":
                      P_hgas= getattr(params_instance, "p_ref_%s"%gas_keys[i])


        #for mass prior
        D = 3.086e+16 * dist
        R = -1.0
        if (r2d2 > 0.):
            R = np.sqrt(r2d2) * D
        g = (10.**logg)/100.
        M = (R**2 * g/(6.67E-11))/1.898E27
        Rj = R / 69911.e3
        #         and  and (-5. < logbeta < 0))

        if (all(invmr > -12.0) and all(invmr < 0.0) and (np.sum(10.**(invmr)) < 1.0)
            and np.all(gas_profile > -25.0) and np.all(gas_profile < 0.0)
            and np.log10(press[0]) <= P_hgas <= 2.4
            and all(pcover > 0.) and (np.sum(pcover) == 1.0)
            and  metscale[0] <= mh <= metscale[-1]
            and  coscale[0] <= co <= coscale[-1]
            and  0.0 < logg < 6.0
            and  1.0 < M < 80
            and  0. < r2d2 < 1.
            and 0.1 < scale1 < 10.0
            and 0.1 < scale2 < 10.0
            and  0.5 < Rj < 2.0
            and -0.01 < dlam < 0.01
            and (min(T) > 1.0) and (max(T) < 6000.)
            and (0 < gam < 5000)
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
        a1,a2,P1,P3,T3=intemp[:]
        T = np.empty([press.size])
        T[:] = -100.

        gas_profile=-1
        P_hgas=0
        if count_N>0:
            gas_profile = np.full((count_N, press.size), -1.0)
            gas_profile_index =0
            for i in range(len(gastype_values)):
                if  gastype_values[i]=="N":
                    P_gas= getattr(params_instance, "p_ref_%s"%gas_keys[i])
                    gas_alpha= getattr(params_instance, "alpha_%s"%gas_keys[i])
                    t_gas= getattr(params_instance, gas_keys[i])
                    if (0. < gas_alpha < 1. and -12.0 < t_gas < 0.0  and np.log10(press[0]) <= P_gas <= 2.4):
                        gas_profile[gas_profile_index,:]=gas_nonuniform.non_uniform_gas(press,P_gas,t_gas,gas_alpha)
                    else:
                        gas_profile[gas_profile_index,:]=-30
                    gas_profile_index+=1

                if  gastype_values[i]=="H":
                      P_hgas= getattr(params_instance, "p_ref_%s"%gas_keys[i])
    

        if (0. < a1 < 1. and 0. < a2 < 1.0
            and T3 > 0.0 and P3 >= P1 and P1 >= np.log10(press[0])
            and P3 <= 5):
            T = TPmod.set_prof(proftype,junkP,press,intemp) #theta[:-2]  the pressure  deck  of  nonuniform gas
        #for mass prior
        D = 3.086e+16 * dist
        R = -1.0
        if (r2d2 > 0.):
            R = np.sqrt(r2d2) * D
        g = (10.**logg)/100.
        M = (R**2 * g/(6.67E-11))/1.898E27
        Rj = R / 69911.e3
        #         and  and (-5. < logbeta < 0))
        if (all(invmr > -12.0) and all(invmr < 0.0) and (np.sum(10.**(invmr)) < 1.0)
            and np.all(gas_profile > -25.0) and np.all(gas_profile < 0.0)
            and np.log10(press[0]) <= P_hgas <= 2.4
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
        a1,a2,P1,P2,P3,T3=intemp[:]

        gas_profile=-1
        P_hgas=0
        if count_N>0:
            gas_profile = np.full((count_N, press.size), -1.0)
            gas_profile_index =0
            for i in range(len(gastype_values)):
                if  gastype_values[i]=="N":
                    P_gas= getattr(params_instance, "p_ref_%s"%gas_keys[i])
                    gas_alpha= getattr(params_instance, "alpha_%s"%gas_keys[i])
                    t_gas= getattr(params_instance, gas_keys[i])
                    if (0. < gas_alpha < 1. and -12.0 < t_gas < 0.0  and np.log10(press[0]) <= P_gas <= 2.4):
                        gas_profile[gas_profile_index,:]=gas_nonuniform.non_uniform_gas(press,P_gas,t_gas,gas_alpha)
                    else:
                        gas_profile[gas_profile_index,:]=-30
                    gas_profile_index+=1

                if  gastype_values[i]=="H":
                      P_hgas= getattr(params_instance, "p_ref_%s"%gas_keys[i])
    

        T = np.empty([press.size])
        T[:] = -100.
        if  (0. < a1 < 1. and 0. < a2 < 1.0
            and T3 > 0.0 and P3 >= P2 and P3 >= P1 and P2 >= np.log10(press[0]) and P1 >= np.log10(press[0])
             and P3 <= 5):
            T = TPmod.set_prof(proftype,junkP,press,intemp) #theta[:-2]  the pressure  deck  of  nonuniform gas ? and P2 >= P1 


        #for mass prior
        D = 3.086e+16 * dist
        R = -1.0
        if (r2d2 > 0.):
            R = np.sqrt(r2d2) * D
        g = (10.**logg)/100.
        M = (R**2 * g/(6.67E-11))/1.898E27
        Rj = R / 69911.e3
        #         and  and (-5. < logbeta < 0))
        if (all(invmr > -12.0) and all(invmr < 0.0) and (np.sum(10.**(invmr)) < 1.0)
            and np.all(gas_profile > -25.0) and np.all(gas_profile < 0.0)
            and np.log10(press[0]) <= P_hgas <= 2.4
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

    elif (proftype == 7):
        Tint,alpha,lndelta,T1,T2,T3=intemp[:] #kappa/grav replaced by delta   kth ～(10^-2,10^-3) [cm^2g-1] g [cm s^-2] 1e4

        gas_profile=-1
        P_hgas=0
        if count_N>0:
            gas_profile = np.full((count_N, press.size), -1.0)
            gas_profile_index =0
            for i in range(len(gastype_values)):
                if  gastype_values[i]=="N":
                    P_gas= getattr(params_instance, "p_ref_%s"%gas_keys[i])
                    gas_alpha= getattr(params_instance, "alpha_%s"%gas_keys[i])
                    t_gas= getattr(params_instance, gas_keys[i])
                    if (0. < gas_alpha < 1. and -12.0 < t_gas < 0.0  and np.log10(press[0]) <= P_gas <= 2.4):
                        gas_profile[gas_profile_index,:]=gas_nonuniform.non_uniform_gas(press,P_gas,t_gas,gas_alpha)
                    else:
                        gas_profile[gas_profile_index,:]=-30
                    gas_profile_index+=1

                if  gastype_values[i]=="H":
                      P_hgas= getattr(params_instance, "p_ref_%s"%gas_keys[i])
    
        delta=np.exp(lndelta)
        Tconnect = (((3/4) * Tint**4) * ((2/3) + (0.1)))**(1/4)
        T = np.empty([press.size])
        T[:] = -100.
        # if  (1. < alpha  < 2. and 0. < delta < 0.1
        #     and T1 > 0.0 and T1 < T2 and T2 < T3 and T3 < Tconnect and Tint >0.0):
        #     T = TPmod.set_prof(proftype,junkP,press,theta[pc+nc:-2]) #inversion not allowed, theta[:-2]  the pressure  deck  of  nonuniform gas
        # put prior on P1 to put it shallower than 100 bar 
        P1 = ((1/delta)**(1/alpha))
        if  (1 < alpha  < 2. and P1 < 100 and P1 > press[0]
            and T1 > 0.0 and T2 > 0.0 and T3 > 0.0 and Tint >0.0):
            T = TPmod.set_prof(proftype,junkP,press,intemp) #allow inversion theta[:-2]  the pressure  deck  of  nonuniform gas

        #for mass prior
        D = 3.086e+16 * dist
        R = -1.0
        if (r2d2 > 0.):
            R = np.sqrt(r2d2) * D
        g = (10.**logg)/100.
        M = (R**2 * g/(6.67E-11))/1.898E27
        Rj = R / 69911.e3
        #         and  and (-5. < logbeta < 0))
        if (all(invmr > -12.0) and all(invmr < 0.0) and (np.sum(10.**(invmr)) < 1.0)
            and np.all(gas_profile > -25.0) and np.all(gas_profile < 0.0)
            and np.log10(press[0]) <= P_hgas <= 2.4
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

    elif (proftype == 77):
        Tint,alpha,lndelta,T1,T2,T3=intemp[1:] #kappa/grav replaced by delta   kth ～(10^-2,10^-3) [cm^2g-1] g [cm s^-2] 1e4

        gas_profile=-1
        P_hgas=0
        if count_N>0:
            gas_profile = np.full((count_N, press.size), -1.0)
            gas_profile_index =0
            for i in range(len(gastype_values)):
                if  gastype_values[i]=="N":
                    P_gas= getattr(params_instance, "p_ref_%s"%gas_keys[i])
                    gas_alpha= getattr(params_instance, "alpha_%s"%gas_keys[i])
                    t_gas= getattr(params_instance, gas_keys[i])
                    if (0. < gas_alpha < 1. and -12.0 < t_gas < 0.0  and np.log10(press[0]) <= P_gas <= 2.4):
                        gas_profile[gas_profile_index,:]=gas_nonuniform.non_uniform_gas(press,P_gas,t_gas,gas_alpha)
                    else:
                        gas_profile[gas_profile_index,:]=-30
                    gas_profile_index+=1

                if  gastype_values[i]=="H":
                      P_hgas= getattr(params_instance, "p_ref_%s"%gas_keys[i])
    

        delta= np.exp(lndelta)
        Tconnect = (((3/4) * Tint**4) * ((2/3) + (0.1)))**(1/4)

        T = np.empty([press.size])
        T[:] = -100.

         #if  (1 < alpha  < 2. and 0. < delta < 0.1
          #   and T1 > 0.0 and T1 < T2 and T2 < T3 and T3 < Tconnect and Tint >0.0):
           #  T = TPmod.set_prof(proftype,junkP,press,theta[pc+nc:])  # no inversion 
        P1 = ((1/delta)**(1/alpha))
        # put prior on P1 to put it shallower than 100 bar   
        if  (1 < alpha  < 2. and P1 < 100. and P1 > press[0]
             and T1 > 0.0 and T2 > 0.0 and T3 > 0.0 and Tint >0.0):
            T = TPmod.set_prof(proftype,junkP,press,intemp) # allow inversion 

        # bits for smoothing in prior
        gam = params_instance.gamma
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

        if (all(invmr > -12.0) and all(invmr < 0.0) and (np.sum(10.**(invmr)) < 1.0)
            and  np.all(gas_profile > -25.0) and np.all(gas_profile < 0.0)
            and np.log10(press[0]) <= P_hgas <= 2.4
            and  all(pcover > 0.) and (np.sum(pcover) == 1.0)
            and  metscale[0] <=  mh <= metscale[-1]
            and  coscale[0] <= co <= coscale[-1]
            and  0.0 < logg < 6.0
            and 1.0 < M < 80.
            and  0. < r2d2 < 1.
            and 0.1 < scale1 < 10.0
            and 0.1 < scale2 < 10.0
            and  0.5 < Rj < 2.0
            # and -250 < vrad < 250
            # and 0. < vsini < 100.0
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
            and (gam > 0.)
            and  (min(T) > 1.0) and (max(T) < 6000.)):
            logbeta = -5.0
            beta=10.**logbeta
            alpha=1.0
            x=gam
            invgamma=((beta**alpha)/math.gamma(alpha)) * (x**(-alpha-1)) * np.exp(-beta/x)
            prprob = (-0.5/gam)*np.sum(diff[1:-1]**2) - 0.5*pp*np.log(gam) + np.log(invgamma)
    
            return prprob
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
        print("sum VMRs = ", np.sum(10.**(invmr)))
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
        if (all(invmr> -12.0)
            and all(invmr < 0.0)
            and (np.sum(10.**(invmr)) < 1.0)
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




def priormap_dic(theta,re_params):

    all_params,all_params_values =utils.get_all_parametres(re_params.dictionary) 
    params_master = namedtuple('params',all_params)
    params_instance = params_master(*theta)

    args_instance=settings.runargs
# Unpack all necessary parameters into local variables
    (   gases_myP,
        chemeq,
        dist,
        cloudtype,
        do_clouds,
        gaslist,
        gasnames,
        gasmass,
        cloudnum,
        inlinetemps,
        coarsePress,
        press,
        inwavenum,
        linelist,
        cia,
        ciatemps,
        use_disort,
        fwhm,
        obspec,
        proftype,
        do_fudge,
        prof,
        do_bff,
        bff_raw,
        ceTgrid,
        metscale,
        coscale
    ) = (
        args_instance.gases_myP,
        args_instance.chemeq,
        args_instance.dist,
        args_instance.cloudtype,
        args_instance.do_clouds,
        args_instance.gaslist,
        args_instance.gasnames,
        args_instance.gasmass,
        args_instance.cloudnum,
        args_instance.inlinetemps,
        args_instance.coarsePress,
        args_instance.press,
        args_instance.inwavenum,
        args_instance.linelist,
        args_instance.cia,
        args_instance.ciatemps,
        args_instance.use_disort,
        args_instance.fwhm,
        args_instance.obspec,
        args_instance.proftype,
        args_instance.do_fudge,
        args_instance.prof,
        args_instance.do_bff,
        args_instance.bff_raw,
        args_instance.ceTgrid,
        args_instance.metscale,
        args_instance.coscale
    )

    phi = np.zeros_like(theta)
    gaslist=list(re_params.dictionary["gas"].keys())
    gastype_values = [info['gastype'] for key, info in re_params.dictionary['gas'].items() if 'gastype' in info]

    gaspara=[]
    for i in range(len(gaslist)):
        gaspara.append(gaslist[i])
        if  gastype_values[i]=='N':
            gaspara.append("p_ref_%s"%gaslist[i])
            gaspara.append("alpha_%s"%gaslist[i])

        if  gastype_values[i]=='H':
            gaspara.append("p_ref_%s"%gaslist[i])
             
    ng=len(gaspara)
    

    if ng==2:
        phi[0] = (theta[0] * (metscale[-1] - metscale[0])) + metscale[0]
        phi[1] = (theta[1] * (coscale[-1] - coscale[0])) + coscale[0]
        
    else:
        rem = 1
        for i in range(0, ng):
            if gaspara[i] in gaslist:
                phi[i] = np.log10(rem) -  (theta[i] * 12.)
                rem = rem - (10**phi[i])
            elif gaspara[i].startswith('p_ref'):
                phi[i]= (theta[i]* \
                             (np.log10(press[-1]) - np.log10(press[0]))) + np.log10(press[0])

            elif gaspara[i].startswith('alpha'):
                phi[i]= theta[i]

    max_mass = 80. # jupiters
    min_mass = 1.0 # jupiters
    min_rad = 0.5 # jupiters
    max_rad = 2.5 # jupiters
    
    
    mass_index=params_instance._fields.index('M')
    
    # this is a simple uniform prior on mass
    # we want to use the radius, to set a mass prior. 
    # this will correlate these parameters??? Yes. which is fine.
    phi[mass_index] = (theta[mass_index] * (max_mass - min_mass)) + min_mass

    # this is if we want log g prior: phi[ng] = theta[ng] * 5.5
    # now we retrieve radius in R_jup 
    R_index=params_instance._fields.index('R')
    R_j = ((max_rad - min_rad)*theta[R_index]) + min_rad
    phi[R_index] = R_j
    
    
    if (fwhm < 0.0):
        if (fwhm == -1 or fwhm == -3 or fwhm == -4 or fwhm == -8):
            s1  = np.where(obspec[0,:] < 2.5)
            s2  = np.where(np.logical_and(obspec[0,:] > 2.5,obspec[0,:] < 5.0))
            s3 =  np.where(obspec[0,:] > 5.)
            # scale parameters here - generous factor 2 either side??
            scale1_index=params_instance._fields.index('scale1')
            scale2_index=params_instance._fields.index('scale2')
            phi[scale1_index] = (theta[scale1_index] * 1.5) + 0.5
            phi[scale2_index] = (theta[scale2_index] * 1.5) + 0.5
            # now dlam
            dlam_index=params_instance._fields.index('dlambda')
            phi[dlam_index] = (theta[dlam_index] * 0.02) - 0.01
            if (do_fudge == 1):
                # These are tolerances for each band due to difference SNRs
                minerr_s1 =np.log10((0.01 *  np.min(obspec[2,s1]))**2.)
                maxerr_s1 =np.log10((100.*np.max(obspec[2,s1]))**2.)
                tolerance_parameter_1_index=params_instance._fields.index('tolerance_parameter_1')
                phi[tolerance_parameter_1_index] = (theta[tolerance_parameter_1_index] * (maxerr_s1 - minerr_s1)) + minerr_s1
                
                minerr_s2 =np.log10((0.01 *  np.min(obspec[2,s2]))**2.)
                maxerr_s2 =np.log10((100.*np.max(obspec[2,s2]))**2.)
                tolerance_parameter_2_index=params_instance._fields.index('tolerance_parameter_2')
                phi[tolerance_parameter_2_index] = (theta[tolerance_parameter_2_index] * (maxerr_s2 - minerr_s2)) + minerr_s2
                
                minerr_s3 =np.log10((0.01 *  np.min(obspec[2,s3]))**2.)
                maxerr_s3 = np.log10((100.*np.max(obspec[2,s3]))**2.)
                tolerance_parameter_3_index=params_instance._fields.index('tolerance_parameter_3')
                phi[tolerance_parameter_3_index] = (theta[tolerance_parameter_3_index] * (maxerr_s3 - minerr_s3)) + minerr_s3

        elif (fwhm == -2):
            s1  = np.where(obspec[0,:] < 2.5)
            s3 =  np.where(obspec[0,:] > 5.)
            # scale parameter
            scale1_index=params_instance._fields.index('scale1')
            phi[scale1_index] = (theta[scale1_index] * 1.5) + 0.5
            # dlam now:
            dlam_index=params_instance._fields.index('dlambda')
            phi[dlam_index] = (theta[dlam_index] * 0.02) - 0.01
            if (do_fudge == 1):
                # These are tolerances for each band due to difference SNR
                minerr_s1 = np.log10((0.01 *  np.min(obspec[2,s1]))**2.)
                maxerr_s1 = np.log10((100.*np.max(obspec[2,s1]))**2.)
                tolerance_parameter_1_index=params_instance._fields.index('tolerance_parameter_1')
                phi[tolerance_parameter_1_index] = (theta[tolerance_parameter_1_index] * (maxerr_s1 - minerr_s1)) + minerr_s1
                
                minerr_s3 = np.log10((0.01 *  np.min(obspec[2,s3]))**2.)
                maxerr_s3 = np.log10((100.*np.max(obspec[2,s3]))**2.)
                tolerance_parameter_2_index=params_instance._fields.index('tolerance_parameter_2')
                phi[tolerance_parameter_2_index] = (theta[tolerance_parameter_2_index] * (maxerr_s3 - minerr_s3)) + minerr_s3

        elif (fwhm == -6):
            ##Geballe NIR CGS4 data
            s1  = np.where(obspec[0,:] < 1.585)
            s3  = np.where(obspec[0,:] > 1.585)
            #not including relative scale factor since data is calibrated to the same photometry
            #dlam:
            dlam_index=params_instance._fields.index('dlambda')
            phi[dlam_index] = (theta[dlam_index] * 0.02) - 0.01
            #Tolerance parameter (only one):
            if (do_fudge==1):
                minerr = np.log10((0.01 *  np.min(obspec[2,:]))**2.)
                maxerr = np.log10((100.*np.max(obspec[2,:]))**2.)
                tolerance_parameter_1_index=params_instance._fields.index('tolerance_parameter_1')
                phi[tolerance_parameter_1_index] = (theta[tolerance_parameter_1_index] * (maxerr - minerr)) + minerr

        elif (fwhm == -7): #Geballe NIR + NIRC + CGS4 MIR
            s1 = np.where(obspec[0, :] < 1.585)
            s2 = np.where(obspec[0, :] > 1.585)
            s3 = np.where(np.logical_and(obspec[0, :] > 2.52, obspec[0, :] < 4.2))  #NIRC
            s4 = np.where(obspec[0, :] > 4.2) #CGS4
            # scale parameter
            scale1_index=params_instance._fields.index('scale1')
            scale2_index=params_instance._fields.index('scale2')
            phi[scale1_index] = (theta[scale1_index] * 1.5) + 0.5
            phi[scale2_index] = (theta[scale2_index] * 1.5) + 0.5
            #dlam
            dlam_index=params_instance._fields.index('dlambda')
            phi[dlam_index] = (theta[dlam_index] * 0.02) - 0.01
            if (do_fudge == 1):
                # These are tolerances for each band due to difference SNRs
                minerr_s1 =np.log10((0.01 *  np.min(obspec[2,s1]))**2.)
                maxerr_s1 =np.log10((100.*np.max(obspec[2,s1]))**2.)
                tolerance_parameter_1_index=params_instance._fields.index('tolerance_parameter_1')
                phi[tolerance_parameter_1_index] = (theta[tolerance_parameter_1_index] * (maxerr_s1 - minerr_s1)) + minerr_s1
                
                minerr_s2 =np.log10((0.01 *  np.min(obspec[2,s2]))**2.)
                maxerr_s2 =np.log10((100.*np.max(obspec[2,s2]))**2.)
                tolerance_parameter_2_index=params_instance._fields.index('tolerance_parameter_2')
                phi[tolerance_parameter_2_index] = (theta[tolerance_parameter_2_index] * (maxerr_s2 - minerr_s2)) + minerr_s2
                
                minerr_s3 =np.log10((0.01 *  np.min(obspec[2,s3]))**2.)
                maxerr_s3 = np.log10((100.*np.max(obspec[2,s3]))**2.)
                tolerance_parameter_3_index=params_instance._fields.index('tolerance_parameter_3')
                phi[tolerance_parameter_3_index] = (theta[tolerance_parameter_3_index] * (maxerr_s3 - minerr_s3)) + minerr_s3


    else:
        # this just copes with normal, single instrument data
        # so do dlam next
        dlam_index=params_instance._fields.index('dlambda')
        phi[dlam_index] = (theta[dlam_index] * 0.02) - 0.01
        # now fudge
        if (do_fudge == 1):
            # logf here
            minerr =np.log10((0.01 *  np.min(obspec[2,:]))**2.)
            maxerr = np.log10((100.*np.max(obspec[2,:]))**2.)
            tolerance_parameter_1_index=params_instance._fields.index('tolerance_parameter_1')
            phi[tolerance_parameter_1_index] = (theta[tolerance_parameter_1_index] * (maxerr - minerr)) + minerr

 
    if (proftype == 1):

        intemp_keys = list(re_params.dictionary['pt']['params'].keys())
        gam_index=params_instance._fields.index(intemp_keys[0])   
        phi[gam_index] = theta[gam_index] *5000

        tempkeys=intemp_keys[1:]
        for i in range(len(tempkeys)):
            index=params_instance._fields.index(tempkeys[i])                                            
            phi[index] = theta[index] *3999  + 1

    if (proftype == 2):
                   
        alpha1_index=params_instance._fields.index('alpha1')
        alpha2_index=params_instance._fields.index('alpha2')
        logP1_index=params_instance._fields.index('logP1')
        logP3_index=params_instance._fields.index('logP3')
        T3_index=params_instance._fields.index('T3')
                                                                  
        # a1
        phi[alpha1_index] = 0.25 + (theta[alpha1_index]*0.25)
        # a2
        phi[alpha2_index] = 0.1 + (theta[alpha2_index] * 0.1)
        #P1
        phi[logP1_index] = (theta[logP1_index]* \
                             (np.log10(press[-1]) - np.log10(press[0]))) + np.log10(press[0])
        #P3
        #P3 must be greater than P1
        phi[logP3_index] = (theta[logP3_index] * \
                             (np.log10(press[-1]) - phi[logP1_index])) + phi[logP1_index]
        #T3
        phi[T3_index] = (theta[T3_index] * 3000.) + 1500.0


    elif (proftype == 3):
                                               
                                               
        alpha1_index=params_instance._fields.index('alpha1')
        alpha2_index=params_instance._fields.index('alpha2')
        logP1_index=params_instance._fields.index('logP1')
        logP2_index=params_instance._fields.index('logP2')
        logP3_index=params_instance._fields.index('logP3')
        T3_index=params_instance._fields.index('T3')
                                               

        # a1
        phi[alpha1_index] = 0.25 + (theta[alpha1_index]*0.25)
        # a2
        phi[alpha2_index] = 0.1 * (theta[alpha2_index] * 0.1)
        #P3 in press[0]--press[-1]
        phi[logP3_index] = (theta[logP3_index] * (np.log10(press[-1]) - np.log10(press[0]))) + np.log10(press[0])                                        
                                               
        # press[0]<P1<P3
        phi[logP1_index] = (theta[logP1_index]* (phi[logP3_index] - np.log10(press[0]))) + np.log10(press[0])
        ## press[0]<P2<P3
        phi[logP2_index] = (theta[logP2_index]* (phi[logP3_index] - np.log10(press[0]))) + np.log10(press[0])
       
        #T3
        phi[T3_index] = (theta[T3_index] * 3000.) + 1500.

    
    elif (proftype == 7):
                                                   
        Tint_index=params_instance._fields.index('Tint')
        alpha_index=params_instance._fields.index('alpha')
        lndelta_index=params_instance._fields.index('lndelta')
        T1_index=params_instance._fields.index('T1')
        T2_index=params_instance._fields.index('T2')
        T3_index=params_instance._fields.index('T3')
                                                                   
                                               
       # Tint - prior following Molliere+2020
        phi[Tint_index] = 300 + (theta[Tint_index] * 2000)
        # alpha, between 1 and 2
        phi[alpha_index] = theta[alpha_index] + 1. 
        # lndlelta
        plen = np.log10(press[-1]) - np.log10(press[0])
        pmax=phi[alpha_index]*plen 
        p_diff=np.log(0.1)-phi[alpha_index]*np.log10(press[-1])
        phi[lndelta_index] = theta[lndelta_index]*pmax+p_diff                                           

        # T1
        phi[T1_index] = 10. + (theta[T1_index] * 4000)
        # T2
        phi[T2_index] = 10. + (theta[T2_index] * 4000)
        # T3
        phi[T3_index] = 10.+ (theta[T3_index] * 4000)

    elif (proftype == 77):
                                                   
        Tint_index=params_instance._fields.index('Tint')
        alpha_index=params_instance._fields.index('alpha')
        lndelta_index=params_instance._fields.index('lndelta')
        T1_index=params_instance._fields.index('T1')
        T2_index=params_instance._fields.index('T2')
        T3_index=params_instance._fields.index('T3')
                                                                   
                                               
       # Tint - prior following Molliere+2020
        phi[Tint_index] = 300 + (theta[Tint_index] * 2000)
        # alpha, between 1 and 2
        phi[alpha_index] = theta[alpha_index] + 1. 
        # lndlelta
        plen = np.log10(press[-1]) - np.log10(press[0])
        pmax=phi[alpha_index]*plen 
        p_diff=np.log(0.1)-phi[alpha_index]*np.log10(press[-1])
        phi[lndelta_index] = theta[lndelta_index]*pmax+p_diff                                           

        # T1
        phi[T1_index] = 10. + (theta[T1_index] * 4000)
        # T2
        phi[T2_index] = 10. + (theta[T2_index] * 4000)
        # T3
        phi[T3_index] = 10.+ (theta[T3_index] * 4000)

        intemp_keys = list(re_params.dictionary['pt']['params'].keys())
        gam_index=params_instance._fields.index(intemp_keys[0])   
        phi[gam_index] = theta[gam_index] *5000


    npatches = do_clouds.size
    # only really ready for 2 patches here
    if (npatches > 1):
        fcld_index=params_instance._fields.index('fcld')
        phi[fcld_index] = theta[fcld_index]

    if (cloudtype.size > cloudtype.shape[1]):
        nclouds = cloudtype.shape[1]
    else:
        nclouds = cloudtype.size

    if np.all(do_clouds!= 0):
        cloudlist=list(re_params.dictionary["cloud"]["patch 1"].keys())

        for cloudname in cloudlist:
            if cloudname=='grey cloud deck':
            # 'cloudnum': 99,'cloudtype':2,
                logp_gcd_index=params_instance._fields.index('logp_gcd')
                dp_gcd_index=params_instance._fields.index('dp_gcd')
                #cloud top
                phi[logp_gcd_index] = \
                    (theta[logp_gcd_index] *(np.log10(press[-1]) \
                                    - np.log10(press[0])))\
                                    + np.log10(press[0])
                # cloud height
                phi[dp_gcd_index] = theta[dp_gcd_index] * 7.
                        
            elif cloudname=='grey cloud slab':
            # 'cloudnum': 99,'cloudtype':1,
                tau_gcs_index=params_instance._fields.index('tau_gcs')
                logp_gcs_index=params_instance._fields.index('logp_gcs')
                dp_gcs_index=params_instance._fields.index('dp_gcs')
                # cloud tau
                phi[tau_gcs_index] = theta[tau_gcs_index]*100.
                #cloud base
                phi[logp_gcs_index] = \
                    (theta[logp_gcs_index] *(np.log10(press[-1]) \
                                        - np.log10(press[0]))) \
                                        + np.log10(press[0])
                # cloud height
                phi[dp_gcs_index] = theta[dp_gcs_index] *\
                    (phi[logp_gcs_index] - np.log10(press[0]))
                                    
        
            elif cloudname=='powerlaw cloud deck':
            # 'cloudnum': 89,'cloudtype':2,
                logp_pcd_index=params_instance._fields.index('logp_pcd')
                dp_pcd_index=params_instance._fields.index('dp_pcd')
                alpha_pcd_index=params_instance._fields.index('alpha_pcd') 
                #cloud top
                phi[logp_pcd_index] = \
                            (theta[logp_pcd_index] *(np.log10(press[-1]) \
                                            - np.log10(press[0]))) \
                                            + np.log10(press[0])
                # cloud height
                phi[dp_pcd_index] = theta[dp_pcd_index] * 7.
                # power law
                phi[alpha_pcd_index] = (theta[alpha_pcd_index] * 20.) - 10.


            elif 'Mie scattering cloud deck' in cloudname:
            #   'cloudnum': cloudnum,'cloudtype':2,

                cloudspecies=cloudname.split('--')[1].strip()
                logp_pcd_index=params_instance._fields.index('logp_mcd_%s'%cloudspecies)
                dp_pcd_index=params_instance._fields.index('dp_mcd_%s'%cloudspecies)
                #cloud base
                phi[logp_pcd_index] = \
                    (theta[logp_pcd_index] *(np.log10(press[-1]) \
                                    - np.log10(press[0])))\
                                    + np.log10(press[0])

                
                # cloud height
                phi[dp_pcd_index] = theta[dp_pcd_index] * 7.

                particle_dis=re_params.dictionary["cloud"]["patch 1"]["particle_dis"]
                if  particle_dis=="hansan": 
                    hansan_a_mcd_index=params_instance._fields.index('hansan_a_mcd_%s'%cloudspecies)
                    hansan_b_mcd_index=params_instance._fields.index('hansan_b_mcd_%s'%cloudspecies)                                               
                    # particle effective radius
                    phi[hansan_a_mcd_index] = (theta[hansan_a_mcd_index] * 6.) - 3.
                    # particle spread
                    phi[hansan_b_mcd_index] = theta[hansan_b_mcd_index]
                elif particle_dis=="log_normal":
                    mu_mcd_index=params_instance._fields.index('mu_mcd_%s'%cloudspecies)
                    sigma_mcd_index=params_instance._fields.index('sigma_mcd_%s'%cloudspecies)
                    # particle effective radius
                    phi[mu_mcd_index] = (theta[mu_mcd_index] * 6.) - 3.
                    # particle spread
                    phi[mu_mcd_index] = theta[mu_mcd_index]                                                 

            elif cloudname=='power law cloud slab':
                    # 'cloudnum': 89, 'cloudtype':1,
                    tau_pcs_index=params_instance._fields.index('tau_pcs')
                    logp_pcs_index=params_instance._fields.index('logp_pcs')
                    dp_pcs_index=params_instance._fields.index('dp_pcs')
                    alpha_pcs_index=params_instance._fields.index('alpha_pcs')

                    # cloud tau
                    phi[tau_pcs_index] = theta[tau_pcs_index]*100.
                    #cloud base
                    phi[logp_pcs_index] = \
                        (theta[logp_pcs_index]*\
                            (np.log10(press[-1]) - np.log10(press[0]))) \
                            + np.log10(press[0])
                    # cloud height
                    phi[dp_pcs_index] = \
                        theta[dp_pcs_index] * (phi[logp_pcs_index] \
                                            - np.log10(press[0]))
                    # power law
                    phi[alpha_pcs_index] = (theta[alpha_pcs_index] * 20.) - 10.


            elif 'Mie scattering cloud slab' in cloudname:

                cloudspecies=cloudname.split('--')[1].strip()

                tau_mcs_index=params_instance._fields.index('tau_mcs_%s'%cloudspecies)
                logp_mcs_index=params_instance._fields.index('logp_mcs_%s'%cloudspecies)
                dp_mcs_index=params_instance._fields.index('dp_mcs_%s'%cloudspecies)
                # cloud tau
                phi[tau_mcs_index] = theta[tau_mcs_index]*100.
                #cloud base
                phi[logp_mcs_index] = \
                    (theta[logp_mcs_index] *(np.log10(press[-1]) \
                                        - np.log10(press[0]))) \
                                        + np.log10(press[0])
                # cloud height
                phi[dp_mcs_index] = theta[dp_mcs_index] * \
                    (phi[logp_mcs_index] - np.log10(press[0]))

                particle_dis=re_params.dictionary["cloud"]["patch 1"]["particle_dis"]
                                                                        
                if particle_dis=="hansan": 
                    hansan_a_mcs_index=params_instance._fields.index('hansan_a_mcs_%s'%cloudspecies)
                    hansan_b_mcs_index=params_instance._fields.index('hansan_b_mcs_%s'%cloudspecies)                                               
                    # particle effective radius
                    phi[hansan_a_mcs_index] = (theta[hansan_a_mcs_index] * 6.) - 3.
                    # particle spread
                    phi[hansan_b_mcs_index] = theta[hansan_b_mcs_index]
                elif particle_dis=="log_normal":
                    mu_mcs_index=params_instance._fields.index('mu_mcs_%s'%cloudspecies)
                    sigma_mcs_index=params_instance._fields.index('sigma_mcs_%s'%cloudspecies)
                    # particle effective radius
                    phi[mu_mcs_index] = (theta[mu_mcs_index] * 6.) - 3.
                    # particle spread
                    phi[mu_mcs_index] = theta[mu_mcs_index]


    return phi
                                                                     

def lnlike(theta,re_params):


    all_params,all_params_values =utils.get_all_parametres(re_params.dictionary) 
    params_master = namedtuple('params',all_params)
    params_instance = params_master(*theta)
    args_instance=settings.runargs

# Unpack all necessary parameters into local variables
    (   gases_myP,
        chemeq,
        dist,
        cloudtype,
        do_clouds,
        gaslist,
        gasnames,
        gasmass,
        cloudnum,
        inlinetemps,
        coarsePress,
        press,
        inwavenum,
        linelist,
        cia,
        ciatemps,
        use_disort,
        fwhm,
        obspec,
        proftype,
        do_fudge,
        prof,
        do_bff,
        bff_raw,
        ceTgrid,
        metscale,
        coscale
    ) = (
        args_instance.gases_myP,
        args_instance.chemeq,
        args_instance.dist,
        args_instance.cloudtype,
        args_instance.do_clouds,
        args_instance.gaslist,
        args_instance.gasnames,
        args_instance.gasmass,
        args_instance.cloudnum,
        args_instance.inlinetemps,
        args_instance.coarsePress,
        args_instance.press,
        args_instance.inwavenum,
        args_instance.linelist,
        args_instance.cia,
        args_instance.ciatemps,
        args_instance.use_disort,
        args_instance.fwhm,
        args_instance.obspec,
        args_instance.proftype,
        args_instance.do_fudge,
        args_instance.prof,
        args_instance.do_bff,
        args_instance.bff_raw,
        args_instance.ceTgrid,
        args_instance.metscale,
        args_instance.coscale
    )

    # get the spectrum
    # for MCMC runs we don't want diagnostics
    gnostics = 0
    shiftspec, photspec,tauspec,cfunc = modelspec(theta,re_params,args_instance,gnostics)


    # Get the scaling factors for the spectra. What is the FWHM? Negative number: preset combination of instruments
    if (fwhm < 0.0):
        if (fwhm == -1 or fwhm == -3 or fwhm == -4):
            scale1 =  params_instance.scale1 
            scale2 =  params_instance.scale2
            if (do_fudge == 1):
                logf = [params_instance.tolerance_parameter_1,params_instance.tolerance_parameter_2,params_instance.tolerance_parameter_3] #theta[ng+5:ng+8]
            else:
                # This is a place holder value so the code doesn't break
                logf = np.log10(0.1*(max(obspec[2,10::3]))**2)
        elif (fwhm == -2):
            scale1 = params_instance.r2d2
            if (do_fudge == 1):
                logf = [params_instance.tolerance_parameter_1,params_instance.tolerance_parameter_2]
            else:
                # This is a place holder value so the code doesn't break
                logf = np.log10(0.1*(max(obspec[2,10::3]))**2)
        elif (fwhm == -5):
            if (do_fudge == 1):
                logf = params_instance.tolerance_parameter_1
            else:
                # This is a place holder value so the code doesn't break
                logf = np.log10(0.1*(max(obspec[2,10::3]))**2)
        elif (fwhm == -6):
            if (do_fudge == 1):
                logf = params_instance.tolerance_parameter_1
            else:
                # This is a place holder value so the code doesn't break
                logf = np.log10(0.1*(max(obspec[2,10::3]))**2)
    else:
        if (do_fudge == 1):
            logf = params_instance.tolerance_parameter_1
        else:
            # This is a place holder value so the code doesn't break
            logf = np.log10(0.1*(max(obspec[2,10::3]))**2)

    modspec = np.array([shiftspec[0,::-1],shiftspec[1,::-1]])

    # If we've set a value for FWHM that we're using...
    if (fwhm > 0.00 and fwhm < 1.00):
        # this is a uniform FWHM in microns
        spec = conv_uniform_FWHM(obspec,modspec,fwhm)
        if (do_fudge == 1):
            s2=obspec[2,:]**2 + 10.**logf
        else:
            s2 = obspec[2,:]**2

        lnLik=-0.5*np.sum((((obspec[1,:] - spec[:])**2) / s2) + np.log(2.*np.pi*s2))       
    elif (fwhm > 10.00):
        # this is a uniform resolving power R.
        Res = fwhm
        spec = conv_uniform_R(obspec,modspec,Res)
        if (do_fudge == 1):
            s2=obspec[2,:]**2 + 10.**logf
        else:
            s2 = obspec[2,:]**2

        lnLik=-0.5*np.sum((((obspec[1,:] - spec[:])**2) / s2) + np.log(2.*np.pi*s2))
    elif (fwhm == 0.0):
        # Use convolution for Spex
        spec = prism_non_uniform(obspec,modspec,3.3)
        if (do_fudge == 1):
            s2=obspec[2,::3]**2 + 10.**logf
        else:
            s2 = obspec[2,::3]**2

        lnLik=-0.5*np.sum((((obspec[1,::3] - spec[::3])**2) / s2) + np.log(2.*np.pi*s2))
    elif (fwhm == 1.0):
        # Use convolution for JWST-NIRSpec PRISM
        spec = prism_non_uniform(obspec,modspec,2.2)
        if (do_fudge == 1):
            s2=obspec[2,::2]**2 + 10.**logf
        else:
            s2 = obspec[2,::2]**2

        lnLik=-0.5*np.sum((((obspec[1,::2] - spec[::2])**2) / s2) + np.log(2.*np.pi*s2))
    elif (fwhm == 2.0):
        # combo of JWST-NIRSpec PRISM + G395H grism
        # single scaling & single fudge factor
        spec = np.zeros_like(obspec[0,:])
        # first convolution for JWST-NIRSpec PRISM
        or1  = np.where(obspec[0,:] < 2.9)
        spec[or1] = prism_non_uniform(obspec[:,or1],modspec,2.2)
        # now 1st grism bit
        dL = 0.0015
        or2  = np.where(np.logical_and(obspec[0,:] > 2.9,obspec[0,:] < 3.69))
        spec[or2] =  conv_uniform_FWHM(obspec[:,or2],modspec,dL)
        # a bit more prism
        or3 = np.where(np.logical_and(obspec[0,:] > 3.69,obspec[0,:] < 3.785))
        spec[or3] = prism_non_uniform(obspec[:,or3],modspec,2.2)
        # 2nd bit of grism
        or4 = np.where(np.logical_and(obspec[0,:] > 3.785,obspec[0,:] < 5.14))
        spec[or4] =  conv_uniform_FWHM(obspec[:,or4],modspec,dL)
        # the rest of prism
        or5 = np.where(obspec[0,:] > 5.14)
        spec[or5] = prism_non_uniform(obspec[:,or5],modspec,2.2)
        if (do_fudge == 1):
            s2=obspec[2,:]**2 + 10.**logf
        else:
            s2 = obspec[2,:]**2

        lnLik=-0.5*np.sum((((obspec[1,:] - spec[:])**2) / s2) + np.log(2.*np.pi*s2))

        
            
    elif (fwhm < 0.0):
        lnLik = 0.0
        # This is for multi-instrument cases
        # -1: spex + akari + IRS
        # -2: spex + IRS
        # -3: spex + Lband + IRS
        if (fwhm == -1):

            # Spex
            or1  = np.where(obspec[0,:] < 2.5)
            spec1 = prism_non_uniform(obspec[:,or1],modspec,3.3)

            # AKARI IRC
            # dispersion constant across order 0.0097um
            # R = 100 at 3.6um for emission lines
            # dL ~constant at 3.6 / 120
            dL = 0.03
            or2 = np.where(np.logical_and(obspec[0,:] > 2.5,obspec[0,:] < 5.0))
            spec2 = scale1 * conv_uniform_FWHM(obspec[:,or2],modspec,dL)

            # Spitzer IRS
            # R roughly constant within orders, and orders both appear to
            # have R ~ 100
            R = 100.0
            or3 = np.where(obspec[0,:] > 5.0)
            spec3 = scale2 * conv_uniform_R(obspec[:,or3],modspec,R)

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
            or1  = np.where(obspec[0,:] < 2.5)
            spec1 = prism_non_uniform(obspec[:,or1],modspec,3.3)

            # Spitzer IRS
            # R roughly constant within orders, and orders both appear to
            # have R ~ 100
            R = 100.0
            or3 = np.where(obspec[0,:] > 5.0)
            spec3 = scale1 * conv_uniform_R(obspec[:,or3],modspec,R)

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
            or1  = np.where(obspec[0,:] < 2.5)
            spec1 = prism_non_uniform(obspec[:,or1],modspec,3.3)

            # Mike Cushing supplied L band R = 425
            # dispersion constant across order 0.0097um
            # R = 425
            R = 425
            or2 = np.where(np.logical_and(obspec[0,:] > 2.5,obspec[0,:] < 5.0))
            spec2 = scale1 * conv_uniform_R(obspec[:,or2],modspec,R)

            # Spitzer IRS
            # R roughly constant within orders, and orders both appear to
            # have R ~ 100
            R = 100.0
            or3 = np.where(obspec[0,:] > 5.0)
            spec3 = scale2 * conv_uniform_R(obspec[:,or3],modspec,R)

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
            
        elif (fwhm == -4):
            # This is spex + GNIRS L band R = 600 + IRS
            # Spex
            or1  = np.where(obspec[0,:] < 2.5)
            spec1 = prism_non_uniform(obspec[:,or1],modspec,3.3)

            # Katelyn Allers spectrum of GNIRS R = 600
            # R = 600 @ 3.5um linearly increading across order
            # i.e. FWHM - 0.005833
            dL = 0.005833
            #dL = 0.0097

            or2 = np.where(np.logical_and(obspec[0,:] > 2.5,obspec[0,:] < 5.0))
            spec2 = scale1 * conv_uniform_FWHM(obspec[:,or2],modspec,dL)

            # Spitzer IRS
            # R roughly constant within orders, and orders both appear to
            # have R ~ 100
            R = 100.0
            #mr3 = np.where(modspec[0,:] > 5.0)
            or3 = np.where(obspec[0,:] > 5.0)
            spec3 = scale2 * conv_uniform_R(obspec[:,or3],modspec,R)

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

        elif (fwhm == -5):
            # This is JWST NIRSpec + MIRI MRS no scaling + 1 fudge
            join = np.array([0.,5.1,5.7,7.59,11.6,13.4,15.49,18.01,20.0])
            pix = np.array([2.2,1.9,2.0,2.2,2.4,3.1,3.0,3.3])

            # Now we just work through the Prism +MRS orders,
            # using mid point in overlap regions
            # divided into chunk based on fwhm of res element in pixels
            spec = np.zeros_like(obspec[0,:])
                                 
            for i in range(0,pix.size):
                bit = np.where(np.logical_and(obspec[0,:] > join[i],obspec[0,:] < join[i+1]))
                spec[bit] = prism_non_uniform(obspec[:,bit],modspec,pix[i])

         
            if (do_fudge == 1):
                s2 = obspec[2,:]**2 + 10.**logf
            else:
                s2 = obspec[2,:]**2

            lnLik=-0.5*np.sum((((obspec[1,:] - spec)**2) / s2) + np.log(2.*np.pi*s2))

        elif (fwhm == -6):
            # This is UKIRT orders 1 and 2 based on Geballe 1996 cuts 
            # Second Order                           
            # R ~ 780 x Lambda (linear increase across order)
            # Order 2 (0.95 - 1.40 um)
            # FWHM ~ 1.175/780 = 0.001506    
            dL1 = 0.001506
            or1  = np.where(obspec[0,:] < 1.585)

            spec1 = conv_uniform_FWHM(obspec[:,or1],modspec,dL1)

            # First Order                            
            # R ~ 390 x Lambda (linear increase across order)
            # Order 1 (1.30 - 5.50 um) 
            # FWHM ~ 3.4/390 = 0.008717
            dL2 = 0.008717
            or2 = np.where(obspec[0,:] > 1.585)

            spec2 = conv_uniform_FWHM(obspec[:,or2],modspec,dL2)

            if (do_fudge == 1):
                s1 = obspec[2,or1]**2 + 10.**logf
                s3 = obspec[2,or2]**2 + 10.**logf
            else:
                s1 = obspec[2,or1]**2
                s3 = obspec[2,or2]**2


            lnLik1=-0.5*np.sum((((obspec[1,or1[0][::7]] - spec1[::7])**2) / s1[0][::7]) + np.log(2.*np.pi*s1[0][::7]))
            lnLik3=-0.5*np.sum((((obspec[1,or2[0][::3]] - spec2[::3])**2) / s3[0][::3]) + np.log(2.*np.pi*s3[0][::3]))
            lnLik = lnLik1 + lnLik3

        elif (fwhm == -7):
            #This is CGS4 NIR + NIRC Lband * CGS4 Mband
            # CGS4 Second order R = 780xLambda
            dL1 = 0.001506
            or1 = np.where(obspec[0, :] < 1.585)
            spec1 = conv_uniform_FWHM(obspec[:, or1], modspec, dL1)

            # CGS4 First order R = 390xLambda
            dL2 = 0.008717
            or2 = np.where(np.logical_and(obspec[0, :] > 1.585, obspec[0, :] < 2.52))
            spec2 = conv_uniform_FWHM(obspec[:, or2], modspec, dL2)

            # Oppenheimer 1998 NIRC L band spectrum
            ###EDIT### Central wavelength @ 3.492 with FWHM=1.490 for lw band
            # Using R=164
            #dL3 = 0.0213
            R=164.0
            or3 = np.where(np.logical_and(obspec[0, :] > 2.52, obspec[0, :] < 4.15))
            spec3 = scale1 * conv_uniform_R(obspec[:, or3], modspec, R)

            # CGS4 M band
            # Order 1 using 1".2 slit, 75 line/mm grating, 150 mm focal length camera
            ###EDIT### R=400xLambda
            dL4 = 0.0085
            or4 = np.where(obspec[0, :] > 4.15)
            spec4 = scale2 * conv_uniform_FWHM(obspec[:, or4], modspec, dL4)

            if (do_fudge == 1):
                s1 = obspec[2, or1] ** 2 + 10. ** logf[0]
                s2 = obspec[2, or2] ** 2 + 10. ** logf[0]
                s3 = obspec[2, or3] ** 2 + 10. ** logf[1]
                s4 = obspec[2, or4] ** 2 + 10. ** logf[2]
            else:
                s1 = obspec[2, or1] ** 2
                s2 = obspec[2, or2] ** 2
                s3 = obspec[2, or3] ** 2
                s4 = obspec[2, or4] ** 2  
  
            lnLik1 = -0.5 * np.sum((((obspec[1, or1[0][::7]] - spec1[::7]) ** 2) / s1[0][::7]) + np.log(2. * np.pi * s1[0][::7]))
            lnLik2 = -0.5 * np.sum((((obspec[1, or2[0][::3]] - spec2[::3]) ** 2) / s2[0][::3]) + np.log(2. * np.pi * s2[0][::3]))
            lnLik3 = -0.5 * np.sum((((obspec[1, or3] - spec3) ** 2) / s3) + np.log(2. * np.pi * s3))
            lnLik4 = -0.5 * np.sum((((obspec[1, or4] - spec4) ** 2) / s4) + np.log(2. * np.pi * s4))

            lnLik = lnLik1 + lnLik2 + lnLik3 + lnLik4

    if np.isnan(lnLik):
        lnLik = -np.inf

    samplemode=re_params.samplemode.lower()
    if samplemode=="multinest":
        if proftype==1 or proftype==77:
            intemp_keys = list(re_params.dictionary['pt']['params'].keys())
            intemp = np.array([getattr(params_instance, key) for key in intemp_keys])
            T=intemp[1:]

            if proftype==77:
                T = np.empty([press.size])
                T[:] = -100.
                delta= np.exp(params_instance.lndelta)
                alpha=params_instance.alpha
                P1 = ((1/delta)**(1/alpha))
                # put prior on P1 to put it shallower than 100 bar   
                if  (1 < alpha  < 2. and P1 < 100. and P1 > press[0]
                    and params_instance.T1 > 0.0 and params_instance.T2 > 0.0 and params_instance.T3 > 0.0 and params_instance.Tint >0.0):
                    T = TPmod.set_prof(proftype,junkP,press,intemp) # allow inversion 

            # bits for smoothing in prior
            gam = params_instance.gamma

            if (gam>0 and  (min(T) > 1.0) and (max(T) < 6000.)):
                diff=np.roll(T,-1)-2.*T+np.roll(T,1)
                pp=len(T)
                logbeta = -5.0
                beta=10.**logbeta
                alpha=1.0
                x=gam
                invgamma=((beta**alpha)/math.gamma(alpha)) * (x**(-alpha-1)) * np.exp(-beta/x)
                prprob = (-0.5/gam)*np.sum(diff[1:-1]**2) - 0.5*pp*np.log(gam) + np.log(invgamma)

                lnLik+=prprob
        
    return lnLik


def modelspec(theta,re_params,args_instance,gnostics):


    all_params,all_params_values =utils.get_all_parametres(re_params.dictionary) 
    params_master = namedtuple('params',all_params)
    params_instance = params_master(*theta)

    # Unpack all necessary parameters into local variables
    (   gases_myP,
        chemeq,
        dist,
        cloudtype,
        do_clouds,
        gaslist,
        gasnames,
        gasmass,
        cloudnum,
        inlinetemps,
        coarsePress,
        press,
        inwavenum,
        linelist,
        cia,
        ciatemps,
        use_disort,
        fwhm,
        obspec,
        proftype,
        do_fudge,
        prof,
        do_bff,
        bff_raw,
        ceTgrid,
        metscale,
        coscale
    ) = (
        args_instance.gases_myP,
        args_instance.chemeq,
        args_instance.dist,
        args_instance.cloudtype,
        args_instance.do_clouds,
        args_instance.gaslist,
        args_instance.gasnames,
        args_instance.gasmass,
        args_instance.cloudnum,
        args_instance.inlinetemps,
        args_instance.coarsePress,
        args_instance.press,
        args_instance.inwavenum,
        args_instance.linelist,
        args_instance.cia,
        args_instance.ciatemps,
        args_instance.use_disort,
        args_instance.fwhm,
        args_instance.obspec,
        args_instance.proftype,
        args_instance.do_fudge,
        args_instance.prof,
        args_instance.do_bff,
        args_instance.bff_raw,
        args_instance.ceTgrid,
        args_instance.metscale,
        args_instance.coscale
    )
        
    nlayers = press.size
    if chemeq == 0:
        gas_keys = re_params.dictionary['gas'].keys()
        gas_keys=list(gas_keys)
        invmr=np.array([getattr(params_instance, key) for key in gas_keys])


    else:
        mh  = params_instance.mh
        co =  params_instance.co

        mfit = interp1d(metscale,gases_myP,axis=0)
        gases_myM = mfit(mh)
        cfit = interp1d(coscale,gases_myM,axis=0)
        invmr = cfit(co)


    if re_params.samplemode=='multinest':
        M= params_instance.M
        R= params_instance.R* 69911e3
        GM = (6.67E-11 * M*1.898e27)
        logg = np.log10(100.* GM / R**2.)
        D = (dist + (np.random.randn()*args_instance.dist_err)) * 3.086e16
        # D = dist * 3.086e16
        R2D2 = R**2. / D**2.

    
    if (fwhm < 0.0):
        if (fwhm == -1 or fwhm == -3 or fwhm == -4):
            if re_params.samplemode=='mcmc':
                r2d2 = [params_instance.r2d2,params_instance.scale1,params_instance.scale2]  #theta[ng+1:ng+4]
            dlam = params_instance.dlambda
            # if (do_fudge == 1):
            #     logf =[params_instance.tolerance_parameter_1,params_instance.tolerance_parameter_2,params_instance.tolerance_parameter_3] #theta[ng+5:ng+8]
            # else:
            #     # This is a place holder value so the code doesn't break
            #     logf = np.log10(0.1*(max(obspec[2,10::3]))**2)

        elif (fwhm == -2):
            if re_params.samplemode=='mcmc':
                r2d2 = [params_instance.r2d2,params_instance.scale1]  #theta[ng+1:ng+3]
            dlam = params_instance.dlambda
            # if (do_fudge == 1):
            #     logf =[params_instance.tolerance_parameter_1,params_instance.tolerance_parameter_2] # theta[ng+4:ng+6]
            # else:
            #     # This is a place holder value so the code doesn't break
            #     logf = np.log10(0.1*(max(obspec[2,10::3]))**2)

        elif (fwhm == -5):
            if re_params.samplemode=='mcmc':
                r2d2 = params_instance.r2d2
            dlam = params_instance.dlambda
            # if (do_fudge == 1):
            #     logf = params_instance.tolerance_parameter_1
            # else:
            #     # This is a place holder value so the code doesn't break
            #     logf = np.log10(0.1*(max(obspec[2,10::3]))**2)

        elif (fwhm == -6):
            if re_params.samplemode=='mcmc':
                r2d2 = params_instance.r2d2
            dlam = params_instance.dlambda
            # if (do_fudge == 1):
            #     logf = params_instance.tolerance_parameter_1
            # else:
            #     # This is a place holder value so the code doesn't break                                                                      
            #     logf = np.log10(0.1*(max(obspec[2,10::3]))**2)

    else:
        if re_params.samplemode=='mcmc':
            r2d2 = params_instance.r2d2
        dlam = params_instance.dlambda
        # if (do_fudge == 1):
        #     logf = params_instance.tolerance_parameter_1
        # else:
        #     # This is a place holder value so the code doesn't break
        #     logf = np.log10(0.1*(max(obspec[2,10::3]))**2)


        
    npatches = do_clouds.size
    if (npatches > 1):
        prat =  params_instance.fcld
        pcover = np.array([prat,(1.-prat)])

    else:
        pcover = 1.0
        
    # use correct unpack method depending on situation

    if ((npatches > 1) and np.all(do_clouds != 0)):
        cloudparams = cloud_dic.unpack_patchy(re_params,params_instance,cloudtype,cloudnum,do_clouds)
    else:
        cloudparams = cloud_dic.unpack_default(re_params,params_instance,cloudtype,cloudnum,do_clouds)

    ndim = len(theta)

    intemp_keys = list(re_params.dictionary['pt']['params'].keys())
    intemp = np.array([getattr(params_instance, key) for key in intemp_keys])

    if (proftype == 1):
        gam = params_instance.gamma
        intemp=intemp[1:]

    elif (proftype == 77):
        gam = params_instance.gamma
        intemp=intemp[1:]
    elif (proftype == 2 or proftype ==3 or proftype ==7):
        intemp=intemp
    elif (proftype == 9):
        intemp = prof

    else:
        raise ValueError("not valid profile type %proftype" % (char, string))

    # set the profile
    temp = TPmod.set_prof(proftype,coarsePress,press,intemp)

    ngas = len(gaslist)
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

        tmpvmr = np.empty(ngas,dtype='d')
        if (gaslist[len(gaslist)-1] == 'na' and gaslist[len(gaslist)-2] == 'k'):
            tmpvmr[0:(ngas-2)] = invmr[0:-1]
            tmpvmr[ngas-2] = np.log10(10.**invmr[-1] / (alkratio+1.)) # K
            tmpvmr[ngas-1] = np.log10(10.**invmr[-1] * (alkratio / (alkratio+1.))) # Na
        elif (gaslist[len(gaslist)-1] == 'cs'):
            #f values are ratios between Na and (K+Cs) and K and Cs respectively
            f1 = 1.348
            f2 = 8912.5
            tmpvmr[0:(ngas-3)] = invmr[0:-1]
            tmpvmr[ngas-1] = np.log10(10.**invmr[-1] / ((f1+1)*(f2+1))) # Cs
            tmpvmr[ngas-2] = np.log10(10.**invmr[-1] * (f1 /(f1+1)) ) # Na
            tmpvmr[ngas-3] = np.log10(10.**invmr[-1] - 10.**tmpvmr[ngas-2] - 10.**tmpvmr[ngas-1]) #K
        else:
            tmpvmr[0:ngas] = invmr[0:ngas]
            
        for i in range(0,ngas):
            logVMR[i,:] = tmpvmr[i]

        # # set H- in the high atmosphere    
        # if (gaslist[gasnum.size-1] == 'hmins'):
        #     logVMR[ngas-1,0:P_hmins_index]=tmpvmr[ngas-1]
        #     logVMR[ngas-1,P_hmins_index:]=-100

        # # set non-uniform gas profile 
        
        gastype_values = [info['gastype'] for key, info in re_params.dictionary['gas'].items() if 'gastype' in info]
            
        for i in range(len(gastype_values)):
            if  gastype_values[i]=="N":
                P_gas= getattr(params_instance, "p_ref_%s"%gas_keys[i])
                gas_alpha= getattr(params_instance, "alpha_%s"%gas_keys[i])
                t_gas= getattr(params_instance, gas_keys[i])
                gas_profile=gas_nonuniform.non_uniform_gas(press,P_gas,t_gas,gas_alpha)
                logVMR[i]=gas_profile

            elif gastype_values[i]=="H":

                p_gas= getattr(params_instance, "p_ref_%s"%gas_keys[i])
                P_gas = 10**p_gas

                if np.size(np.where((press>=P_gas))[0])==0:   #may return null array
                    p_gas_index=np.size(press)-1
                else:
                    p_gas_index=np.where((press>=P_gas))[0][0] 

                logVMR[i,0:p_gas_index]=tmpvmr[i]
                logVMR[i,p_gas_index:]=-100





    # now need to translate cloudparams in to cloud profile even
    # if do_clouds is zero..

    cloudprof,cloudrad,cloudsig = cloud_dic.atlas(do_clouds,cloudnum,cloudtype,cloudparams,press)
    cloudprof = np.asfortranarray(cloudprof,dtype = 'float64')
    cloudrad = np.asfortranarray(cloudrad,dtype = 'float64')
    cloudsig = np.asfortranarray(cloudsig,dtype = 'float64')
    pcover = np.asfortranarray(pcover,dtype = 'float32')
    cloudnum = np.asfortranarray(cloudnum,dtype='i')
    do_clouds = np.asfortranarray(do_clouds,dtype = 'i')

    # get r2d2 sorted for multi-instruments
    if re_params.samplemode=='mcmc':
        logg = params_instance.logg
        if (fwhm < 0.0):
            if (fwhm == -1 or fwhm == -3 or fwhm == -4):
                R2D2 = r2d2[0]
                scale1 = r2d2[1]
                scale2 = r2d2[2]
            elif (fwhm == -2):
                R2D2 = r2d2[0]
                scale1 = r2d2[1]
            elif (fwhm == -5):
                R2D2 = r2d2
            elif (fwhm == -6):
                R2D2 = r2d2
        else:
            R2D2 = r2d2


    # Now get the BFF stuff sorted
    if (chemeq == 0 and do_bff == 1):
        for gas in range(0,3):
            for i in range(0,nlayers):
                tfit = InterpolatedUnivariateSpline(ceTgrid,bff_raw[:,i,gas],k=1)
                bff[gas,i] = tfit(temp[i])

        if (gaslist[len(gaslist)-1] == 'hmins'):
            bff[2,0:p_gas_index] = -50

    bff = np.asfortranarray(bff, dtype='float64')
    press = np.asfortranarray(press,dtype='float32')
    temp = np.asfortranarray(temp,dtype='float64')
    logVMR = np.asfortranarray(logVMR,dtype='float64')
    
    
    # print(invmr)
    # print(cloudprof,cloudrad,cloudsig)
    # print(temp)


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
    outspec,tmpclphotspec,tmpophotspec,cf = forwardmodel.marv(temp,logg,R2D2,gasnames,gasmass,logVMR,pcover,do_clouds,cloudnum,cloudrad,cloudsig,cloudprof,inlinetemps,press,inwavenum,linelist,cia,ciatemps,use_disort,clphot,ophot,make_cf,do_bff,bff)
        

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

    # shiftspec[1,:][np.isnan(shiftspec[1,:])] = 1e-20

    return shiftspec, cloud_phot_press,other_phot_press,cfunc




