#!/usr/bin/env python

""" Module of bits to plug into Brewster """
import math
import gc
import numpy as np
import scipy as sp
import forwardmodel
import cloud
import TPmod
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.convolution import convolve, convolve_fft
from astropy.convolution import Gaussian1DKernel
from mikesconv import instrument_non_uniform
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



def lnlike(intemp, invmr, pcover, cloudtype, cloudparams, r2d2, logg, dlam, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,logf,proftype,do_fudge,do_bff,bff_raw,bfTgrid):

    # Hard code nlayers
    nlayers = press.shape[0]
    # set the profile
    temp = TPmod.set_prof(proftype,coarsePress,press,intemp)

    # get the ngas for forward model (ngas, not ng
    if (gasnum[gasnum.size-1] == 21):
        ngas = invmr.shape[0] + 1
    elif (gasnum[gasnum.size-1] == 23):
        ngas = invmr.shape[0] + 2
    else:
        ngas = invmr.shape[0]
    # now loop through gases and get VMR for model
    # check if its a fixed VMR or a profile
    # VMR is log10(VMR) !!!
    logVMR = np.empty((ngas,nlayers),dtype='d')
    alkratio = 16.2 #  from Asplund et al (2009)
    if invmr.size > invmr.shape[0]:
        # this case is a profile
        # now sort Na and K
        tmpvmr = np.empty((ngas,nlayers),dtype='d')
        if (gasnum[gasnum.size-1] == 21):
            tmpvmr[0:(ngas-2),:] = invmr[0:(ngas-2),:]
            tmpvmr[ngas-2,:] = np.log10(10.**invmr[ngas-2,:] / (alkratio+1.)) # K
            tmpvmr[ngas-1,:] = np.log10(10.**invmr[ngas-2,:] * (alkratio / (alkratio+1.))) # Na                                
        elif (gasnum[gasnum.size-1] == 23):
            #f values are ratios between Na and (K+Cs) and K and Cs respectively
            f1 = 1.348
            f2 = 8912.5
            tmpvmr[0:(ngas-3),:] = invmr[0:(ngas-3),:]
            tmpvmr[ngas-1,:] = np.log10(10.**invmr[ngas-3,:] / ((f1+1)*(f2+1))) # Cs
            tmpvmr[ngas-2,:] = np.log10(10.**invmr[ngas-3,:] * (f1 /(f1+1)) ) # Na
            tmpvmr[ngas-3,:] = np.log10(10.**invmr[ngas-3,:] - 10.**tmpvmr[ngas-2,:] - 10.**tmpvmr[ngas-1,:]) #K 
        else:
            tmpvmr[0:ngas,:] = invmr[0:ngas,:]
            
        for i in range(0,ngas):
            vfit = sp.interpolate.splrep(np.log10(coarsepress),tmpvmr[i,:],s=0)
            logVMR[i,:] = sp.interpolate.splev(np.log10(press),vfit,der=0)
    else:
        # This caseis fixed VMR
        # now sort Na and K
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

    # Now get the BFF stuff sorted
    bff = np.zeros([3,nlayers],dtype="float64") 
    if (do_bff == 1):
        for gas in range(0,3):
            for i in range(0,nlayers):
                tfit = InterpolatedUnivariateSpline(bfTgrid,bff_raw[:,i,gas],k=1) 
                bff[gas,i] = 10.**tfit(temp[i])

    bff = np.asfortranarray(bff, dtype='float64')
    press = np.asfortranarray(press,dtype='float32')
    temp = np.asfortranarray(temp,dtype='float64')
    logVMR = np.asfortranarray(logVMR,dtype='float64')
    # Set pspec and tspec as we don't need these in the emcee run
    tspec = 0
    pspec = 0
    
    # now we can call the forward model
    outspec,photspec,tauspec = forwardmodel.marv(temp,logg,r2d2,gasnum,logVMR,pcover,do_clouds,cloudnum,cloudrad,cloudsig,cloudprof,inlinetemps,press,inwavenum,linelist,cia,ciatemps,use_disort,pspec,tspec,do_bff,bff)
    # Trim to length where it is defined.
    nwave = inwavenum.size
    trimspec = np.zeros((2,nwave),dtype='d')
    trimspec[:,:] = outspec[:,:nwave]
    #print trimspec
    # now shift wavelen by delta_lambda
    shiftspec = np.empty_like(trimspec)
    shiftspec[0,:] =  trimspec[0,:] + dlam
    shiftspec[1,:] =  trimspec[1,:]
 
    # length and interval for later
    #wlen = shiftspec.shape[1]
    #wint =  shiftspec[0,0] - shiftspec[0,wlen-1]

    # convolve with instrumental profile
    # start by setting up kernel
    # First step is finding the array index length of the FWHM
    #disp = wint / wlen
    #gwidth = int((((fwhm / disp) // 2) * 2) +1)

    # needs to be odd
    # now get the kernel and convolve
    #gauss = Gaussian1DKernel(gwidth)
    #cspec = convolve(shiftspec[1,:],gauss,boundary='extend')
    #spec = np.array([shiftspec[0,::-1],cspec[::-1]])
    
    # rebin to observed dispersion
    #wfit = sp.interpolate.splrep(spec[0,:],spec[1,:],s=0)
    #modspec = sp.interpolate.splev(obspec[0,:],wfit,der=0)
    
    # Below is method for rebinning using conserve flux method
    #    oblen = obspec.shape[1]
    #    modspec = np.empty((2,oblen),dtype='d')
    #    modspec[1,:] =  rebinspec(spec[0,:], spec[1,:], obspec[0,:])
    # get log-likelihood
    # We've lifted this from Mike's code, below is original from emcee docs
    # Just taking every 3rd point to keep independence

    # Use Mike's convolution
    wno = 1e4 / shiftspec[0,:]
    modspec = instrument_non_uniform(obspec[0,:],wno,shiftspec[1,:])

    if (do_fudge == 1):
        s2=obspec[2,::3]**2 + 10.**logf
    else:
        s2 = obspec[2,::3]**2
    lnLik=-0.5*np.sum((((obspec[1,::3] - modspec[::3])**2) / s2) + np.log(2.*np.pi*s2))


    return lnLik
    #chi2 log likelihood--can modify this
    #invsigma2 = 1.0/((obspec[2,::3])**2 + modspec[1,::3]**2 * np.exp(2*lnf))
    #return -0.5*(np.sum((obspec[1,::3] - modspec[1,::3])**2 * invsigma2 - np.log(invsigma2)))
    
    
def lnprob(theta,dist,cloudtype, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge,prof,do_bff,bff_raw,bfTgrid):
    
    if (gasnum[gasnum.size-1] == 21):
        ng = gasnum.size - 1
    elif (gasnum[gasnum.size-1] == 23):
        ng = gasnum.size -2
    else:
        ng = gasnum.size
        
    invmr = theta[0:ng]
    logg = theta[ng]
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
    
    if ((npatches > 1) and np.all(do_clouds > 0)):
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


    
    # now check against the priors, if not beyond them, run the likelihood
    lp = lnprior(theta,obspec,dist,proftype,press,do_clouds,gasnum,cloudnum,cloudtype,do_fudge)
    if not np.isfinite(lp):
        return -np.inf
    # else run the likelihood
    lnlike_value = lnlike(intemp, invmr,pcover, cloudtype,cloudparams,r2d2, logg, dlam, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,logf,proftype,do_fudge,do_bff,bff_raw,bfTgrid)

    lnprb = lp+lnlike_value
    if np.isnan(lnprb):
        lnprb = -np.inf
    return lnprb


def lnprior(theta,obspec,dist,proftype,press,do_clouds,gasnum,cloudnum,cloudtype,do_fudge):
    # set up the priors here
    if (gasnum[gasnum.size-1] == 21):
        ng = gasnum.size - 1
    elif (gasnum[gasnum.size-1] == 23):
        ng = gasnum.size - 2
    else:
        ng = gasnum.size
    invmr = theta[0:ng]
    logg = theta[ng]
    r2d2 = theta[ng+1]
    dlam = theta[ng+2]
    if (do_fudge == 1):
        logf = theta[ng+3]
        pc = ng + 4
    elif (do_fudge == 0):
        logf = np.log10(0.1*(max(obspec[2,10::3]))**2)
        pc = ng + 3
    npatches = do_clouds.size
    if (npatches > 1):
        prat = theta[ng+4]
        pcover = np.array([prat,(1.-prat)])
        pc = pc + 1
    else:
        pcover = np.array([0.5,0.5])
        
    # use correct unpack method depending on situation
    
    if ((npatches > 1) and np.all(do_clouds > 0)):
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
    cloud_dens0 =  np.empty_like(cloud_tau0)
    rg = np.empty_like(cloud_tau0)
    rsig = np.empty_like(cloud_tau0)

    if (sum(do_clouds) >= 1):    
        for i in range(0,npatches):
            if (do_clouds[i] == 1):
                for j in range (0, nclouds):
                    if (cloudnum[i,j] == 99):
                        if (cloudtype[i,j] == 1):
                            cloud_tau0[i,j] = cloudparams[0,i,j]
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = cloudparams[2,i,j]
                            cloud_bot[i,j] = cloud_top[i,j] + cloud_height[i,j]
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = 0.0
                            cloud_dens0[i,j] = -20.0
                            rg[i,j] = 1.0
                            rsig[i,j] = 0.1
                        elif (cloudtype[i,j] == 2):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[press.size - 1])
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = cloudparams[2,i,j]
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = 0.0
                            cloud_dens0[i,j] = -20.0
                            rg[i,j] = 1.0
                            rsig[i,j] = 0.1
                        elif (cloudtype[i,j] == 3):
                            cloud_tau0[i,j] = cloudparams[0,i,j]
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = 0.005
                            cloud_bot[i,j] = cloud_top[i,j] + cloud_height[i,j]
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = 0.0
                            cloud_dens0[i,j] = -20.0
                            rg[i,j] = 1.0
                            rsig[i,j] = 0.1
                        elif (cloudtype[i,j] == 4):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[press.size - 1])
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = 0.005
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = 0.0
                            cloud_dens0[i,j] = -20.0
                            rg[i,j] = 1.0
                            rsig[i,j] = 0.1
                        elif (cloudtype[i,j] == 0):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[press.size-1])
                            cloud_top[i,j] = np.log10(press[0])
                            cloud_height[i,j] = 0.1
                            w0[i,j] = +0.5
                            taupow[i,j] = 0.0
                            cloud_dens0[i,j] = -20.
                            rg[i,j] =  1.0
                            rsig[i,j] = 0.1
                    elif (cloudnum[i,j] == 89):
                        if (cloudtype[i,j] == 1):
                            cloud_tau0[i,j] = cloudparams[0,i,j]
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = cloudparams[2,i,j]
                            cloud_bot[i,j] = cloud_top[i,j] + cloud_height[i,j]
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = cloudparams[4,i,j]
                            cloud_dens0[i,j] = -20.0
                            rg[i,j] = 1.0
                            rsig[i,j] = 0.1
                        elif (cloudtype[i,j] == 2):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[press.size - 1])
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = cloudparams[2,i,j]
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = cloudparams[4,i,j]
                            cloud_dens0[i,j] = -20.0
                            rg[i,j] = 1.0
                            rsig[i,j] = 0.1
                        elif (cloudtype[i,j] == 3):
                            cloud_tau0[i,j] = cloudparams[0,i,j]
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = 0.005
                            cloud_bot[i,j] = cloud_top[i,j] + cloud_height[i,j]
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = cloudparams[4,i,j]
                            cloud_dens0[i,j] = -20.0
                            rg[i,j] = 1.0
                            rsig[i,j] = 0.1
                        elif (cloudtype[i,j] == 4):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[press.size - 1])
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = 0.005
                            w0[i,j] = cloudparams[3,i,j]
                            taupow[i,j] = cloudparams[4,i,j]
                            cloud_dens0[i,j] = -20.0
                            rg[i,j] = 1.0
                            rsig[i,j] = 0.1
                        elif (cloudtype[i,j] == 0):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[press.size-1])
                            cloud_top[i,j] = np.log10(press[0])
                            cloud_height[i,j] = 0.1
                            w0[i,j] = +0.5
                            taupow[i,j] = 0.0
                            cloud_dens0[i,j] = -20.
                            rg[i,j] =  1.0
                            rsig[i,j] = 0.1
                    else:
                        if (cloudtype[i,j] == 1):
                            cloud_tau0[i,j] = 1.0
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = cloudparams[2,i,j]
                            cloud_bot[i,j] = cloud_top[i,j] + cloud_height[i,j]
                            w0[i,j] = 0.5
                            taupow[i,j] = 0.0
                            cloud_dens0[i,j] = cloudparams[0,i,j]
                            rg[i,j] = cloudparams[3,i,j]
                            rsig[i,j] = cloudparams[4,i,j]
                        elif (cloudtype[i,j] == 2):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[press.size-1])
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = cloudparams[2,i,j]
                            w0[i,j] = +0.5
                            taupow[i,j] =0.0
                            cloud_dens0[i,j] = cloudparams[0,i,j]
                            rg[i,j] =  cloudparams[3,i,j]
                            rsig[i,j] =  cloudparams[4,i,j]
                        elif (cloudtype[i,j] == 3):
                            cloud_tau0[i,j] = 1.0
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = 0.005
                            cloud_bot[i,j] = cloud_top[i,j] + cloud_height[i,j]
                            w0[i,j] = 0.5
                            taupow[i,j] = 0.0
                            cloud_dens0[i,j] = cloudparams[0,i,j]
                            rg[i,j] = cloudparams[3,i,j]
                            rsig[i,j] = cloudparams[4,i,j]
                        elif (cloudtype[i,j] == 4):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[press.size-1])
                            cloud_top[i,j] = cloudparams[1,i,j]
                            cloud_height[i,j] = 0.005
                            w0[i,j] = +0.5
                            taupow[i,j] =0.0
                            cloud_dens0[i,j] = cloudparams[0,i,j]
                            rg[i,j] =  cloudparams[3,i,j]
                            rsig[i,j] =  cloudparams[4,i,j]
                        elif (cloudtype[i,j] == 0):
                            cloud_tau0[i,j] = 1.0
                            cloud_bot[i,j] = np.log10(press[press.size-1])
                            cloud_top[i,j] = np.log10(press[0])
                            cloud_height[i,j] = 0.1
                            w0[i,j] = +0.5
                            taupow[i,j] = 0.0
                            cloud_dens0[i,j] = -20.
                            rg[i,j] =  1.0
                            rsig[i,j] = 0.1
    else:
        cloud_tau0[:,:] = 1.0
        cloud_bot[:,:] = np.log10(press[press.size-1])
        cloud_top[:,:] = np.log10(press[0])
        cloud_height[:,:] = 0.1
        w0[:,:] = +0.5
        taupow[:,:] = 0.0
        cloud_dens0[:,:] = -20.
        rg[:,:] =  1.0
        rsig[:,:] = 0.1
                

    if (proftype == 1):
        gam = theta[pc+nc]
        T = theta[pc+nc:]
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
            and  0.0 < logg < 6.0 
            and 1.0 < M < 80. 
            and  0. < r2d2 < 1.
            and  0.5 < Rj < 2.0
            and -0.01 < dlam < 0.01 
            and (min(T) > 1.0) and (max(T) < 6000.) 
            and (gam > 0.)
            and ((0.01*np.min(obspec[2,:]**2)) < 10.**logf
                 < (100.*np.max(obspec[2,:]**2)))
            and (np.all(cloud_tau0 >= 0.0))
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
            and np.all(-25.0 < cloud_dens0)
            and np.all(cloud_dens0 < -11.0)
            and np.all(rg > 0.01)
            and np.all(rsig > 0.001)):
                 
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
        junkP = np.ones([13])
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
            and  0.0 < logg < 6.0 
            and 1.0 < M < 80. 
            and  0. <= r2d2 < 1.
            and  0.5 < Rj < 2.0
            and -0.01 < dlam < 0.01 
            and ((0.01*np.min(obspec[2,:]**2)) < 10.**logf
                 < (100.*np.max(obspec[2,:]**2)))
            and (np.all(cloud_tau0 >= 0.0))
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
            and np.all(-25.0 < cloud_dens0)
            and np.all(cloud_dens0 < -11.0)
            and (np.all(rg > 0.01))
            and (np.all(rsig > 0.001))
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
            and all(pcover > 0.) and (np.sum(pcover) == 1.0)
            and  0.0 < logg < 6.0 
            and 1.0 < M < 80. 
            and  0. < r2d2 < 1.
            and  0.5 < Rj < 2.0
            and -0.01 < dlam < 0.01 
            and ((0.01*np.min(obspec[2,:]**2)) < 10.**logf
                 < (100.*np.max(obspec[2,:]**2)))
            and all(cloud_tau0 >= 0.0)
            and (cloud_top < cloud_bot <= np.log10(press[press.size-1])).all
            and (np.log10(press[0]) <= cloud_top < cloud_bot).all
            and all(0. < cloud_height < 7.0)
            and all(0.0 < w0 < 1.0)
            and (all(-10.0 < taupow < +10.0))
            and  all(-25.0 < cloud_dens0 < -11.0)
            and all(rg > 0.01)
            and all(rsig > 0.001)
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
        print "Rj = ", Rj
        print "M = ", M
        print "logg = ", logg
        print "R2D2 = ", r2d2
        print "dlam = ", dlam
        print "VMRs = ", invmr
        print "ng = ", ng
        print "sum VMRs = ", np.sum(10.**(invmr[0:ng]))
        print "clouddens = ", cloud_dens0
        print "cloud_top = ", cloud_top
        print "cloud_bot = ", cloud_bot
        print "rg = ", rg
        print "rsig = ", rsig
        print "cloud_tau0 = ", cloud_tau0
        print "w0 = ", w0
        print "pcover = ", pcover
        print "sum(pcover) = ", np.sum(pcover)
        if (all(invmr[0:ng] > -12.0) and all(invmr[0:ng] < 0.0) and (np.sum(10.**(invmr[0:ng])) < 1.0)
            and  all(pcover > 0.) and (np.sum(pcover) == 1.0)
            and  0.0 < logg < 6.0 
            and 1.0 < M < 100. 
            and  0. <= r2d2 < 1.
            and  0.5 < Rj < 2.0
            and -0.01 < dlam < 0.01 
            and ((0.01*np.min(obspec[2,:]**2)) < 10.**logf
                 < (100.*np.max(obspec[2,:]**2)))
            and (np.all(cloud_tau0 >= 0.0))
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
            and np.all(-25.0 < cloud_dens0)
            and np.all(cloud_dens0 < -11.0)
            and (np.all(rg > 0.01))
            and (np.all(rsig > 0.001))):
            return 0.0
        return -np.inf
