from __future__ import print_function
import pickle
import numpy as np
import emcee
import os
from rotBroadInt import rot_int_cmj as rotBroad


def get_endchain(runname,fin,results_path='./'):
    if (fin == 1):
        pic = results_path+runname+".pk1"
        sampler = pickle_load(pic)
        nwalkers = sampler.chain.shape[0]
        niter = sampler.chain.shape[1]
        ndim = sampler.chain.shape[2]
        flatprobs = sampler.lnprobability[:,:].reshape((-1))
        max_like = flatprobs[np.argmax(flatprobs)]
        print("maximum likelihood = ", max_like)
        flatendchain = sampler.chain[:,niter-2000:,:].reshape((-1,ndim))
        if (emcee.__version__ == '3.0rc2'):
            flatendprobs = sampler.lnprobability[niter-2000:,:].reshape((-1))
        else:
            flatendprobs = sampler.lnprobability[:, niter-2000:].reshape((-1))
        theta_max_end = flatendchain[np.argmax(flatendprobs)]
        max_end_like = np.amax(flatendprobs)
        print("maximum likelihood in final 2K iterations= ", max_end_like)
        print("Mean autocorrelation time: {0:.3f} steps"
              .format(np.mean(sampler.get_autocorr_time(discard=0,c=10,quiet=True))))

    elif(fin ==0):
        pic = results_path+runname+"_snapshot.pic"
        chain,probs = pickle_load(pic)
        nwalkers = chain.shape[0]
        ntot = chain.shape[1]
        ndim = chain.shape[2]
        niter = int(np.count_nonzero(chain) / (nwalkers*ndim))
        flatprobs = probs[:,:].reshape((-1))
        max_like = flatprobs[np.argmax(probs)]
        print("Unfinished symphony. Number of successful iterations = ", niter)
        print("maximum likelihood = ", max_like)
        flatendchain = chain[:,(niter-2000):niter,:].reshape((-1,ndim))
        if (emcee.__version__ == '3.0rc2'):
            flatendprobs = probs[niter-2000:,:].reshape((-1))
        else:
            flatendprobs = probs[:, niter-2000:].reshape((-1))
        theta_max_end = flatendchain[np.argmax(flatendprobs)]
        max_end_like = np.amax(flatendprobs)
        print("maximum likelihood in final 2K iterations= ", max_end_like)
    else:
        print("File extension not recognised")
        stop
        
    return flatendchain, flatendprobs,ndim


def proc_spec(shiftspec,theta,fwhm,chemeq,gaslist,obspec):
    import numpy as np
    import scipy as sp
    from bensconv import prism_non_uniform
    from bensconv import conv_uniform_R
    from bensconv import conv_uniform_FWHM

    if chemeq == 0:
        if (gaslist[len(gaslist)-1] == 'Na'):
            ng = len(gaslist) - 1
        elif (gaslist[len(gaslist)-1] == 'Cs'):
            ng = len(gaslist) -2
        else:
            ng = len(gaslist)
            invmr = theta[0:ng]
        
    else:
        ng = 2

    if (fwhm < 0.0):
        if (fwhm == -1 or fwhm == -3 or fwhm == -4 or fwhm == -7):
            scale1 = theta[ng+2]
            scale2 = theta[ng+3]
        elif (fwhm == -2 or fwhm == -8):
            scale1 = theta[ng+2]
        elif (fwhm == -6):
            scale1 = 1.0

    elif (fwhm == 3.0):
        vrad = theta[ng+2]
        vsini = theta[ng+3]
        
    modspec = np.array([shiftspec[0,::-1],shiftspec[1,::-1]])

    # If we've set a value for FWHM that we're using... 
    if (fwhm > 0.00 and fwhm < 1.00):
        # this is a uniform FWHM in microns
        
        outspec = conv_uniform_FWHM(obspec,modspec,fwhm)

    elif (fwhm > 10.00):
        # this is a uniform resolving power R.
        Res = fwhm
        outspec = conv_uniform_R(obspec,modspec,Res)

    elif (fwhm == 0.0):
        # Use Mike's convolution for Spex
        outspec = prism_non_uniform(obspec,modspec,3.3)
    elif (fwhm == 1.0):
        # Use convolution for JWST/NIRSpec Prism
        outspec = prism_non_uniform(obspec,modspec,2.2)
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
        outspec = spec

    elif (fwhm == 3.0):
        # JWST NIRSpec G395H data can resolve vsini, so we do it here
        rotspec = rotBroad(modspec[0],modspec[1],vsini)
        modspec[1,:] = rotspec

        # JWST NIRSpec G395H data with 2.2 pixels per resolving element
        # Use convolution for Spex
        outspec = prism_non_uniform(obspec,modspec,2.2)

        
    elif (fwhm < 0.0):
        # This is for multi-instrument cases
        # -1: spex + akari + IRS
        # -2: spex + IRS
        # -3: spex + Lband + IRS
        if (fwhm == -1):

            # Spex
            mr1 = np.where(modspec[0,:] < 2.5)
            or1  = np.where(obspec[0,:] < 2.5)
            spec1 = prism_non_uniform(obspec[:,or1],modspec,3.3)

            # AKARI IRC
            # dispersion constant across order 0.0097um
            # R = 100 at 3.6um for emission lines
            # dL ~constant at 3.6 / 120
            dL = 0.03
            mr2 = np.where(np.logical_and(modspec[0,:] > 2.5,modspec[0,:] < 5.0))
            or2 = np.where(np.logical_and(obspec[0,:] > 2.5,obspec[0,:] < 5.0))
            spec2 = scale1 * conv_uniform_FWHM(obspec[:,or2],modspec,dL)

            # Spitzer IRS
            # R roughly constant within orders, and orders both appear to
            # have R ~ 100
            R = 100.0
            mr3 = np.where(modspec[0,:] > 5.0)
            or3 = np.where(obspec[0,:] > 5.0)
            spec3 = scale2 * conv_uniform_R(obspec[:,or3],modspec,R)

            outspec = np.array(np.concatenate((spec1,spec2,spec3),axis=0))

        elif (fwhm == -2):
            # This is just spex + IRS
            # Spex
            mr1 = np.where(modspec[0,:] < 2.5)
            or1  = np.where(obspec[0,:] < 2.5)
            spec1 = prism_non_uniform(obspec[:,or1],modspec,3.3)

 
            # Spitzer IRS
            # R roughly constant within orders, and orders both appear to
            # have R ~ 100
            R = 100.0
            mr3 = np.where(modspec[0,:] > 5.0)
            or3 = np.where(obspec[0,:] > 5.0)
            spec3 = scale1 * conv_uniform_R(obspec[:,or3],modspec,R)

            outspec = np.array(np.concatenate((spec1,spec3),axis=0))
            
        elif (fwhm == -3):
            # This is spex + Mike Cushing's L band R = 425 + IRS
            # Spex
            mr1 = np.where(modspec[0,:] < 2.5)
            or1  = np.where(obspec[0,:] < 2.5)
            spec1 = prism_non_uniform(obspec[:,or1],modspec,3.3)

             # Mike Cushing supplied L band R = 425 
            # dispersion constant across order 0.0097um
            # R = 425
            R = 425
            mr2 = np.where(np.logical_and(modspec[0,:] > 2.5,modspec[0,:] < 5.0))
            or2 = np.where(np.logical_and(obspec[0,:] > 2.5,obspec[0,:] < 5.0))
            spec2 = scale1 * conv_uniform_R(obspec[:,or2],modspec,R)

            # Spitzer IRS
            # R roughly constant within orders, and orders both appear to
            # have R ~ 100
            R = 100.0
            mr3 = np.where(modspec[0,:] > 5.0)
            or3 = np.where(obspec[0,:] > 5.0)
            spec3 = scale2 * conv_uniform_R(obspec[:,or3],modspec,R)

            outspec =  np.array(np.concatenate((spec1,spec2,spec3),axis=0))

        elif (fwhm == -4):
            # This is spex  + GNIRS L band R = 600 + IRS 
            # Spex
            mr1 = np.where(modspec[0,:] < 2.5)
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
            mr3 = np.where(modspec[0,:] > 5.0)
            or3 = np.where(obspec[0,:] > 5.0)
            spec3 = scale2 * conv_uniform_R(obspec[:,or3],modspec,R)

            outspec =  np.array(np.concatenate((spec1,spec2,spec3),axis=0))

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

            outspec = spec
            

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

            outspec =  np.array(np.concatenate((spec1,spec2),axis=0))

        elif (fwhm == -7):
            #This is CGS4 NIR + NIRC Lband + CGS4 Mband
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
            # dL3 = 0.0213
            R = 164.0
            or3 = np.where(np.logical_and(obspec[0, :] > 2.52, obspec[0, :] < 4.15))
            spec3 = scale1 * conv_uniform_R(obspec[:, or3], modspec, R)

            # CGS4 M band
            # Order 1 using 1".2 slit, 75 line/mm grating, 150 mm focal length camera
            ###EDIT### R=400xLambda
            dL4 = 0.0085
            or4 = np.where(obspec[0, :] > 4.15)
            spec4 = scale2 * conv_uniform_FWHM(obspec[:, or4], modspec, dL4)

            outspec = np.array(np.concatenate((spec1, spec2, spec3, spec4), axis=0))

        elif (fwhm == -8):
            # This is NIRSpec + MIRI, no order shifts, just instrument
            # NIRSpec
            R = 2700
            or1  = np.where(obspec[0,:] < 5.0)
            spec1 = conv_uniform_R(obspec[:,or1],modspec,R)

            #MIRI MRS roughly constant R = 2700
            R = 2700
            or3 = np.where(obspec[0,:] > 5.0)
            spec3 = scale1 * conv_uniform_R(obspec[:,or3],modspec,R)

            outspec = np.array(np.concatenate((spec1,spec3),axis=0))
            

    return outspec
            

class MacOSFile(object):

    def __init__(self, f):
        self.f = f

    def __getattr__(self, item):
        return getattr(self.f, item)

    def read(self, n):
        # print("reading total_bytes=%s" % n, flush=True)
        if n >= (1 << 31):
            buffer = bytearray(n)
            idx = 0
            while idx < n:
                batch_size = min(n - idx, 1 << 31 - 1)
                # print("reading bytes [%s,%s)..." % (idx, idx + batch_size), end="", flush=True)
                buffer[idx:idx + batch_size] = self.f.read(batch_size)
                # print("done.", flush=True)
                idx += batch_size
            return buffer
        return self.f.read(n)

    def write(self, buffer):
        n = len(buffer)
        print("writing total_bytes=%s..." % n, flush=True)
        idx = 0
        while idx < n:
            batch_size = min(n - idx, 1 << 31 - 1)
            print("writing bytes [%s, %s)... " % (idx, idx + batch_size), end="", flush=True)
            self.f.write(buffer[idx:idx + batch_size])
            print("done.", flush=True)
            idx += batch_size


def pickle_dump(obj, file_path):
    with open(file_path, "wb") as f:
        return pickle.dump(obj, MacOSFile(f), protocol=pickle.HIGHEST_PROTOCOL)


def pickle_load(file_path):
    with open(file_path, "rb") as f:
        return pickle.load(MacOSFile(f))
