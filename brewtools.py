from __future__ import print_function
import pickle
import numpy as np
import os

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
        flatendprobs = sampler.lnprobability[niter-2000:,:].reshape((-1))
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
        flatendprobs = probs[(niter-2000):niter,:].reshape((-1))
        theta_max_end = flatendchain[np.argmax(flatendprobs)]
        max_end_like = np.amax(flatendprobs)
        print("maximum likelihood in final 2K iterations= ", max_end_like)
    else:
        print("File extension not recognised")
        stop
        
    return flatendchain, flatendprobs,ndim


def proc_spec(shiftspec,theta,fwhm,chemeq,gasnum,obspec):
    import numpy as np
    import scipy as sp
    from bensconv import spex_non_uniform
    from bensconv import conv_uniform_R
    from bensconv import conv_uniform_FWHM

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
        elif (fwhm == -2):
            scale1 = theta[ng+2]
 
    modspec = np.array([shiftspec[0,::-1],shiftspec[1,::-1]])

    # If we've set a value for FWHM that we're using... 
    if (fwhm > 0.00 and fwhm < 1.00):
        # this is a uniform FWHM in microns
        
        outspec = conv_uniform_FWHM(obspec,modspec,fwhm)
   

    elif (fwhm == 0.0):
        # Use Mike's convolution for Spex
        outspec = spex_non_uniform(obspec,modspec)
             
    elif (fwhm < 0.0):
        # This is for multi-instrument cases
        # -1: spex + akari + IRS
        # -2: spex + IRS
        # -3: spex + Lband + IRS
        if (fwhm == -1):

            # Spex
            mr1 = np.where(modspec[0,:] < 2.5)
            or1  = np.where(obspec[0,:] < 2.5)
            spec1 = spex_non_uniform(obspec[:,or1],modspec)

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
            spec1 = spex_non_uniform(obspec[:,or1],modspec)

 
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
            spec1 = spex_non_uniform(obspec[:,or1],modspec)

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
