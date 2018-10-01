def get_endchain(runname,fin):
    import pickle as pickle
    import numpy as np
    if (fin == 1):
        pic = runname+".pk1"
        with open(pic, 'rb') as input:
            sampler = pickle.load(input) 
        nwalkers = sampler.chain.shape[0]
        niter = sampler.chain.shape[1]
        ndim = sampler.chain.shape[2]
        flatprobs = sampler.lnprobability[:,:].reshape((-1))
        max_like = flatprobs[np.argmax(flatprobs)]
        print "maximum likelihood = ", max_like
        flatendchain = sampler.chain[:,niter-2000:,:].reshape((-1,ndim))
        flatendprobs = sampler.lnprobability[:,niter-2000:].reshape((-1))
        theta_max_end = flatendchain[np.argmax(flatendprobs)]
        max_end_like = np.amax(flatendprobs)
        print "maximum likelihood in final 2K iterations= ", max_end_like

    elif(fin ==0):
        pic = runname+"_snapshot.pic"
        with open(pic, 'rb') as input:
            chain,probs = pickle.load(input) 
        nwalkers = chain.shape[0]
        ntot = chain.shape[1]
        ndim = chain.shape[2]
        niter = np.count_nonzero(chain) / (nwalkers*ndim)
        flatprobs = probs[:,:].reshape((-1))
        max_like = flatprobs[np.argmax(probs)]
        print "Unfinished symphony. Number of successful iterations = ", niter
        print "maximum likelihood = ", max_like
        flatendchain = chain[:,(niter)-2000:niter].reshape((-1,ndim))
        flatendprobs = probs[:,(niter)-2000:niter].reshape((-1))
        theta_max_end = flatendchain[np.argmax(flatendprobs)]
        max_end_like = np.amax(flatendprobs)
        print "maximum likelihood in final 2K iterations= ", max_end_like
    else:
        print "File extension not recognised"
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
 

    # If we've set a value for FWHM that we're using... 
    if (fwhm > 0.00 and fwhm < 1.00):
        # this is a uniform FWHM in microns
        
        outspec = conv_uniform_FWHM(shiftspec,obspec,fwhm)
   

    elif (fwhm == 0.0):
        # Use Mike's convolution for Spex
        wno = 1e4 / shiftspec[0,:]
        outspec = spex_non_uniform(obspec[0,:],wno,shiftspec[1,:])
             
    elif (fwhm < 0.0):
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

            outspec = np.array(np.concatenate((spec1,spec2,spec3),axis=0))

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

            outspec = np.array(np.concatenate((spec1,spec3),axis=0))
            
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

            outspec =  np.array(np.concatenate((spec1,spec2,spec3),axis=0))

    return outspec
            

