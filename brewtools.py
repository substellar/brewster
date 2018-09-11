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


