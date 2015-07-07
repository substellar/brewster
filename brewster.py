#!/usr/bin/env python

"""This is Brewster: the golden retriever of smelly atmospheres"""


import numpy as np
import scipy as sp
import emcee
import testkit
import cPickle as pickle

__author__ = "Ben Burningham"
__copyright__ = "Copyright 2015 - Ben Burningham"
__credits__ = ["Ben Burningham"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ben Burningham"
__email__ = "burninghamster@gmail.com"
__status__ = "Development"


# set up the model arguments the drop these into theta(state vector) or runargs

intemp = np.asfortranarray(np.loadtxt("16temps.dat",dtype='d'))
r2d2 = 1.
logg = 4.5
dlam = 0.
w1 = 1.0
w2 = 2.5
pcover = 1.0
do_clouds = 0
# cloudparams is structured array with 5 entries
# each one has a patch*cloud entries
cloudparams = np.ones(5)
# 5 entries in cloudparams for simple slab model are:
# 0) log10(number density)
# 1) top layer id (or pressure)
# 2) base ID (these are both in 61 layers)
# 3) rg
# 4) rsig
cloudparams[0,0,0] = -20.
cloudparams[0,0,1] = 10
cloudparams[0,0,0] = 12
cloudparams[0,0,0] = 1e-4
cloudparams[0,0,0] = 1e-5
# hardwired gas and cloud IDs
gasnum = np.asfortranarray(np.array([1,2,20],dtype='i'))
cloudnum = np.array([1],dtype='i')

# hardwired FWHM of data
fwhm = 0.005


# get the observed spectrum
obspec = np.asfortranarray(np.loadtxt("sim_spectrum.dat",dtype='d',unpack='true'))

runargs = w1,w2,intemp, pcover, cloudparams, r2d2, logg, dlam, do_clouds,gasnum,cloudnum,fwhm,obspec

# now set up the EMCEE stuff
ndim, nwalkers = 3, 100
p0 = -10.* np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=runargs)
# run the sampler
sampler.run_mcmc(p0, 500)

# get rid of problematic bit of sampler object
del sampler.__dict__['pool']

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

save_object(sampler,'retrieval_result.pk1')



