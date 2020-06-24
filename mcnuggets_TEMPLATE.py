#!/usr/bin/env python

""" McNuggets: the post-processing tool for brewster"""
from __future__ import print_function

from builtins import str
from builtins import range
import numpy as np
import scipy as sp
import testkit
import ciamod
import TPmod
import nugbits_TEMPLATE as nb
import settings
import os
import gc
import sys
import pickle
import time
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from mpi4py import MPI

__author__ = "Ben Burningham"
__copyright__ = "Copyright 2015 - Ben Burningham"
__credits__ = ["Ben Burningham","The EMCEE DOCS"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ben Burningham"
__email__ = "burninghamster@gmail.com"
__status__ = "Development"

def split(container, count):
#    """
#    Simple function splitting a container into equal length chunks.
#    Order is not preserved but this is potentially an advantage depending on
#    the use case.
#    """
    return [container[_i::count] for _i in range(count)]


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

# set up parallel bits
    # Start by getting the chain, then set up the run arguments, then loop for Teff etc
settings.init()

COMM = MPI.COMM_WORLD


# only the first rank has to do all the crap...
runname "YOUR RUNNAME HERE"

# using the Burrow's alkalis (via Mike)?
malk = 0

#Are we testing?
testrun = 0

# length of test??
testlen = 30

# OK finish?
fin = 1

# what's the error on the distance?
sigDist = 0.1

# What's the error on the photometry used to flux calibrate the data?
sigPhot = 0.02

# Where are the pickles stored?
outdir = "/beegfs/car/bb/"

# which opacity set did we use?
xlist = "gaslistR10K.dat"

# Where are the cross sections?
# give the full path
xpath = "/beegfs/car/bb/Linelists/"

# that's all the input.. .off we go...
rawsamples, probs,ndim = nb.get_endchain(outdir+runname,fin)



slen = rawsamples.shape[0]

if (testrun == 1):
    settings.samples = rawsamples[slen-testlen:,:]
    if (COMM.rank == 0):
        print('testing with '+str(testlen)+' samples')
        print('full run would be '+str(slen))
else:
    settings.samples = rawsamples
    
slen = settings.samples.shape[0]
    
#print(((str(slen), ' samples for post production')))
    
    
# get the run arguments
gases_myP,chemeq,dist,cloudtype, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge,prof,do_bff,bff_raw,ceTgrid,metscale,coscale = nb.getargs(outdir+runname)

# If we're corrcting the distance (e.g for Gaia update),
# put it here, or or comment out 
# dist = 11.56

# But we want to a much wider wavelength range...
w1 = 0.7
w2 = 20.0

# So, we'll use gasnum from runargs^^ to replace the opacity arrays from the
# retrieval with new ones covering our new wavelength range. 

with open(xlist, 'r') as infile:
    lines=infile.readlines()[1:]
gaslist = []
for i in gasnum:
    for j in range(0,len(lines)-1):
        index = lines[j].split()[0]
        if (int(index) == i):
            gaslist.append(lines[j].split()[1])

# Now we'll get the opacity files into an array
inlinetemps,inwavenum,linelist,gasnum,nwave = testkit.get_opacities(gaslist,w1,w2,press,xpath,xlist,malk)

# Get the cia bits
tmpcia, ciatemps = ciamod.read_cia("CIA_DS_aug_2015.dat",inwavenum)
cia = np.asfortranarray(np.empty((4,ciatemps.size,nwave)),dtype='float32')
cia[:,:,:] = tmpcia[:,:,:nwave] 
ciatemps = np.asfortranarray(ciatemps, dtype='float32')


# now write runargs to global variable in settings:

settings.runargs = gases_myP,chemeq,dist,cloudtype, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge,prof,do_bff,bff_raw,ceTgrid,metscale,coscale




# Collect whatever has to be done in a list. Here we'll just collect a list of
# numbers. Only the first rank has to do this.
if COMM.rank == 0:
    jobs = list(range(slen))
    t1 = time.process_time()

    # Split into however many cores are available.
    jobs = split(jobs, COMM.size)
else:
    jobs = None

# Scatter jobs across cores.
jobs = COMM.scatter(jobs, root=0)

# Now each rank just does its jobs and collects everything in a results list.
# Make sure to not use super big objects in there as they will be pickled to be
# exchanged over MPI.
results = []

for job in jobs:
    results.append(nb.teffRM(settings.samples[job,0:ndim],sigDist,sigPhot))

# Gather results on rank 0.
results = MPI.COMM_WORLD.gather(results, root=0)

if COMM.rank == 0:
    # Flatten list of lists.
    results = [_i for tmp in results for _i in tmp]

    print("writing results to samplus")
    runtime = time.process_time() - t1 
    print('time taken = '+str(runtime))
    samplus = np.array(results)

    save_object(samplus,outdir+runname+'_postprod.pk1')

