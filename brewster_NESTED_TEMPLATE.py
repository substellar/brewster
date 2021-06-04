#!/usr/bin/env python

"""This is Brewster: the golden retriever of smelly atmospheres"""
from __future__ import print_function
from __future__ import division

from builtins import str
from builtins import range
import multiprocessing
import time
import numpy as np
import scipy as sp
import pymultinest as mn
import nestkit
import ciamod
import TPmod
import settings
import os
import gc
import sys
import pickle
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
import mpi4py



__author__ = "Ben Burningham"
__copyright__ = "Copyright 2015 - Ben Burningham"
__credits__ = ["Ben Burningham"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ben Burningham"
__email__ = "burninghamster@gmail.com"
__status__ = "Development"


# This module set up the model arguments the drop these into
# theta(state vector) or runargs

# First get data and parameters for object

# Give the run name
runname = "G570D_NestTest"

# get the observed spectrum
obspec = np.asfortranarray(np.loadtxt("G570D_2MHcalib.dat",dtype='d',unpack='true'))

# Now the wavelength range
w1 = 1.0
w2 = 2.5

# FWHM of data in microns(WE DON'T USE THIS FOR SPEX DATA. SET TO 0.0)
# -1 = spex + AKARI + spitzer
# -2 = spex + spitzer
fwhm = 0

# DISTANCE (in parsecs)
dist = 5.882
dist_err = 0.005

# How many patches & clouds do we want??
# Must be at least 1 of each, but can turn off cloud below
npatches = 1
nclouds = 1

# set up array for setting patchy cloud answers
do_clouds = np.zeros([npatches],dtype='i')

# Which patchdes are cloudy
do_clouds[:] = 0

# set up cloud detail arrays
cloudnum = np.zeros([npatches,nclouds],dtype='i')
cloudtype =np.zeros([npatches,nclouds],dtype='i')

# Now fill cloud details. What kind of clouds and shape are they?
# Cloud types
# 1:  slab cloud
# 2: deep thick cloud , we only see the top
# 3: slab with fixed thickness log dP = 0.005 (~1% height)
# 4: deep thick cloud with fixed height log dP = 0.005
# In both cases the cloud properties are density, rg, rsig for real clouds
# and dtau, w0, and power law for cloudnum = 89 or 99 for grey
cloudnum[:,0] = 1
cloudtype[:,0] = 1

#cloudnum[:,1] = 10
#cloudtype[:,1] = 1

# second patch turn off top cloud
#cloudnum[1,0] = 5
#cloudtype[1,0] = 1

# Are we assuming chemical equilibrium, or similarly precomputed gas abundances?
# Or are we retrieving VMRs (0)
chemeq = 0

# Are we doing H- bound-free, free-free continuum opacities?
# (Is the profile going above 3000K in the photosphere?)
do_bff = 0

# Set the profile type. If we're using a fixed one. Give the file name
proftype = 1
pfile = "t1700g1000f3.dat"


# set up pressure grids in bar cos its intuitive
logcoarsePress = np.linspace(-4.0, 2.4, 5)
logfinePress = np.arange(-4.0, 2.4, 0.1)
# forward model wants pressure in bar
#logcoarsePress = np.arange(-4.0, 3.0, 0.5)
coarsePress = pow(10,logcoarsePress)
press = pow(10,logfinePress)


# Where are the cross sections?
# give the full path
xpath = "/nobackup/bburning/Linelists/"
xlist = 'gaslistRox.dat'


# now the cross sections

# Now list the gases.
# If Na is after K, at the end of the list, alkalis will be tied
# together at Asplund solar ratio. See Line at al (2015)
# Else if K is after Na, they'll be separate

gaslist = ['h2o','ch4','co','co2','nh3','h2s','K','Na']

ngas = len(gaslist)

# some switches for alternatives... 
# Use Mike's (Burrows') alkali opacities?
malk = 1

# Is this a test or restart?
runtest = 0

# Are we writing the arguments to a pickle?
# Set = 1 for write and exit (no run); = 2 for write and continue
# option 2 may cause a memory issue and crash a production run
make_arg_pickle = 2

# Where is the output going?
outdir = "/nobackup/bburning/"

# Names for the final output files:
finalout = runname+".pk1"

# Are we using DISORT for radiative transfer?
# (HINT: Not in this century)
use_disort = 0 

# use the fudge factor?
do_fudge = 1


prof = np.full(5,100.)
if (proftype == 9):
    modP,modT = np.loadtxt(pfile,skiprows=1,usecols=(1,2),unpack=True)
    tfit = InterpolatedUnivariateSpline(np.log10(modP),modT,k=1)
    prof = tfit(logcoarsePress)


# Now we'll get the opacity files into an array
inlinetemps,inwavenum,linelist,gasnum,nwave = nestkit.get_opacities(gaslist,w1,w2,press,xpath,xlist,malk)

# Get the cia bits
tmpcia, ciatemps = ciamod.read_cia("CIA_DS_aug_2015.dat",inwavenum)
cia = np.asfortranarray(np.empty((4,ciatemps.size,nwave)),dtype='float32')
cia[:,:,:] = tmpcia[:,:,:nwave] 
ciatemps = np.asfortranarray(ciatemps, dtype='float32')

# grab BFF and Chemical grids
bff_raw,ceTgrid,metscale,coscale,gases_myP = nestkit.sort_bff_and_CE(chemeq,"chem_eq_tables_P3K.pic",press,gaslist)


settings.init()
settings.runargs = gases_myP,chemeq,dist,dist_err,cloudtype,do_clouds,gasnum,gaslist,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge, prof,do_bff,bff_raw,ceTgrid,metscale,coscale




# Write the arguments to a pickle if needed
if (make_arg_pickle > 0):
    pickle.dump(settings.runargs,open(outdir+runname+"_runargs.pic","wb"))
    if( make_arg_pickle == 1):
        sys.exit()

    
# put it all together in the sampler..

n_params = nestkit.countdims(settings.runargs)

#print(n_params)

result = mn.solve(LogLikelihood=nestkit.lnlike, Prior=nestkit.priormap,n_dims=n_params,n_live_points=6000, outputfiles_basename=outdir+runname, verbose=True)

print()
print('evidence: %(logZ).1f +- %(logZerr).1f' % result)
print()
print('parameter values:')
for col in zip(result['samples'].transpose()):
	print('%.3f +- %.3f' % (col.mean(), col.std()))

