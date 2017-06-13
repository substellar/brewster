#!/usr/bin/env python

"""This is Brewster: the golden retriever of substellar atmospheres"""

import multiprocessing
import time
import numpy as np
import scipy as sp
import emcee
import testkit
import ciamod
import TPmod
import os
import gc
import sys
import pickle
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from emcee.utils import MPIPool





__author__ = "Ben Burningham"
__copyright__ = "Copyright 2015 - Ben Burningham"
__credits__ = ["Ben Burningham"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Ben Burningham"
__email__ = "burninghamster@gmail.com"
__status__ = "beta"


# This module set up the model arguments the drop these into
# theta(state vector) or runargs

# First get data and parameters for object

# Give the run name
runname = "D1425_pow"

# get the observed spectrum
# this spectrum should be trimmed to the w1 w2 coverage you want to use 
obspec = np.asfortranarray(np.loadtxt("D1425_2MassJcalib.dat",dtype='d',unpack='true'))

# Now the wavelength range
# set this here. Make sure you're doing what you think you're doing...
w1 = 0.8
w2 = 2.5

# FWHM of data in microns(WE DON'T USE THIS FOR SPEX DATA. SET TO 0.0)
fwhm = 0.00

# DISTANCE (in parsecs)
dist = 11.58

# How many patches & clouds do we want??
# Must be at least 1 of each, but can turn off cloud below
npatches = 1
nclouds = 1

# set up array for setting patchy cloud answers
do_clouds = np.zeros([npatches],dtype='i')

# Which patches are cloudy
do_clouds[:] = 1

# set up cloud detail arrays
cloudnum = np.zeros([npatches,nclouds],dtype='i')
cloudtype =np.zeros([npatches,nclouds],dtype='i')

# Now fill cloud details. What kind of clouds and shape are they?
# Cloud types
# 1:  slab cloud
# 2: deep thick cloud , we only see the top
# 3: slab with fixed thickness log dP = 0.005 (i.e. fits within a single layer)
# 4: deep thick cloud with fixed height (i.e. appears in a single layer)
# In both cases the cloud properties are density, rg, rsig for real clouds
# and dtau, w0, and power law for cloudnum = 89 or 99 for grey
cloudnum[:,0] = 89
cloudtype[:,0] = 3

#cloudnum[:,1] = 99
#cloudtype[:,1] = 4



# Are we doing H- bound-free, free-free continuum opacities?
# (Is the profile going above 3000K in the photosphere?)
do_bff = 0

# Set the profile type. If we're using a fixed one. Give the file name
proftype = 2
pfile = "t1700g1000f3.dat"


# set up pressure grids in bar cos its intuitive
logcoarsePress = np.arange(-4.0, 2.5, 0.53)
logfinePress = np.arange(-4.0, 2.4, 0.1)
# forward model wants pressure in bar
#logcoarsePress = np.arange(-4.0, 3.0, 0.5)
coarsePress = pow(10,logcoarsePress)
press = pow(10,logfinePress)


# Where are the cross sections?
# give the full path
xpath = "/nobackup/bburning/Linelists/"

# now the cross sections

# Now list the gases.
# If Na is after K, at the end of the list, alkalis will be tied
# together at Asplund solar ratio. See Line at al (2015)
# Else if K is after Na, they'll be separate

gaslist = ['h2o','co','tio','vo','crh','feh','na','k']

ngas = len(gaslist)

# some switches for alternative cross sections
# Use Mike's Alkalis?
malk = 0
# Use Mike's CH4?
mch4 = 0

# now set up the EMCEE stuff
# How many dimensions???  Count them up in the p0 declaration. Carefully
ndim  = 21

# How many walkers we running?
nwalkers = ndim * 16

# How many burn runs do we want before reset?
nburn = 10000

# How many iterations are we running?
niter = 30000

# Is this a test?
runtest = 0


# Where is the output going?
outdir = "/nobackup/bburning/"

# Names for the final output files:

# full final sampler with likelihoods, chain, bells and whistles
finalout = runname+".pk1"


# periodic dumps/snapshots
# just the chain
chaindump = runname+"_last_nc.pic"
# The whole thing w/ probs
picdump = runname+"_snapshot.pic"

# Names for status file runtimes
statfile = "status_ball"+runname+".txt"
rfile = "runtimes_"+runname+".dat"

# Are we using DISORT for radiative transfer?
# (HINT: Not in this century)
use_disort = 0 

# use the fudge factor?
do_fudge = 1

# get first guess at scale factor
r2d2 = (7e7 / (dist*3.086e+16))**2.

# If we want fresh guess set to 0, total inherit the previous set 1
# inherit plus randomise the VMRs. 2. See below to enter this filename
fresh = 0
p0 = np.empty([nwalkers,ndim])
if (fresh == 0):
    p0[:,0] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 3.5 # H2O
    p0[:,1] = (np.random.randn(nwalkers).reshape(nwalkers)) - 3.5 # CO
    p0[:,2] = (np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # tio
    p0[:,3] = (np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # VO
    p0[:,4] = (np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # crh
    p0[:,5] = (np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # FeH
    p0[:,6] = (np.random.randn(nwalkers).reshape(nwalkers)) - 5.5 # Na
    p0[:,7] = (np.random.randn(nwalkers).reshape(nwalkers)) - 5.5 # K   
    p0[:,8] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) + 4.5
    p0[:,9] = (0.1 * np.random.randn(nwalkers).reshape(nwalkers) * r2d2) + r2d2
    p0[:,10] = np.random.randn(nwalkers).reshape(nwalkers) * 0.001
    p0[:,11] = np.log10((np.random.rand(nwalkers).reshape(nwalkers) * (max(obspec[2,:]**2)*(0.1 - 0.01))) + (0.01*min(obspec[2,10::3]**2)))
    #Cloud bits
    p0[:,12] = np.random.rand(nwalkers).reshape(nwalkers)
    p0[:,13] = -2. + 1.5*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,14] = np.random.randn(nwalkers).reshape(nwalkers)    
    p0[:,15] = np.random.randn(nwalkers).reshape(nwalkers)
    # And now the T-P params
    p0[:,16] = 0.39 + 0.1*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,17] = 0.14 +0.05*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,18] = -1.5 + np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,19] = 3.0 + np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,20] = 4500. + (500.*  np.random.randn(nwalkers).reshape(nwalkers))

    for i in range (0,nwalkers):
        while True:
            Tcheck = TPmod.set_prof(proftype,coarsePress,press,p0[i,16:])
            if (min(Tcheck) > 1.0):
                break
            else:
                p0[i,16] = 0.39 + 0.01*np.random.randn()
                p0[i,17] = 0.14 + 0.01*np.random.randn()
                p0[i,18] = -1.2 + 0.2*np.random.randn()
                p0[i,19] = 2. + 0.2*np.random.randn()
                p0[i,20] = 4000. + (200.*  np.random.randn())
                
    
if (fresh != 0):
    fname=runname+"_snapshot.pic"
    pic=pickle.load(open(fname,'rb'))
    p0=pic





prof = np.full(13,100.)
if (proftype == 9):
    modP,modT = np.loadtxt(pfile,skiprows=1,usecols=(1,2),unpack=True)
    tfit = InterpolatedUnivariateSpline(np.log10(modP),modT,k=1)
    prof = tfit(logcoarsePress)



# Now we'll get the opacity files into an array
totgas = 24
gasdata = []
with open('gaslist.dat') as fa:
     for line_aa in fa.readlines()[1:totgas+1]:
        line_aa = line_aa.strip()
        gasdata.append(line_aa.split())
    
    
list1 = []    
for i in range(0,ngas):
    for j in range(0,totgas):
            if (gasdata[j][1].lower() == gaslist[i].lower()):
                list1.append(gasdata[j])

if (malk == 1):
    for i in range (0,ngas):    
        list1[i] = [w.replace('K_xsecs.pic', 'K_Mike_xsecs.pic') for w in list1[i]]
        list1[i] = [w.replace('Na_xsecs.pic', 'Na_Mike_xsecs.pic') for w in list1[i]]

if (mch4 ==1):
    for i in range (0,ngas):    
        list1[i] = [w.replace('CH4_xsecs.pic', 'CH4_Mike_xsecs.pic') for w in list1[i]]
    

lists = [xpath+i[3] for i in list1[0:ngas]]
gasnum = np.asfortranarray(np.array([i[0] for i in list1[0:ngas]],dtype='i'))


# get the basic framework from water list
rawwavenum, inpress, inlinetemps, inlinelist = pickle.load( open(xpath+'/H2O_xsecs.pic', "rb" ) )

wn1 = 10000./w2
wn2 = 10000. / w1
inwavenum = np.asfortranarray(rawwavenum[np.where(np.logical_not(np.logical_or(rawwavenum[:] > wn2, rawwavenum[:] < wn1)))],dtype='float64')
ntemps = inlinetemps.size
npress= press.size
nwave = inwavenum.size
r1 = np.amin(np.where(np.logical_not(np.logical_or(rawwavenum[:] > wn2, rawwavenum[:] < wn1))))
r2 = np.amax(np.where(np.logical_not(np.logical_or(rawwavenum[:] > wn2, rawwavenum[:] < wn1))))

# Here we are interpolating the linelist onto our fine pressure scale.
# pickles have linelist as 4th entry....
linelist = (np.ones([ngas,npress,ntemps,nwave],order='F')).astype('float64', order='F')
for gas in range (0,ngas):
    inlinelist= pickle.load( open(lists[gas], "rb" ) )[3]
    for i in range (0,ntemps):
        for j in range (r1,r2+1):
            pfit = interp1d(np.log10(inpress),np.log10(inlinelist[:,i,j]))
            linelist[gas,:,i,(j-r1)] = np.asfortranarray(pfit(np.log10(press)))

linelist[np.isnan(linelist)] = -50.0



# Get the cia bits
tmpcia, ciatemps = ciamod.read_cia("CIA_DS_aug_2015.dat",inwavenum)
cia = np.asfortranarray(np.empty((4,ciatemps.size,nwave)),dtype='float32')
cia[:,:,:] = tmpcia[:,:,:nwave] 
ciatemps = np.asfortranarray(ciatemps, dtype='float32')

# Sort out the BFF opacity stuff:

intab = np.loadtxt("2015_06_1060grid_feh_00_co_10.txt",skiprows=1)
test = np.array(intab)
test2 = test.reshape(60,18,36)
Pgrid = test2[20:21,:,1].reshape(18)
Tgrid =  test2[:,10:11,0].reshape(60)
abunds= test2[:,:,2:]
nlayers = press.shape[0]
nabpress = 18
nabtemp = 60
nabgas = 34
ab_myP = np.empty([nabtemp,nlayers,nabgas])
for gas in range (0,nabgas):
    for i in range (0,nabtemp):
            pfit = InterpolatedUnivariateSpline(Pgrid,np.log10(abunds[i,:,gas]),k=1)
            ab_myP[i,:,gas] = pfit(np.log10(press))
            
bff_raw = np.empty([nabtemp,nlayers,3])
bff_raw[:,:,0] = ab_myP[:,:,0]
bff_raw[:,:,1] = ab_myP[:,:,2]
bff_raw[:,:,2] = ab_myP[:,:,4]

bfTgrid = Tgrid

runargs = dist, cloudtype,do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge, prof,do_bff,bff_raw,bfTgrid

    
# Now we set up the MPI bits
pool=MPIPool(loadbalance=True)
if not pool.is_master():
	pool.wait()
	sys.exit()


# put it all together in the sampler..


sampler = emcee.EnsembleSampler(nwalkers, ndim, testkit.lnprob, args=(runargs),pool=pool)
#'''
# run the sampler
print "running the sampler"
clock = np.empty(80000)
k=0
times = open(rfile,"w")
times.close()
if (runtest == 0):
    pos,prob,state = sampler.run_mcmc(p0,nburn)
    sampler.reset()
    p0 = pos
for result in sampler.sample(p0,iterations=niter):
    clock[k] = time.clock()
    if (k > 1):
        tcycle = clock[k] - clock[k-1]
        times = open(rfile,"a")
        times.write("*****TIME FOR CYCLE*****")
        times.write(str(tcycle))
        times.close()
    k=k+1
    position = result[0]
    f = open(statfile, "w")
    f.write("****Iteration*****")
    f.write(str(k))
    f.write("****Reduced Chi2*****")
    f.write(str(result[1] * -2.0/(obspec.shape[1] / 3.0)))
    f.write("****Accept Fraction*****")
    f.write(str(sampler.acceptance_fraction))
    f.write("*****Values****")
    f.write(str(result[0]))
    f.close()
    if (k==10 or k==1000 or k==1500 or k==2000 or k==2500 or k==3000 or k==3500 or k==4000 or k==4500 or k==5000 or k==6000 or k==7000 or k==8000 or k==9000 or k==10000 or k==11000 or k==12000 or k==15000 or k==18000 or k==21000 or k==25000) or k == 30000 or k == 35000 or k == 40000 or k == 45000 or k == 50000 or k == 55000 or k == 60000 or k == 65000:
        chain=sampler.chain
	lnprob=sampler.lnprobability
	output=[chain,lnprob]
	pickle.dump(output,open(outdir+picdump,"wb"))
	pickle.dump(chain[:,k-1,:], open(chaindump,'wb'))


# get rid of problematic bit of sampler object
del sampler.__dict__['pool']

def save_object(obj, filename):
    with open(filename, "wb") as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

pool.close()

save_object(sampler,outdir+finalout)


