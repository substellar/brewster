#!/usr/bin/env python

"""This is Brewster: the golden retriever of smelly atmospheres"""

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
__version__ = "0.1"
__maintainer__ = "Ben Burningham"
__email__ = "burninghamster@gmail.com"
__status__ = "Development"

# Now we set up the MPI bits
pool=MPIPool(loadbalance=True)
if not pool.is_master():
    pool.wait()
    sys.exit()

# This module set up the model arguments the drop these into
# theta(state vector) or runargs

# First get data and parameters for object

# Give the run name
runname = "2M2224_FeEnst2c2P_AlkTrim"

# get the observed spectrum
obspec = np.asfortranarray(np.loadtxt("2M2224_11_15_trim.dat",dtype='d',unpack='true'))

# Now the wavelength range
w1 = 1.1
w2 = 15.

# FWHM of data in microns(WE DON'T USE THIS FOR SPEX DATA. SET TO 0.0)
fwhm = -1

# DISTANCE (in parsecs)
dist = 11.35

# How many patches & clouds do we want??
# Must be at least 1 of each, but can turn off cloud below
npatches = 2
nclouds = 2

# set up array for setting patchy cloud answers
do_clouds = np.zeros([npatches],dtype='i')

# Which patchdes are cloudy
do_clouds[:] = 1

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
cloudnum[:,0] = 5
cloudtype[:,0] = 1

cloudnum[:,1] = 2
cloudtype[:,1] = 2

# second patch turn off top cloud
cloudnum[1,0] = 5
cloudtype[1,0] = 0

# Are we assuming chemical equilibrium, or similarly precomputed gas abundances?
# Or are we retrieving VMRs (0)
chemeq = 0

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
xpath = "/car-data/bb/Linelists/"

# now the cross sections

# Now list the gases.
# If Na is after K, at the end of the list, alkalis will be tied
# together at Asplund solar ratio. See Line at al (2015)
# Else if K is after Na, they'll be separate

gaslist = ['h2o','co','ch4','tio','vo','crh','feh','k','na']

ngas = len(gaslist)

# some switches for alternative cross sections
# Use Mike's Alkalis?
malk = 0
# Use Mike's CH4?
mch4 = 1

# now set up the EMCEE stuff
# How many dimensions???  Count them up in the p0 declaration. Carefully
ndim  = 31

# How many walkers we running?
nwalkers = ndim * 16

# How many burn runs do we want before reset?
nburn = 3000

# How many iterations are we running?
niter = 30000

# Is this a test or restart?
runtest = 0

# Are we writing the arguments to a pickle?
# Set = 1 for write and exit (no run); = 2 for write and continue
# option 2 may cause a memory issue and crash a production run
make_arg_pickle = 0

# Where is the output going?
outdir = "/car-data/bb/"

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


# If we want fresh guess set to 0, total inherit the previous set 1
# inherit plus randomise the VMRs. 2. See below to enter this filename
fresh = 0
p0 = np.empty([nwalkers,ndim])
if (fresh == 0):
    p0[:,0] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 3.5 # H2O
    p0[:,1] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 4.0 # CO
    p0[:,2] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 4.0 # CH4
    p0[:,3] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # TiO
    p0[:,4] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # VO
    p0[:,5] = (1.0*np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # CrH
    p0[:,6] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # FeH
    p0[:,7] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 5.5 # Na+K
    p0[:,8] = np.random.rand(nwalkers).reshape(nwalkers) + 4.0
    p0[:,9] =  1.0e-19 + (np.random.rand(nwalkers).reshape(nwalkers) * 5.e-21)
    p0[:,10] =  1.0 + (np.random.randn(nwalkers).reshape(nwalkers) * 0.1)
    p0[:,11] =  1.0 + (np.random.randn(nwalkers).reshape(nwalkers) * 0.1)    
    p0[:,12] = np.random.randn(nwalkers).reshape(nwalkers) * 0.001
    p0[:,13] = np.log10((np.random.rand(nwalkers).reshape(nwalkers) * (max(obspec[2,:]**2)*(0.1 - 0.01))) + (0.01*min(obspec[2,10::3]**2)))
    p0[:,14] = np.log10((np.random.rand(nwalkers).reshape(nwalkers) * (max(obspec[2,:]**2)*(0.1 - 0.01))) + (0.01*min(obspec[2,10::3]**2)))
    p0[:,15] = np.log10((np.random.rand(nwalkers).reshape(nwalkers) * (max(obspec[2,:]**2)*(0.1 - 0.01))) + (0.01*min(obspec[2,10::3]**2)))
    # some cloud bits now. two clouds, thin first, then deck, both power
    p0[:,16] = np.random.rand(nwalkers).reshape(nwalkers) # covering fraction
    p0[:,17] = np.random.rand(nwalkers).reshape(nwalkers)    
    p0[:,18] = -2. + np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,19] = np.random.rand(nwalkers).reshape(nwalkers)    
    p0[:,20] = -1. + 0.1*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,21] = 0.2*np.random.rand(nwalkers).reshape(nwalkers)
    # 2nd deep cloud
    p0[:,22] = 0.5 + 0.1*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,23] = np.random.rand(nwalkers).reshape(nwalkers)       
    p0[:,24] = np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,25] = np.random.rand(nwalkers).reshape(nwalkers)
    # And now the T-P params
    p0[:,26] = 0.39 + 0.1*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,27] = 0.14 +0.05*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,28] = -1.2 + 0.2*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,29] = 2.25+ 0.2*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,30] = 4200. + (500.*  np.random.randn(nwalkers).reshape(nwalkers))
    for i in range (0,nwalkers):
        while True:
            Tcheck = TPmod.set_prof(proftype,coarsePress,press,p0[i,ndim-5:])
            if (min(Tcheck) > 1.0):
                break
            else:
                p0[i,ndim-5] = 0.39 + 0.01*np.random.randn()
                p0[i,ndim-4] = 0.14 + 0.01*np.random.randn()
                p0[i,ndim-3] = -1.2 + 0.2*np.random.randn()
                p0[i,ndim-2] = 2. + 0.2*np.random.randn()
                p0[i,ndim-1] = 4200. + (200.*  np.random.randn())

    
if (fresh != 0):
    fname=chaindump
    pic=pickle.load(open(fname,'rb'))
    p0=pic
    if (fresh == 2):
        for i in range(0,9):
            p0[:,i] = (np.random.rand(nwalkers).reshape(nwalkers)*0.5) + p0[:,i]




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


# Sort out the BFF opacity stuff and chemical equilibrium tables:
metscale,coscale,Tgrid,Pgrid,gasnames,abunds = pickle.load( open( "chem_eq_tables.pic", "rb" ) )
nlayers = press.shape[0]
nabpress = Pgrid.size
nabtemp = Tgrid.size
nabgas = abunds.shape[4]
nmet = metscale.size
nco = coscale.size



bff_raw = np.zeros([nabtemp,nlayers,3])
gases_myP = np.zeros([nmet,nco,nabtemp,nlayers,ngas+3])
gases = np.zeros([nmet,nco,nabtemp,nabpress,ngas+3])

if (chemeq == 0):
    # Just want the ion fractions for solar metallicity in this case
    ab_myP = np.empty([nabtemp,nlayers,nabgas])
    i1 = np.where(metscale == 0.0)
    i2 = np.where(coscale == 1.0)
    for gas in range (0,nabgas):
        for i in range (0,nabtemp):
            pfit = InterpolatedUnivariateSpline(Pgrid,np.log10(abunds[i1[0],i2[0],i,:,gas]),k=1)
            ab_myP[i,:,gas] = pfit(np.log10(press))
            
            bff_raw = np.zeros([nabtemp,nlayers,3])
            bff_raw[:,:,0] = ab_myP[:,:,0]
            bff_raw[:,:,1] = ab_myP[:,:,2]
            bff_raw[:,:,2] = ab_myP[:,:,4]

else:
    # In this case we need the rows for the gases we're doing and ion fractions
    gases[:,:,:,:,0] = abunds[:,:,:,:,0]
    gases[:,:,:,:,1] = abunds[:,:,:,:,2]
    gases[:,:,:,:,2] = abunds[:,:,:,:,4]
    nmatch = 0 
    for i in range(0,ngas):
        for j in range(0,nabgas):
            if (gasnames[j].lower() == gaslist[i].lower()):
                gases[:,:,:,:,i+3] = abunds[:,:,:,:,j]
                nmatch = nmatch + 1
    if (nmatch != ngas):
        print "you've requested a gas that isn't in the Vischer table. Please chaeck and try again."
        exit
    
    for i in range(0,nmet):
        for j in range(0,nco):
            for k in range(0,ngas+3):
                for l in range(0,nabtemp):
                    pfit = InterpolatedUnivariateSpline(Pgrid,np.log10(gases[i,j,l,:,k]),k=1)
                    gases_myP[i,j,l,:,k] = pfit(np.log10(press))
    

    
            
ceTgrid = Tgrid

runargs = gases_myP,chemeq,dist, cloudtype,do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge, prof,do_bff,bff_raw,ceTgrid,metscale,coscale

    

# Write the arguments to a pickle if needed
if (make_arg_pickle > 0):
    pickle.dump(runargs,open(outdir+runname+"_runargs.pic","wb"))
    if( make_arg_pickle == 1):
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
    if (k==10 or k==1000 or k==1500 or k==2000 or k==2500 or k==3000 or k==3500 or k==4000 or k==4500 or k==5000 or k==6000 or k==7000 or k==8000 or k==9000 or k==10000 or k==11000 or k==12000 or k==15000 or k==18000 or k==19000 or k==20000 or k==21000 or k==22000 or k==23000 or k==24000 or k==25000 or k==26000 or k==27000 or k==28000 or k==29000 or k == 30000 or k == 35000 or k == 40000 or k == 45000 or k == 50000 or k == 55000 or k == 60000 or k == 65000):
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


