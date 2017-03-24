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


# Bit to fix CPU affinity after numpy import

#if __name__ == '__main__':
#    pool_size = multiprocessing.cpu_count()
#    os.system('taskset -cp 0-%d %s' % (pool_size, os.getpid()))



# set up the model arguments the drop these into theta(state vector) or runargs

# set up pressure grids in bar cos its intuitive
logcoarsePress = np.arange(-4.0, 2.5, 0.53)
#logcoarsePress = np.arange(-4.0, 3.0, 0.5)
coarsePress = pow(10,logcoarsePress)
logfinePress = np.arange(-4.0, 2.4, 0.1)
finePress = pow(10,logfinePress)
# forward model wants pressure in bar
press = finePress


w1 = 0.8
w2 = 2.4

dist = 11.35
# hardwired FWHM of data in microns
fwhm = 0.005

npatches = 1
nclouds = 1


do_clouds = np.array([1],dtype='i')

do_bff = 0

# CURRENTLY ONLY COPE WITH ONE CLOUDY PATCH.
#SO MAKE ALL CLOUD PARAMETERS THE SAME FOR EASE OF PROCESSING 

cloudnum = np.zeros([npatches,nclouds],dtype='i')
cloudnum[:,:] = 89
cloudtype = np.asfortranarray(np.ones([npatches,nclouds]),dtype='i')
cloudtype[:,:] = 2

use_disort = 0 

# use the fudge factor?
do_fudge = 1

# Set the profile type
proftype = 2

prof = np.full(13,100.)
if (proftype == 9):
    modP,modT = np.loadtxt("t1700g1000f3.dat",skiprows=1,usecols=(1,2),unpack=True)
    tfit = InterpolatedUnivariateSpline(np.log10(modP),modT,k=1)
    prof = tfit(logcoarsePress)


# now the linelist
# Set up number of gases, and point at the lists. see gaslist.dat
gasnum = np.asfortranarray(np.array([1,4,11,21,20],dtype='i'))
ngas = gasnum.size
lists = ["/nobackup/bburning/Linelists/H2O_xsecs.pic","/nobackup/bburning/Linelists/co_xsecs.pic","/nobackup/bburning/Linelists/feh_xsecs.pic","/nobackup/bburning/Linelists/Na_xsecs.pic","/nobackup/bburning/Linelists/K_xsecs.pic"]
# get the basic framework from water list
rawwavenum, inpress, inlinetemps, inlinelist = pickle.load( open('/nobackup/bburning/Linelists/H2O_xsecs.pic', "rb" ) )

wn1 = 10000./w2
wn2 = 10000. / w1
inwavenum = np.asfortranarray(rawwavenum[np.where(np.logical_not(np.logical_or(rawwavenum[:] > wn2, rawwavenum[:] < wn1)))],dtype='float64')
ntemps = inlinetemps.size
npress= finePress.size
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
            linelist[gas,:,i,(j-r1)] = np.asfortranarray(pfit(np.log10(finePress)))

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
nlayers = finePress.shape[0]
nabpress = 18
nabtemp = 60
nabgas = 34
ab_myP = np.empty([nabtemp,nlayers,nabgas])
for gas in range (0,nabgas):
    for i in range (0,nabtemp):
            pfit = InterpolatedUnivariateSpline(Pgrid,np.log10(abunds[i,:,gas]),k=1)
            ab_myP[i,:,gas] = pfit(np.log10(finePress))
            
bff_raw = np.empty([nabtemp,nlayers,3])
bff_raw[:,:,0] = ab_myP[:,:,0]
bff_raw[:,:,1] = ab_myP[:,:,2]
bff_raw[:,:,2] = ab_myP[:,:,4]

bfTgrid = Tgrid



# get the observed spectrum
obspec = np.asfortranarray(np.loadtxt("2M2224_mkoJcalib_trim.dat",dtype='d',unpack='true'))


runargs = dist, cloudtype,do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge, prof,do_bff,bff_raw,bfTgrid

# now set up the EMCEE stuff

ndim  = 18 #((ngas-1) + 9 + 5)
nwalkers = ndim * 32
#int(((ndim * ndim) // 2) * 2)


# If we want fresh guess set to 0, total inherit the previous set 1, inherit plus randomise the VMRs. 2.
fresh = 0
p0 = np.empty([nwalkers,ndim])
if (fresh == 0):
    p0[:,0] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 3.5 # H2O
    p0[:,1] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 4.0 # CO
    p0[:,2] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # FeH
    p0[:,3] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 5.5 # Na
    p0[:,4] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 5.5 # K
    p0[:,5] = np.random.rand(nwalkers).reshape(nwalkers) + 4.2
    p0[:,6] =  1.0e-20 + np.random.rand(nwalkers).reshape(nwalkers) * 5.e-20
    p0[:,7] = np.random.randn(nwalkers).reshape(nwalkers) * 0.001
    p0[:,8] = np.log10((np.random.rand(nwalkers).reshape(nwalkers) * (max(obspec[2,:]**2)*(0.1 - 0.01))) + (0.01*min(obspec[2,10::3]**2)))
    # some cloud bits now. We're just doing grey cloud, tau so need pressure of top where plus cloud height (in dex), SSA, don't need GG
#    p0[:,12] = 0.5* np.random.rand(nwalkers).reshape(nwalkers)
    p0[:,9] = 1.0 + 0.1*np.random.randn(nwalkers).reshape(nwalkers) 
    p0[:,10] = np.random.rand(nwalkers).reshape(nwalkers) 
    p0[:,11] = np.random.rand(nwalkers).reshape(nwalkers)
    p0[:,12] = -1. * np.random.rand(nwalkers).reshape(nwalkers)
    # And now the T-P params
    p0[:,13] = 0.39 + 0.1*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,14] = 0.14 +0.05*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,15] = -1.2 + 0.2*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,16] = 2.25+ 0.2*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,17] = 4200. + (500.*  np.random.randn(nwalkers).reshape(nwalkers))

    for i in range (0,nwalkers):
        while True:
            Tcheck = TPmod.set_prof(proftype,coarsePress,press,p0[i,13:])
            if (min(Tcheck) > 1.0):
                break
            else:
                p0[i,13] = 0.39 + 0.01*np.random.randn()
                p0[i,14] = 0.14 + 0.01*np.random.randn()
                p0[i,15] = -1.2 + 0.2*np.random.randn()
                p0[i,16] = 2. + 0.2*np.random.randn()
                p0[i,17] = 4200. + (200.*  np.random.randn())

    
if (fresh != 0):
    fname='MCMC_last_pow.pic'
    pic=pickle.load(open(fname,'rb'))
    p0=pic
    if (fresh == 2):
        for i in range(0,9):
            p0[:,i] = (np.random.rand(nwalkers).reshape(nwalkers)*0.5) + p0[:,i]


    
# Now we set up the MPI bits
pool=MPIPool(loadbalance=True)
if not pool.is_master():
	pool.wait()
	sys.exit()



sampler = emcee.EnsembleSampler(nwalkers, ndim, testkit.lnprob, args=(runargs),pool=pool)
#'''
# run the sampler
print "running the sampler"
#sampler.run_mcmc(p0, 100)
clock = np.empty(60000)
k=0
times = open("runtimes_5g.dat","w")
times.close()
pos,prob,state = sampler.run_mcmc(p0,5000)
sampler.reset()
for result in sampler.sample(pos,iterations=30000):
    clock[k] = time.clock()
    if (k > 1):
        tcycle = clock[k] - clock[k-1]
        times = open("runtimes_5g.dat","a")
        times.write("*****TIME FOR CYCLE*****")
        times.write(str(tcycle))
        times.close()
    k=k+1
    position = result[0]
    f = open("status_ball_5g.txt", "w")
    f.write("****Iteration*****")
    f.write(str(k))
    f.write("****Reduced Chi2*****")
    f.write(str(result[1] * -2.0/(obspec.shape[1] / 3.0)))
    f.write("****Accept Fraction*****")
    f.write(str(sampler.acceptance_fraction))
    f.write("*****Values****")
    f.write(str(result[0]))
    f.close()
    if (k==10 or k==1000 or k==1500 or k==2000 or k==2500 or k==3000 or k==3500 or k==4000 or k==4500 or k==5000 or k==6000 or k==7000 or k==8000 or k==9000 or k==10000 or k==11000 or k==12000 or k==15000 or k==18000 or k==21000 or k==25000) or k == 30000 or k == 35000 or k == 40000 or k == 45000 or k == 50000:
        chain=sampler.chain
	lnprob=sampler.lnprobability
	output=[chain,lnprob]
	pickle.dump(output,open("/nobackup/bburning/MCMC_pow_thick5g.pic","wb"))
	pickle.dump(chain[:,k-1,:], open('MCMC_last_powthick5g.pic','wb'))


chain=sampler.chain
lnprob=sampler.lnprobability
output=[chain,lnprob]
pickle.dump(output,open("/nobackup/bburning/2m2224_5gas.pic","wb"))
pickle.dump(chain[:,-1,:], open('MCMC_last_powthick5g.pic','wb'))



# get rid of problematic bit of sampler object
del sampler.__dict__['pool']

def save_object(obj, filename):
    with open(filename, "wb") as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

pool.close()

save_object(sampler,"/nobackup/bburning/2M2224_5gas.pk1")




