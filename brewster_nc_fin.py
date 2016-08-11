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


do_clouds = np.array([0],dtype='i')


# CURRENTLY ONLY COPE WITH ONE CLOUDY PATCH.
#SO MAKE ALL CLOUD PARAMETERS THE SAME FOR EASE OF PROCESSING 

cloudnum = np.zeros([npatches,nclouds],dtype='i')
cloudnum[:,:] = 99
cloudtype = np.asfortranarray(np.array([1]),dtype='i')

use_disort = 0 

# use the fudge factor?
do_fudge = 1

# Set the profile type
proftype = 2

# now the linelist
# Set up number of gases, and point at the lists. see gaslist.dat
ngas = 11
gasnum = np.asfortranarray(np.array([1,4,5,7,8,9,10,11,12,20,21],dtype='i'))
lists = ["/nobackup/bburning/Linelists/H2O_xsecs.pic","/nobackup/bburning/Linelists/co_xsecs.pic","/nobackup/bburning/Linelists/co2_xsecs.pic","/nobackup/bburning/Linelists/tio_xsecs.pic","/nobackup/bburning/Linelists/vo_xsecs.pic","/nobackup/bburning/Linelists/cah_xsecs.pic","/nobackup/bburning/Linelists/crh_xsecs.pic" ,"/nobackup/bburning/Linelists/feh_xsecs.pic","/nobackup/bburning/Linelists/mgh_xsecs.pic","/nobackup/bburning/Linelists/K_Mike_xsecs.pic","/nobackup/bburning/Linelists/Na_Mike_xsecs.pic"]
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



# get the observed spectrum
obspec = np.asfortranarray(np.loadtxt("2M2224_mkoJcalib_trim.dat",dtype='d',unpack='true'))



# place holder values for cloudparams
cloudparams = np.ones([5],dtype='d')
cloudparams[0] = 1.0
cloudparams[1] = 1.0
cloudparams[2] = 0.5
cloudparams[3] = 0.5
cloudparams[4] = 0.0


runargs = dist, cloudtype,cloudparams,do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge

# now set up the EMCEE stuff

ndim  = 19 #((ngas-1) + 9 + 5)
nwalkers = ndim * 16
#int(((ndim * ndim) // 2) * 2)


# If we want fresh guess set to 0, total inherit the previous set 1, inherit plus randomise the VMRs. 2.
fresh = 0
p0 = np.empty([nwalkers,ndim])
if (fresh == 0):
    p0[:,0] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 3.5 # H2O
    p0[:,1] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 4.0 # CO
    p0[:,2] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # CO2
    p0[:,3] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # TiO
    p0[:,4] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # VO
    p0[:,5] = (1.0*np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # CaH
    p0[:,6] = (1.0*np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # CrH
    p0[:,7] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # FeH
    p0[:,8] = (1.0*np.random.randn(nwalkers).reshape(nwalkers)) - 8.0 # MgH
    p0[:,9] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 5.5 # Na+K
    p0[:,10] = np.random.rand(nwalkers).reshape(nwalkers) + 4.2
    p0[:,11] =  1.0e-20 + np.random.rand(nwalkers).reshape(nwalkers) * 5.e-20
    p0[:,12] = np.random.randn(nwalkers).reshape(nwalkers) * 0.001
    p0[:,13] = np.log10((np.random.rand(nwalkers).reshape(nwalkers) * (max(obspec[2,:]**2)*(0.1 - 0.01))) + (0.01*min(obspec[2,10::3]**2)))
    # some cloud bits now. We're just doing grey cloud, tau so need pressure of top where plus cloud height (in dex), SSA, don't need GG
    # And now the T-P params
    p0[:,14] = 0.6+ 0.1*np.random.rand(nwalkers).reshape(nwalkers)
    p0[:,15] = 0.1 +0.05*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,16] = 0.2+ 0.05*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,17] = 2.+ 0.2*np.random.randn(nwalkers).reshape(nwalkers)
    p0[:,18] = 4000. + (1000.*  np.random.rand(nwalkers).reshape(nwalkers))

    for i in range (0,nwalkers):
        while True:
            Tcheck = TPmod.set_prof(proftype,coarsePress,press,p0[i,14:])
            if (min(Tcheck) > 1.0):
                break
            else:
                p0[i,14] = 0.6+ 0.1*np.random.rand()
                p0[i,15] = 0.1+ 0.05*np.random.randn()
                p0[i,16] = 0.2+ 0.05*np.random.randn()
                p0[i,17] = 2. + 0.2*np.random.randn()
                p0[i,18] = 4000. + (1000.*  np.random.rand())

    
if (fresh != 0):
    fname='MCMC_last.pic'
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
times = open("runtimes_nc.dat","w")
times.close()
pos,prob,state = sampler.run_mcmc(p0,10000)
sampler.reset()
for result in sampler.sample(pos, iterations=30000):
    clock[k] = time.clock()
    if (k > 1):
        tcycle = clock[k] - clock[k-1]
        times = open("runtimes_nc.dat","a")
        times.write("*****TIME FOR CYCLE*****")
        times.write(str(tcycle))
        times.close()
    k=k+1
    position = result[0]
    f = open("status_ball_nc.txt", "w")
    f.write("****Iteration*****")
    f.write(str(k))
    f.write("****Reduced Chi2*****")
    f.write(str(result[1] * -2.0/(obspec.shape[1] / 3.0)))
    f.write("****Accept Fraction*****")
    f.write(str(sampler.acceptance_fraction))
    f.write("*****Values****")
    f.write(str(result[0]))
    f.close()
    if (k==10 or k==1000 or k==1500 or k==2000 or k==2500 or k==3000 or k==3500 or k==4000 or k==4500 or k==5000 or k==6000 or k==7000 or k==8000 or k==9000 or k==10000 or k==11000 or k==12000 or k==15000 or k==18000 or k==21000 or k==25000) or k == 30000 or k == 35000 or k == 40000:
        chain=sampler.chain
	lnprob=sampler.lnprobability
	output=[chain,lnprob]
	pickle.dump(output,open("/nobackup/bburning/2M2224_MCMC_nc_fin.pic","wb"))
	pickle.dump(chain[:,k-1,:], open('MCMC_last.pic','wb'))


chain=sampler.chain
lnprob=sampler.lnprobability
output=[chain,lnprob]
pickle.dump(output,open("/nobackup/bburning/2M2224_MCMC_nc_fin.pic","wb"))
pickle.dump(chain[:,-1,:], open('MCMC_last.pic','wb'))



# get rid of problematic bit of sampler object
del sampler.__dict__['pool']

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

pool.close()

save_object(sampler,'/nobackup/bburning/2M2224_mikeConv_nc_fudge.pk1')
#save_object(sampler,'570D_BTretrieval_result.pk1')


