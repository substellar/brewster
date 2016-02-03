#!/usr/bin/env python

"""This is Brewster: the golden retriever of smelly atmospheres"""

import multiprocessing
import time
import numpy as np
import scipy as sp
import emcee
import testkit
import ciamod
import os
import gc
import sys
import pickle
from scipy.io.idl import readsav
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
coarsePress = 1000.* pow(10,logcoarsePress)
logfinePress = np.arange(-4.0, 2.4, 0.1)
finePress = 1000.* pow(10,logfinePress)
# forward model wants pressure in mbar
press = finePress
nprof = coarsePress.size


#r2d2 = 1.
#logg = 5.0
#dlam = 0.
w1 = 1.0
w2 = 10.0
npatches = 1
pcover = np.array([npatches],dtype='f']
pcover[:] = 1.0
do_clouds = np.array([npatches],dtype='i']
do_clouds[:] = 1
cloudnum = np.array([npatches],dtype='i')
cloudnum[:] = 3
cloudtype = 2
use_disort = 0 

# hardwired FWHM of data in microns
fwhm = 0.005

# now the linelist
# Set up number of gases, and point at the lists. see gaslist.dat
ngas = 4
gasnum = np.asfortranarray(np.array([1,2,6,3],dtype='i'))
lists = ["/nobackup/bburning/Linelists/H2O_xsecs.pic","/nobackup/bburning/Linelists/ch4_xsecs.pic","/nobackup/bburning/Linelists/nh3_xsecs.pic","/nobackup/bburning/Linelists/h2s_xsecs.pic"]
# get the basic framework from water list
rawwavenum, inpress, inlinetemps, inlinelist = pickle.load( open('/nobackup/bburning/Linelists/H2O_xsecs.pic', "rb" ) )
inpress=1000.*inpress
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

# cloudparams is structured array with 5 entries
# each one has a patch*cloud entries
# ^^^^  this won't work with EMCEE
#cloudparams = np.ones(5)
# 5 entries in cloudparams for simple slab model are:
# 0) log10(number density)
# 1) top layer id (or pressure)
# 2) base ID (these are both in 61 layers)
# 3) rg
# 4) rsig
#cloudparams[0] = -20.
#cloudparams[1] = 10
#cloudparams[2] = 12
#cloudparams[3] = 1e-4
#cloudparams[4] = 1e-5
# hardwired gas and cloud IDs


# Get the cia bits
tmpcia, ciatemps = ciamod.read_cia("CIA_DS_aug_2015.dat",inwavenum)
cia = np.asfortranarray(np.empty((4,ciatemps.size,nwave)),dtype='float32')
cia[:,:,:] = tmpcia[:,:,:nwave] 
ciatemps = np.asfortranarray(ciatemps, dtype='float32')



# get the observed spectrum
obspec = np.asfortranarray(np.loadtxt("G570D_2MHcalib.dat",dtype='d',unpack='true'))


# NEXT LINE ADDS SYTEMATIC TO SPECTRUM FOR Log F retrieval
#obspec[1,:] = obspec[1,:] + 0.1*(min(obspec[2,10::3]**2))


runargs = pcover, cloudtype,do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec

# now set up the EMCEE stuff
ndim, nwalkers = (nprof + (ngas-1) + 10), 540

# If we want fresh guess set to 0, total inherit the previous set 1, inherit plus randomise the VMRs. 2.
fresh = 0
p0 = np.empty([nwalkers,ndim])
if (fresh == 0):
    p0[:,0] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 3.5 # H2O
    p0[:,1] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 3.5 # CH4
    p0[:,2] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 3.6 # NH3
    p0[:,3] = (0.5*np.random.randn(nwalkers).reshape(nwalkers)) - 5.0  #H2S
    p0[:,4] = (3*np.random.randn(nwalkers).reshape(nwalkers))-19. #ref density
    p0[:,5] = 100. * random.rand(nwalkers).reshape(nwalkers) #cloud top
    p0[:,6] = 10. * random.rand(nwalkers).reshape(nwalkers) #scale height
    #  radii in um
    p0[:,7] = 10. * random.rand(nwalkers).reshape(nwalkers) # rg
    p0[:,8] = 1.0 * random.rand(nwalkers).reshape(nwalkers) # rsig
    p0[:,9] = np.random.rand(nwalkers).reshape(nwalkers) + 4.2
    p0[:,10] = 0.95 + (np.random.randn(nwalkers).reshape(nwalkers) * 0.1)
    p0[:,11] = np.random.randn(nwalkers).reshape(nwalkers) * 0.001
    p0[:,12] =  50. + (np.random.randn(nwalkers).reshape(nwalkers))
    p0[:,13] = np.log10((np.random.rand(nwalkers).reshape(nwalkers) * (max(obspec[2,:]**2)*(10. - 0.01))) + (0.01*min(obspec[2,10::3]**2)))
    BTprof = np.loadtxt("BTtemp800_45_13.dat")
    for i in range(0,nwalkers):
        p0[i,14:] = (BTprof - 300.) + (100. * np.random.rand())

    #for i in range(0,13):
    #    p0[:,i+14] = (BTprof[i] - 200.) + (150. * np.random.randn(nwalkers).reshape(nwalkers))
    #toffset = 200.* np.random.randn(nwalkers).reshape(nwalkers)
    #for i in range (11,10+nprof):
    #    p0[:,i] = 750.+toffset + ((coarsePress[i-11]/100.)**1.1)
    #p0[:,23] = 3000. + 200 * np.random.rand(nwalkers).reshape(nwalkers)

if (fresh != 0):
    fname='MCMC_last.pic'
    pic=pickle.load(open(fname,'rb'))
    p0=pic
    if (fresh == 2):
        for i in range(0,7):
            p0[:,i] = (np.random.rand(nwalkers).reshape(nwalkers)*0.5) + p0[:,i]


    
# Now we set up the MPI bits
pool=MPIPool(loadbalance=True)
if not pool.is_master():
	pool.wait()
	sys.exit(0)

sampler = emcee.EnsembleSampler(nwalkers, ndim, testkit.lnprob, args=(runargs),pool=pool)
#'''
# run the sampler
print "running the sampler"
#sampler.run_mcmc(p0, 100)
clock = np.empty(30000)
k=0
times = open("runtimes.dat","w")
times.close()
for result in sampler.sample(p0, iterations=12000):
    clock[k] = time.clock()
    if (k > 1):
        tcycle = clock[k] - clock[k-1]
        times = open("runtimes.dat","a")
        times.write("*****TIME FOR CYCLE*****")
        times.write(str(tcycle))
        times.close()
    k=k+1
    position = result[0]
    f = open("status_ball.txt", "w")
    f.write("****Iteration*****")
    f.write(str(k))
    f.write("****Reduced Chi2*****")
    f.write(str(result[1] * -2.0/(obspec.shape[1] / 3.0)))
    f.write("****Accept Fraction*****")
    f.write(str(sampler.acceptance_fraction))
    f.write("*****Values****")
    f.write(str(result[0]))
    f.close()
    if (k==10 or k==1000 or k==1500 or k==2000 or k==2500 or k==3000 or k==3500 or k==4000 or k==4500 or k==5000 or k==6000 or k==7000 or k==8000 or k==9000 or k==10000 or k==11000 or k==12000 or k==15000 or k==18000 or k==21000 or k==25000):
        chain=sampler.chain
	lnprob=sampler.lnprobability
	output=[chain,lnprob]
	pickle.dump(output,open("/nobackup/bburning/MCMC.pic","wb"))
	pickle.dump(chain[:,k-1,:], open('MCMC_last.pic','wb'))


chain=sampler.chain
lnprob=sampler.lnprobability
output=[chain,lnprob]
pickle.dump(output,open("/nobackup/bburning/MCMC.pic","wb"))
pickle.dump(chain[:,-1,:], open('MCMC_last.pic','wb'))



# get rid of problematic bit of sampler object
del sampler.__dict__['pool']

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

pool.close()

save_object(sampler,'/nobackup/bburning/570D_pickle_retrieval_result.pk1')
#save_object(sampler,'570D_BTretrieval_result.pk1')


