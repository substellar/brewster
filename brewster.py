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
logcoarsePress = np.arange(-4.0, 2.5, 0.5)
coarsePress = 1000.* pow(10,logcoarsePress)
logfinePress = np.arange(-4.0, 2.5, 0.1)
finePress = 1000.* pow(10,logfinePress)
# forward model wants pressure in mbar
press = finePress
nprof = coarsePress.size
# This bit is all to cover reading in Mike's profile - we don't need this.
#array = pickle.load(open("test_H2H2_H2He_CIA_H2O.pic", "rb")) 
#leveltemp = array[0]
#levelpress = array[1]
#mikespec = np.array([array[2],array[3]],dtype='f')
#mikespec[0] = 10000.0 / mikespec[0]
#mikepress = np.empty(levelpress.size - 1,dtype='float64')
#miketemp = np.empty(leveltemp.size -1, dtype='float64')
#for i in range(0,mikepress.size):
#    mikepress[i] = np.sqrt(levelpress[i] * levelpress[i+1])
#mtfit = interp1d(np.log10(levelpress),leveltemp)
#miketemp = mtfit(np.log10(mikepress))
#tfit = interp1d(np.log10(mikepress),miketemp,bounds_error=False,fill_value=miketemp[miketemp.size-1])
#temp = tfit(np.log10(finePress))
#intemp = temp

# now the linelist
# Set up number of gases, and point at the lists. see gaslist.dat
ngas = 5
gasnum = np.asfortranarray(np.array([1,2,20,4,5],dtype='i'))
lists = ["/nobackup/bburning/Linelists/xsecarrH2O_1wno_500_10000.save","/nobackup/bburning/Linelists/xsecarrCH4_1wno_500_10000.save","/nobackup/bburning/Linelists/xsecarrK_new_1wno_500_10000_02.save","/nobackup/bburning/Linelists/xsecarrCO_1wno_500_10000_02.save","/nobackup/bburning/Linelists/xsecarrCO2_1wno_500_10000_02.save" ]
# get the basic framework from water list
x=readsav('/nobackup/bburning/Linelists/xsecarrH2O_1wno_500_10000.save')
inlinelist=x.xsecarr  #3D array with Nwavenubmers x Ntemps x Npressure
inlinetemps=np.asfortranarray(x.t,dtype='float64')
inpress=1000.*x.p
inwavenum=x.wno
ntemps = inlinetemps.size
npress= finePress.size
nwave = inwavenum.size
# Here we are interpolating the linelist onto our fine pressure scale. 
linelist = (np.ones([ngas,npress,ntemps,nwave],order='F')).astype('float64', order='F')
for gas in range (0,ngas):
    inlinelist=readsav(lists[gas]).xsecarr
    for i in range (0,ntemps):
        for j in range (0,nwave):
            pfit = interp1d(np.log10(inpress),np.log10(inlinelist[:,i,j]))
            linelist[gas,:,i,j] = np.asfortranarray(pfit(np.log10(finePress)))

r2d2 = 1.
logg = 4.5
dlam = 0.
w1 = 1.05
w2 = 5.0
pcover = 1.0
do_clouds = 0
use_disort = 0 
# cloudparams is structured array with 5 entries
# each one has a patch*cloud entries
cloudparams = np.ones(5)
# 5 entries in cloudparams for simple slab model are:
# 0) log10(number density)
# 1) top layer id (or pressure)
# 2) base ID (these are both in 61 layers)
# 3) rg
# 4) rsig
cloudparams[0] = -20.
cloudparams[1] = 10
cloudparams[2] = 12
cloudparams[3] = 1e-4
cloudparams[4] = 1e-5
# hardwired gas and cloud IDs
cloudnum = np.array([1],dtype='i')

# Get the cia bits
cia, ciatemps = ciamod.read_cia("CIA_DS_aug_2015.dat",inwavenum)
cia = np.asfortranarray(cia, dtype='float32')
ciatemps = np.asfortranarray(ciatemps, dtype='float32')

# hardwired FWHM of data
fwhm = 0.005
#fixvmrs = -8.0

# get the observed spectrum
obspec = np.asfortranarray(np.loadtxt("5gas_spectrum.dat",dtype='d',unpack='true'))

# NEXT LINE ADDS SYTEMATIC TO SPECTRUM FOR Log F retrieval
#obspec[1,:] = obspec[1,:] + 0.1*(min(obspec[2,10::3]**2))


runargs = w1,w2, pcover, cloudparams,r2d2,logg, dlam, do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec

# now set up the EMCEE stuff
ndim, nwalkers = 20, 96
p0 = np.empty([nwalkers,ndim])
p0[:,0] = -1.* np.random.rand(nwalkers).reshape(nwalkers) - 3.0
p0[:,1] = -1.* np.random.rand(nwalkers).reshape(nwalkers) - 3.0
p0[:,2] = -1.* np.random.rand(nwalkers).reshape(nwalkers) - 7.5
p0[:,3] =  -1.* np.random.rand(nwalkers).reshape(nwalkers) - 7.0
p0[:,4] =  -1.* np.random.rand(nwalkers).reshape(nwalkers) - 7.2
#p0[:,5] = np.log10((np.random.rand(nwalkers).reshape(nwalkers) * (min(obspec[2,10::3]**2)*(0.1 - 0.001))) + (0.001*min(obspec[2,10::3]**2)))
p0[:,5] =  50. + (np.random.randn(nwalkers).reshape(nwalkers))
p0[:,6] = (-2. *  np.random.rand(nwalkers).reshape(nwalkers)) - 2.
p0[:,7] = 200. + ( np.random.rand(nwalkers).reshape(nwalkers) * 100. )
for i in range (8,7+nprof):
    p0[:,i] = p0[:,7] + (150.*(i-7))

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
clock = np.empty(20000)
k=0
times = open("runtimes.dat","w")
times.close()
for result in sampler.sample(p0, iterations=10):
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
    f.write(str(result[1]/((obspec.shape[1] - 10.)/ 3.0) *(-2)))
    f.write("****Accept Fraction*****")
    f.write(str(sampler.acceptance_fraction))
    f.write("*****Values****")
    f.write(str(result[0]))
    f.close()
    if (k==100 or k==1000 or k==1500 or k==2000 or k==2500 or k==3000 or k==3500 or k==4000 or k==4500 or k==5000 or k==6000 or k==7000 or k==8000 or k==9000 or k==10000 or k==11000 or k==12000 or k==15000 or k==18000 or k==21000 or k==25000):
        chain=sampler.chain
	lnprob=sampler.lnprobability
	output=[chain,lnprob]
	pickle.dump(output,open("/nobackup/bburning/MCMC.pic","wb"))
	pickle.dump(chain[:,k-1,:], open('/nobackup/bburning/MCMC_last.pic','wb'))


chain=sampler.chain
lnprob=sampler.lnprobability
output=[chain,lnprob]
pickle.dump(output,open("/nobackup/bburning/MCMC.pic","wb"))
pickle.dump(chain[:,-1,:], open('/nobackup/bburning/MCMC_last.pic','wb'))



# get rid of problematic bit of sampler object
del sampler.__dict__['pool']

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

pool.close()

save_object(sampler,'/nobackup/bburning/temp_retrieval_result.pk1')
save_object(sampler,'temp_retrieval_result.pk1')


