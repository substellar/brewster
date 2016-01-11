#!/usr/bin/env python


""" McNuggets: the post-processing tool for brewster"""

import numpy as np
import scipy as sp
import forwardmodel
import ciamod
import pickle
from scipy.io.idl import readsav
from scipy import interpolate
from scipy.interpolate import interp1d
from mpi4py import MPI

__author__ = "Ben Burningham"
__copyright__ = "Copyright 2015 - Ben Burningham"
__credits__ = ["Ben Burningham","The EMCEE DOCS"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ben Burningham"
__email__ = "burninghamster@gmail.com"
__status__ = "Development"



def teffRM(theta,runargs):
    
    invmr = theta[0:7]
    logf = 0.0 #theta[5]
    #logbeta = theta[6]
    logg = theta[7]
    r2d2 = theta[8]
    dlam = theta[9]
    gam = theta[10]
    intemp = theta[11:]
    pcover, cloudparams,do_clouds,gasnum,cloudnum,inlinetemps,\
        coarsePress,press,inwavenum,linelist,cia,ciatemps,\
        use_disort = runargs


    # interp temp onto finer grid coarsePress => press
    # spline fit with max smoothing
    tfit = sp.interpolate.splrep(np.log10(coarsePress),np.log10(intemp),s=0)
    temp = 10.**(np.asfortranarray(sp.interpolate.splev(np.log10(press),tfit,der=0),dtype='d'))

    # get the ngas
    ngas = invmr.shape[0] + 1
    # Hard code nlayers
    nlayers = press.shape[0]
    # interp temp onto finer grid coarsePress => press
    # spline fit with max smoothing
    tfit = sp.interpolate.splrep(np.log10(coarsePress),intemp,s=0)
    temp = np.asfortranarray(sp.interpolate.splev(np.log10(press),tfit,der=0),dtype='d')
    # now loop through gases and get VMR for model
    # check if its a fixed VMR or a profile
    # VMR is log10(VMR) !!!
    logVMR = np.empty((ngas,nlayers),dtype='d')
    alkratio = 16.2 #  from Asplund et al (2009)
    if invmr.size > invmr.shape[0]:
        # now sort Na and K
        tmpvmr = np.empty((ngas,nlayers),dtype='d')
        tmpvmr[0:(ngas-2),:] = invmr[0:(ngas-2),:]
        tmpvmr[ngas-2,:] = np.log10(10.**invmr[ngas-2,:] / (alkratio+1.))
        tmpvmr[ngas-1,:] = np.log10(10.**invmr[ngas-2,:] * (alkratio / (alkratio+1.)))                                
        for i in range(0,ngas):
            vfit = sp.interpolate.splrep(np.log10(coarsepress),tmpvmr[i,:],s=0)
            logVMR[i,:] = sp.interpolate.splev(np.log10(press),vfit,der=0)
    else:
        # now sort Na and K
        tmpvmr = np.empty(ngas,dtype='d')
        tmpvmr[0:(ngas-2)] = invmr[0:(ngas-2)]
        tmpvmr[ngas-2] = np.log10(10.**invmr[ngas-2] / (alkratio+1.))
        tmpvmr[ngas-1] = np.log10(10.**invmr[ngas-2] * (alkratio / (alkratio+1.)))
        for i in range(0,ngas):                              
            logVMR[i,:] = tmpvmr[i]

    # now need to translate cloudparams in to cloud profile even
    # if do_clouds is zero..
    # 5 entries for cloudparams for simple slab model are:
    # 0) log10(number density)
    # 1) top layer id (or pressure)
    # 2) base ID (these are both in 61 layers)
    # 3) rg
    # 4) rsig
    if (do_clouds == 1):
        npatch = cloudparams.shape[0]
        ncloud = cloudparams.shape[1]
        cloudrad = np.empty((npatch,nlayers,ncloud),dtype='d')
        cloudsig = np.empty_like(cloudrad)
        cloudprof = np.zeros_like(cloudrad)
        ndens= np.reshape(cloudparams['f0'],(npatch,ncloud))
        c1 = np.reshape(cloudparams['f1'],(npatch,ncloud))
        c2 = np.reshape(cloudparams['f2'],(npatch,ncloud))
        rad = np.reshape(cloudparams['f3'],(npatch,ncloud))
        sig = np.reshape(cloudparams['f4'],(npatch,ncloud))
        for i in range(0, npatch):
            for j in range(0, ncloud):
                b1 = c1[i,j] - 1
                b2 = c2[i,j] -1 
                cloudprof[i,b1:b2+1,j] = ndens[i,j]
                cloudrad[i,:,j] = rad[i,j]
                cloudsig[i,:,j] = sig[i,j]        
    else:
        npatch = 1
        ncloud = 1
        cloudrad = np.ones((npatch,nlayers,ncloud),dtype='d')
        cloudsig = np.ones_like(cloudrad)
        cloudprof = np.ones_like(cloudrad)

    # now we can call the forward model
    outspec = forwardmodel.marv(temp,logg,r2d2,gasnum,logVMR,pcover,do_clouds,cloudnum,cloudrad,cloudsig,cloudprof,inlinetemps,press,inwavenum,linelist,cia,ciatemps,use_disort)

    wave = np.array(outspec[0,::-1])
    flux = np.array(outspec[1,::-1])

    # now calculate Fbol by summing the spectrum across its wave bins
    fbol = 0.0
    
    for j in range(1, (wave.size - 1)):
        sbin = ((wave[j] - wave[j-1]) + (wave[j+1] - wave[j])) / 2. 
    
        fbol = (sbin * flux[j]) + fbol

    # now get T_eff
    t_ff = ((fbol/(r2d2 * 5.670367e-8))**(1./4.))

    # and Radius
    parallax = 5.84
    sigpi = 0.03
    sigphot = 0.02
    
    sigR2d2 = sigphot * r2d2 * (-1./2.5)* log(10.)

    sigD = sigpi * 3.086e16
    D = parallax * 3.086e16

    R = np.sqrt(((np.random.randn() * sigR2D2)+ r2d2)) \
        * ((np.random.randn()* sigD) + D)
    
    # and mass

    M = (R**2 * g/(6.67E-11))/1.898E27

    
    return t_ff, R, M



# Start by getting the chain, then set up the run arguments, then loop for Teff etc

# how many samples are we using?

with open('temp_retrieval_result.pk1', 'rb') as input:
    sampler = pickle.load(input) 

ndim = sampler.chain.shape[2]
samples = sampler.chain[:,15000:,:].reshape((-1, ndim))

slen = samples.shape[0]
samplus = np.zeros([slen,ndim+3])

samplus[:,0:ndim] = samples


# set up run arguments

w1 = 1.0
w2 = 20.0
pcover = 1.0
do_clouds = 0
use_disort = 0


# set up pressure grids
logcoarsePress = np.arange(-4.0, 2.5, 0.53)
#logcoarsePress = np.arange(-4.0, 3.0, 0.5)
coarsePress = 1000.* pow(10,logcoarsePress)
logfinePress = np.arange(-4.0, 2.4, 0.1)
finePress = 1000.* pow(10,logfinePress)
press = finePress

# now the linelist
# Set up number of gases, and point at the lists. see gaslist.dat
ngas = 8
gasnum = np.asfortranarray(np.array([1,2,4,5,6,3,20,21],dtype='i'))
lists = ["/nobackup/bburning/Linelists/xsecarrH2O_1wno_500_10000.save","/nobackup/bburning/Linelists/xsecarrCH4_1wno_500_10000.save","/nobackup/bburning/Linelists/xsecarrCO_1wno_500_10000_02.save","/nobackup/bburning/Linelists/xsecarrCO2_1wno_500_10000_02.save" ,"/nobackup/bburning/Linelists/xsecarrNH3_1wno_500_10000.save","/nobackup/bburning/Linelists/xsecarrH2S_1wno_500_10000.save","/nobackup/bburning/Linelists/xsecarrK_new_1wno_500_10000_02.save","/nobackup/bburning/Linelists/xsecarrNa_new_1wno_500_10000_02.save"]
# get the basic framework from water list
x=readsav('/nobackup/bburning/Linelists/xsecarrH2O_1wno_500_10000.save')
inlinelist=x.xsecarr  #3D array with Nwavenubmers x Ntemps x Npressure
inlinetemps=np.asfortranarray(x.t,dtype='float64')
inpress=1000.*x.p
rawwavenum=x.wno
wn1 = 10000./w2
wn2 = 10000. / w1
inwavenum = np.asfortranarray(rawwavenum[np.where(np.logical_not(np.logical_or(rawwavenum[:] > wn2, rawwavenum[:] < wn1)))],dtype='float64')
ntemps = inlinetemps.size
npress= finePress.size
nwave = inwavenum.size
r1 = np.amin(np.where(np.logical_not(np.logical_or(rawwavenum[:] > wn2, rawwavenum[:] < wn1))))
r2 = np.amax(np.where(np.logical_not(np.logical_or(rawwavenum[:] > wn2, rawwavenum[:] < wn1))))

# Here we are interpolating the linelist onto our fine pressure scale. 
linelist = (np.ones([ngas,npress,ntemps,nwave],order='F')).astype('float64', order='F')
for gas in range (0,ngas):
    inlinelist=readsav(lists[gas]).xsecarr
    for i in range (0,ntemps):
        for j in range (r1,r2+1):
            pfit = interp1d(np.log10(inpress),np.log10(inlinelist[:,i,j]))
            linelist[gas,:,i,(j-r1)] = np.asfortranarray(pfit(np.log10(finePress)))

linelist[np.isnan(linelist)] = -50.0

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
tmpcia, ciatemps = ciamod.read_cia("CIA_DS_aug_2015.dat",inwavenum)
cia = np.asfortranarray(np.empty((4,ciatemps.size,nwave)),dtype='float32')
cia[:,:,:] = tmpcia[:,:,:nwave] 
ciatemps = np.asfortranarray(ciatemps, dtype='float32')

runargs = pcover, cloudparams,do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec

# set up parallel bits

COMM = MPI.COMM_WORLD

def split(container, count):
#    """
#    Simple function splitting a container into equal length chunks.
#    Order is not preserved but this is potentially an advantage depending on
#    the use case.
#    """
    return [container[_i::count] for _i in range(count)]


# Collect whatever has to be done in a list. Here we'll just collect a list of
# numbers. Only the first rank has to do this.
if COMM.rank == 0:
    jobs = list(range(slen))
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
    results.append(teffRM(samplus[job,0:ndim],runargs))

# Gather results on rank 0.
results = MPI.COMM_WORLD.gather(results, root=0)

if COMM.rank == 0:
    # Flatten list of lists.
    results = [_i for tmp in results for _i in tmp]

samplus[:,ndim:ndim+2] = results




def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

save_object(samplus,'570D_postproductchain.pk1')

