#!/usr/bin/env python


""" McNuggets: the post-processing tool for brewster"""

import numpy as np
import scipy as sp
import forwardmodel
import ciamod
import TPmod
import cloud
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

    pcover, cloudtype,cloudparams,do_clouds,gasnum,cloudnum,\
        inlinetemps,coarsePress,press,inwavenum,linelist,cia,\
        ciatemps,use_disort,proftype, do_fudge = runargs

    
    if (gasnum[gasnum.size-1] == 21):
        ng = gasnum.size - 1
    elif (gasnum[gasnum.size-1] == 23):
        ng = gasnum.size - 2
    else:
        ng = gasnum.size
        
    invmr = theta[0:ng]
    logg = theta[ng]
    r2d2 = theta[ng+1]
    dlam = theta[ng+2]
    if (do_fudge == 1):
        logf = theta[ng+3]
        pc = ng + 4
    else:
        pc = ng + 3
    nc = 0
    if (do_clouds == 1):
        if ((cloudtype == 2) and (cloudnum == 99)):
            nc = 3
            cloudparams[1:4] = theta[pc:pc+nc]
            cloudparams[4] = 0.0
        elif ((cloudtype == 1) and (cloudnum == 99)):
            nc = 4
            cloudparams[0:4] = theta[pc:pc+nc]
            cloudparams[4] = 0.0
        elif ((cloudtype == 2) and (cloudnum == 89)):
            nc  = 4
            cloudparams[1:5] = theta[pc:pc+nc]
        elif ((cloudtype == 1) and (cloudnum == 89)):
            nc = 5
            cloudparams = theta[pc:pc+nc]
        elif(cloudnum < 50):
            nc = 5
            cloudparams = theta[pc:pc+nc]
        
    if (proftype == 1):
        gam = theta[pc+nc]
        intemp = theta[pc+1+nc:]
    elif (proftype == 2 or proftype ==3):
        intemp = theta[pc+nc:]
    else:
        raise ValueError("not valid profile type %proftype" % (char, string))

    # Hard code nlayers
    nlayers = press.shape[0]
    npatch = cloudparams.shape[0]
    # set the profile
    temp = TPmod.set_prof(proftype,coarsePress,press,intemp)
 
    # get the ngas for forward model (ngas, not ng
    if (gasnum[gasnum.size-1] == 21):
        ngas = invmr.shape[0] + 1
    elif (gasnum[gasnum.size-1] == 23):
        ngas = invmr.shape[0] + 2
    else:
        ngas = invmr.shape[0]
    # now loop through gases and get VMR for model
    # check if its a fixed VMR or a profile
    # VMR is log10(VMR) !!!
    logVMR = np.empty((ngas,nlayers),dtype='d')
    alkratio = 16.2 #  from Asplund et al (2009)
    if invmr.size > invmr.shape[0]:
        # this case is a profile
        # now sort Na and K
        tmpvmr = np.empty((ngas,nlayers),dtype='d')
        if (gasnum[gasnum.size-1] == 21):
            tmpvmr[0:(ngas-2),:] = invmr[0:(ngas-2),:]
            tmpvmr[ngas-2,:] = np.log10(10.**invmr[ngas-2,:] / (alkratio+1.)) # K
            tmpvmr[ngas-1,:] = np.log10(10.**invmr[ngas-2,:] * (alkratio / (alkratio+1.))) # Na                                
        elif (gasnum[gasnum.size-1] == 23):
            #f values are ratios between Na and (K+Cs) and K and Cs respectively
            f1 = 1.348
            f2 = 8912.5
            tmpvmr[0:(ngas-3),:] = invmr[0:(ngas-3),:]
            tmpvmr[ngas-1,:] = np.log10(10.**invmr[ngas-3,:] / ((f1+1)*(f2+1))) # Cs
            tmpvmr[ngas-2,:] = np.log10(10.**invmr[ngas-3,:] * (f1 /(f1+1)) ) # Na
            tmpvmr[ngas-3,:] = np.log10(10.**invmr[ngas-3,:] - 10.**tmpvmr[ngas-2,:] - 10.**tmpvmr[ngas-1,:]) #K 
        else:
            tmpvmr[0:ngas,:] = invmr[0:ngas,:]
            
        for i in range(0,ngas):
            vfit = sp.interpolate.splrep(np.log10(coarsepress),tmpvmr[i,:],s=0)
            logVMR[i,:] = sp.interpolate.splev(np.log10(press),vfit,der=0)
    else:
        # This caseis fixed VMR
        # now sort Na and K
        tmpvmr = np.empty(ngas,dtype='d')
        if (gasnum[gasnum.size-1] == 21):
            tmpvmr[0:(ngas-2)] = invmr[0:(ngas-2)]
            tmpvmr[ngas-2] = np.log10(10.**invmr[ngas-2] / (alkratio+1.)) # K
            tmpvmr[ngas-1] = np.log10(10.**invmr[ngas-2] * (alkratio / (alkratio+1.))) # Na
        elif (gasnum[gasnum.size-1] == 23):
            #f values are ratios between Na and (K+Cs) and K and Cs respectively
            f1 = 1.348
            f2 = 8912.5
            tmpvmr[0:(ngas-3)] = invmr[0:(ngas-3)]
            tmpvmr[ngas-1] = np.log10(10.**invmr[ngas-3] / ((f1+1)*(f2+1))) # Cs
            tmpvmr[ngas-2] = np.log10(10.**invmr[ngas-3] * (f1 /(f1+1)) ) # Na
            tmpvmr[ngas-3] = np.log10(10.**invmr[ngas-3] - 10.**tmpvmr[ngas-2] - 10.**tmpvmr[ngas-1]) #K   
        else:
            tmpvmr[0:ngas] = invmr[0:ngas]
            
        for i in range(0,ngas):                              
            logVMR[i,:] = tmpvmr[i]

    # now need to translate cloudparams in to cloud profile even
    # if do_clouds is zero..
    # 5 entries for cloudparams for simple slab model are:
    # 0) log10(number density / gas number density)
    # 1) top layer id (or pressure)
    # 2) base ID (these are both in 64 layers)
    # 3) rg
    # 4) rsig
    # in the case of a simple mixto cloud (i.e. cloudnum = 99) we have:
    # 0) ndens = dtau
    # 1) top layer ID
    # 2) bottom later ID
    # 3) rg = albedo
    # 4) rsig = asymmetry
    if (npatch > 1 or do_clouds == 1):
        cloudprof,cloudrad,cloudsig = cloud.atlas(do_clouds,cloudnum,cloudtype,cloudparams,press)
        npatch = cloudprof.shape[0]
        ncloud = cloudprof.shape[1]
    else:
        npatch = 1
        ncloud = 1
        cloudrad = np.ones((npatch,nlayers,ncloud),dtype='d')
        cloudsig = np.ones_like(cloudrad)
        cloudprof = np.ones_like(cloudrad)

    cloudprof = np.asfortranarray(cloudprof,dtype = 'float64')
    cloudrad = np.asfortranarray(cloudrad,dtype = 'float64')
    cloudsig = np.asfortranarray(cloudsig,dtype = 'float64')
    pcover = np.asfortranarray(pcover,dtype = 'float32')
    cloudnum = np.asfortranarray(cloudnum,dtype='i')
    do_clouds = np.asfortranarray(do_clouds,dtype = 'i')

        
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
    parallax = 11.35
    sigpi = 0.14
    sigphot = 0.03
    
    sigR2D2 = sigphot * r2d2 * (-1./2.5)* np.log(10.)

    sigD = sigpi * 3.086e16
    D = parallax * 3.086e16

    R = np.sqrt(((np.random.randn() * sigR2D2)+ r2d2)) \
        * ((np.random.randn()* sigD) + D)

    g = (10.**logg)/100.

    # and mass

    M = (R**2 * g/(6.67E-11))/1.898E27
    R = R / 71492e3

    result = np.concatenate((theta,np.array([t_ff, R, M])),axis=0)
    
    return result




# Start by getting the chain, then set up the run arguments, then loop for Teff etc

# how many samples are we using?

with open('/nobackup/bburning/2M2224_powcloud_RFalks_Cs.pk1', 'rb') as input:
    sampler = pickle.load(input) 

nwalkers = sampler.chain.shape[0]
niter = sampler.chain.shape[1]
ndim = sampler.chain.shape[2]

samples = sampler.chain[:,niter-3000:,:].reshape((-1, ndim))

#samples = samples[1500:2500,:]
slen = samples.shape[0]
#samplus = np.zeros([slen,ndim+3])

#samplus[:,0:ndim] = samples


# set up run arguments


# set up pressure grids
logcoarsePress = np.arange(-4.0, 2.5, 0.53)
#logcoarsePress = np.arange(-4.0, 3.0, 0.5)
coarsePress = pow(10,logcoarsePress)
logfinePress = np.arange(-4.0, 2.4, 0.1)
finePress = pow(10,logfinePress)
press = finePress

w1 = 0.7
w2 = 20.0


npatches = 1
nclouds = 1
pcover = np.ones([npatches],dtype='f')
pcover[:] = 1.0
do_clouds = np.zeros([npatches],dtype='i')
do_clouds[:] = 1
cloudnum = np.zeros([npatches,nclouds],dtype='i')
cloudnum[:,:] = 89
cloudtype = np.array([npatches],dtype='i')
cloudtype[:] = 1

use_disort = 0 

# use the fudge factor?
do_fudge = 1

# Set the profile type
proftype = 2

# now the linelist
# Set up number of gases, and point at the lists. see gaslist.dat
ngas = 12
gasnum = np.asfortranarray(np.array([1,4,5,7,8,9,10,11,12,20,21,23],dtype='i'))
lists = ["/nobackup/bburning/Linelists/H2O_xsecs.pic","/nobackup/bburning/Linelists/co_xsecs.pic","/nobackup/bburning/Linelists/co2_xsecs.pic","/nobackup/bburning/Linelists/tio_xsecs.pic","/nobackup/bburning/Linelists/vo_xsecs.pic","/nobackup/bburning/Linelists/cah_xsecs.pic","/nobackup/bburning/Linelists/crh_xsecs.pic" ,"/nobackup/bburning/Linelists/feh_xsecs.pic","/nobackup/bburning/Linelists/mgh_xsecs.pic","/nobackup/bburning/Linelists/K_xsecs.pic","/nobackup/bburning/Linelists/Na_xsecs.pic","/nobackup/bburning/Linelists/Cs_xsecs.pic"]
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


# place holder values for cloudparams
cloudparams = np.ones([5],dtype='d')
cloudparams[0] = -20.
cloudparams[1] = 10
cloudparams[2] = 12
cloudparams[3] = 1e-4
cloudparams[4] = 1e-5


runargs =  pcover, cloudtype,cloudparams,do_clouds,gasnum,cloudnum,\
        inlinetemps,coarsePress,press,inwavenum,linelist,cia,\
        ciatemps,use_disort,proftype,do_fudge

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
    results.append(teffRM(samples[job,0:ndim],runargs))

# Gather results on rank 0.
results = MPI.COMM_WORLD.gather(results, root=0)

if COMM.rank == 0:
    # Flatten list of lists.
    results = [_i for tmp in results for _i in tmp]

samplus = np.array(results)




def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

save_object(samplus,'/nobackup/bburning/2M2224_CsPowCloud_postprod.pk1')

