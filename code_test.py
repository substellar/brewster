# This script tests the forward model in the context of cloudless T8 dwarf
# parameters a drawn from a retrieval on G570D
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mgimg
import matplotlib.colors as colors
import scipy as sp
import numpy as np
import emcee
import testkit
import corner
import pickle as pickle
import forwardmodel
import ciamod
import TPmod
import cloud
import band
import brewtools
from astropy.convolution import convolve, convolve_fft
from astropy.convolution import Gaussian1DKernel
from scipy.io.idl import readsav
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from bensconv import spex_non_uniform
from bensconv import conv_uniform_R
from bensconv import conv_uniform_FWHM



def NoCloud_Tdwarf(xpath):
     # start with the wavelength range
     w1 = 1.0
     w2 = 2.5
     
     
     # How many patches & clouds do we want??
     # Must be at least 1 of each, but can turn off cloud below
     npatches = 1
     nclouds = 1
     
     # set up array for setting patchy cloud answers
     do_clouds = np.zeros([npatches],dtype='i')
     
     # Which patches are cloudy
     do_clouds[:] = 0
     
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
     cloudnum[:,0] = 99
     cloudtype[:,0] = 1
     
     #cloudnum[:,1] = 1
     #cloudtype[:,1] = 2
     
     #The cloud properties are set in cloudparams
     # These are not used in this notebook
     cloudparams = np.zeros([5,npatches,nclouds],dtype='d')
     
     #What's the cloud covering fraction?
     pcover = np.empty(npatches)
     pcover[:] = 1.0
     
     # Are we doing H- bound-free, free-free continuum opacities?
     # (Is the profile going above approx 2000K in the photosphere?)
     do_bff = 0
     
     # are we doing chemical equilibrium?
     chemeq = 0
     
     # We're grabbing an output profile from G570D test to check everything is working right
     # max likelihood values of profile 2 parameters, or profile read in from file
     proftype = 9
     pfile = 'G570D_model_benchmark_PROFILE.dat'
     
     # set up pressure grids in bar cos its intuitive
     logcoarsePress = np.arange(-4.0, 2.5, 0.53)
     logfinePress = np.arange(-4.0, 2.4, 0.1)
     # forward model wants pressure in bar
     #logcoarsePress = np.arange(-4.0, 3.0, 0.5)
     coarsePress = pow(10,logcoarsePress)
     press = pow(10,logfinePress)
     nlayers = press.size
     
     
     # now the cross sections
     
     # Now list the gases.
     # If Na is after K, at the end of the list, alkalis will be tied
     # together at Asplund solar ratio. See Line at al (2015)
     # Else if K is after Na, they'll be separate
     
     gaslist = ['h2o','ch4','co','co2','nh3','h2s','K','Na']
     invmr = np.array([-3.27,-3.36,-7.27,-8.28,-4.73,-8.71,-5.36])
     
     ngas = len(gaslist)
     
     # gravity??
     logg = 4.89
     
     # scale factor r2d2 from distance 1 Rj radius
     r2d2 = 1.50901046e-19 #1.46166705e-19 #(71492e3)**2. / (dist * 3.086e+16)**2.
     
     dlam = 0.00258329 #2.04395999e-03 
     
     
     #some switches for alternative cross sections
     # Use Mike's Alkalis?
     malk = 0
     
     
     # set the profile
     prof = np.full(13,1700.)
     logmodP,modT = np.loadtxt(pfile,skiprows=0,usecols=(0,1),unpack=True)
     tfit = InterpolatedUnivariateSpline(logmodP,modT,k=1)
     prof = tfit(logcoarsePress)
     temp = TPmod.set_prof(proftype,coarsePress,press,prof)
     
     
     cloudprof,cloudrad,cloudsig = cloud.atlas(do_clouds,cloudnum,cloudtype,cloudparams,press)
     
     cloudprof = np.asfortranarray(cloudprof,dtype = 'float64')
     cloudrad = np.asfortranarray(cloudrad,dtype = 'float64')
     cloudsig = np.asfortranarray(cloudsig,dtype = 'float64')
     pcover = np.asfortranarray(pcover,dtype = 'float32')
     cloudnum = np.asfortranarray(cloudnum,dtype='i')
     do_clouds = np.asfortranarray(do_clouds,dtype = 'i')


     # Now we'll get the opacity files into an array
     inlinetemps,inwavenum,linelist,gasnum,nwave = testkit.get_opacities(gaslist,w1,w2,press,xpath,malk)
 

     # Get the cia bits
     tmpcia, ciatemps = ciamod.read_cia("CIA_DS_aug_2015.dat", inwavenum)
     cia = np.asfortranarray(np.empty((4, ciatemps.size, nwave)), dtype='float32')
     cia[:, :, :] = tmpcia[:, :, :nwave]
     ciatemps = np.asfortranarray(ciatemps, dtype='float32')


     # grab BFF and Chemical grids
     bff_raw,ceTgrid,metscale,coscale,gases_myP = testkit.sort_bff_and_CE(chemeq,"chem_eq_tables_P3K.pic",press,gaslist)
     

     logVMR = np.empty((ngas,nlayers),dtype='d')
     if chemeq == 1:
          # this case is a profile
          ng = invmr.shape[2]
          ngas = ng - 3
          logVMR = np.zeros([ngas,nlayers],dtype='d')
          for p in range(0,nlayers):
               for g in range(0,ng):
                    tfit = InterpolatedUnivariateSpline(ceTgrid,invmr[:,p,g])
                    if (g < 3):
                         bff[g,p] = tfit(temp[p])
                    else:
                         logVMR[g-3,p]= tfit(temp[p])
     elif (gasnum[gasnum.size-1] == 21):
          tmpvmr = np.empty(ngas,dtype='d')
          alkratio = 16.2 # from asplund et al (2009)
          tmpvmr[0:(ngas-2)] = invmr[0:(ngas-2)]
          tmpvmr[ngas-2] = np.log10(10.**invmr[ngas-2] / (alkratio+1.)) # K
          tmpvmr[ngas-1] = np.log10(10.**invmr[ngas-2] * (alkratio / (alkratio+1.)))
        
          for i in range(0,ngas):                              
               logVMR[i,:] = tmpvmr[i]
     else:
          for i in range(0,ngas):                              
               logVMR[i,:] = invmr[i]



     bff = np.zeros([3,nlayers],dtype="float64")
     # Now get the BFF stuff sorted
     if (chemeq == 0 and do_bff == 1):
          for gas in range(0,3):
               for i in range(0,nlayers):
                    tfit = InterpolatedUnivariateSpline(ceTgrid,bff_raw[:,i,gas],k=1) 
                    bff[gas,i] = tfit(temp[i])
     bff = np.asfortranarray(bff, dtype='float64')

     press = np.asfortranarray(press,dtype='float32')
     temp = np.asfortranarray(temp,dtype='float64')
     logVMR = np.asfortranarray(logVMR,dtype='float64')


     # switches for tau spec, phot spec and contribution function
     ophot = 0
     clphot = 0
     make_cf = 0
     # other switches; ignore
     use_disort = 0
     # now we can call the forward model
     tmpoutspec,tmpclphotspec,tmpophotspec,cf = forwardmodel.marv(temp,logg,r2d2,gasnum,logVMR,pcover,do_clouds,cloudnum,cloudrad,cloudsig,cloudprof,inlinetemps,press,inwavenum,linelist,cia,ciatemps,use_disort,clphot,ophot,make_cf,do_bff,bff)

     # This bit trims off the unused portion of the array
     trimspec = np.zeros((2,nwave),dtype='d')
     trimspec[:,:] = tmpoutspec[:,:nwave]
     shiftspec = np.empty_like(trimspec)
     shiftspec[0,:] =  trimspec[0,:] + dlam
     shiftspec[1,:] =  trimspec[1,:]



     benchspec = np.loadtxt('No_cloud_800K_model_benchmark_SPEC.dat',skiprows=3,unpack=True)
     wno = 1e4 / shiftspec[0,:]
     modspec = spex_non_uniform(benchspec[0,:],wno,shiftspec[1,:])


     difference_spectrum = modspec / benchspec[1,:]
     print('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-')
     print('------------------------------------------------------------')
     print(' Test for T dwarf case, forward model only, No clouds')
     print('------------------------------------------------------')
     print('mean value of modelspectrum / T benchmark = '+str(np.mean(difference_spectrum)))
     print('std deviation of modelspectrum / T benchmark = '+str(np.std(difference_spectrum)))

     percent_change = np.mean(abs(modspec - benchspec[1,:])/benchspec[1,:])
     
     if (percent_change < 0.01):
          print("less than 1percent difference with T dwarf regime benchmark")
          print('-------------------------------------------------------')
          print('   ')
          return True
     else:
          print("greater than 1 percent difference T dwarf regime with benchmark")
          print('-------------------------------------------------------')
          print('   ')
          return False


def MieClouds_Ldwarf(xpath):

     # Now the wavelength range
     w1 = 1.0
     w2 = 15.

     # FWHM of data in microns(WE DON'T USE THIS FOR SPEX DATA. SET TO 0.0)
     # -1 = spex + AKARI + spitzer
     # -2 = spex + spitzer
     fwhm = -3

     # DISTANCE (in parsecs)
     dist = 11.35

     # How many patches & clouds do we want??
     # Must be at least 1 of each, but can turn off cloud below
     npatches = 1
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
     cloudnum[:,0] = 8
     cloudtype[:,0] = 1
     
     cloudnum[:,1] = 9
     cloudtype[:,1] = 2
     
     # second patch turn off top cloud
     #cloudnum[1,0] = 5
     #cloudtype[1,0] = 1
     
     # Are we assuming chemical equilibrium, or similarly precomputed gas abundances?
     # Or are we retrieving VMRs (0)
     chemeq = 0
     
     # Are we doing H- bound-free, free-free continuum opacities?
     # (Is the profile going above 3000K in the photosphere?)
     do_bff = 1
     
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
     
     # Now list the gases.
     # If Na is after K, at the end of the list, alkalis will be tied
     # together at Asplund solar ratio. See Line at al (2015)
     # Else if K is after Na, they'll be separate
     
     gaslist = ['h2o','co','co2','ch4','tio','vo','crh','feh','k','na']
     
     
     # some switches for alternative cross sections
     # Use Mike's Alkalis?
     malk = 0
     # Use Mike's CH4?
     mch4 = 0

     # Are we using DISORT for radiative transfer?
     # (HINT: Not in this century)
     use_disort = 0 
     
     # use the fudge factor?
     do_fudge = 1

     prof = np.full(13,100.)
     if (proftype == 9):
          modP,modT = np.loadtxt(pfile,skiprows=1,usecols=(1,2),unpack=True)
          tfit = InterpolatedUnivariateSpline(np.log10(modP),modT,k=1)
          prof = tfit(logcoarsePress)



     # Now we'll get the opacity files into an array
     inlinetemps,inwavenum,linelist,gasnum,nwave = testkit.get_opacities(gaslist,w1,w2,press,xpath,malk)
      
     # Get the cia bits
     tmpcia, ciatemps = ciamod.read_cia("CIA_DS_aug_2015.dat",inwavenum)
     cia = np.asfortranarray(np.empty((4,ciatemps.size,nwave)),dtype='float32')
     cia[:,:,:] = tmpcia[:,:,:nwave] 
     ciatemps = np.asfortranarray(ciatemps, dtype='float32')
     
     # grab BFF and Chemical grids
     bff_raw,ceTgrid,metscale,coscale,gases_myP = testkit.sort_bff_and_CE(chemeq,"chem_eq_tables_P3K.pic",press,gaslist)
     
     
     
     # now set theta... this is from 2M2224_cryESlabRustDeck
     theta = np.array([-3.55278369e+00,-2.83012757e+00,-4.31062021e+00,-4.95190596e+00,-9.77059307e+00,-8.85409603e+00,-8.41539430e+00,-7.85745521e+00,-6.57250890e+00,5.46814691e+00,2.68655361e-20,1.07209671e+00,1.11922607e+00,2.83604013e-03,-3.16119701e+01,-3.32775232e+01,-3.46762823e+01,5.42024548e+00,-2.76574938e+00,4.38059949e-01,-5.73919866e-01,8.58329576e-02,8.72374998e-01,4.39392990e+00,-1.96757779e+00,6.24967679e-02,3.45025551e-01,6.78307874e-02,7.56891116e-02,1.71616709e+00,4.88646433e+03])


     benchspec = np.loadtxt('Mie_cloud_1800K_model_benchmark_SPEC.dat',skiprows=3,unpack=True)
      
     runargs = gases_myP,chemeq,dist, cloudtype,do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,benchspec,proftype,do_fudge, prof,do_bff,bff_raw,ceTgrid,metscale,coscale

     gnostics = 0
     shiftspec, clphotspec, ophotspec,cfunc = testkit.modelspec(theta,runargs,gnostics)

     modspec = brewtools.proc_spec(shiftspec,theta,fwhm,chemeq,gasnum,benchspec) 


      

     difference_spectrum = modspec / benchspec[1,:]
     print('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-')
     print('------------------------------------------------------------')
     print(' Test for L dwarf case, plotting via theta and testkit')
     print(' theta taken from 2M2224 case with Mie clouds')
     print(' crystalline enstatite slab + rust deck')
     print(' optical data as at August 2019 in benchmark spectrum')
     print('------------------------------------------------------')
     print('mean value of modelspectrum / L benchmark = '+str(np.mean(difference_spectrum)))

     print('std deviation of modelspectrum / L benchmark = '+str(np.std(difference_spectrum)))

     percent_change = np.mean(abs(modspec - benchspec[1,:])/benchspec[1,:])
     
     if (percent_change < 0.01):
          print("less than 1percent difference with L dwarf regime benchmark")
          print('-------------------------------------------------------')
          return True
     else:
          print("greater than 1 percent difference with L dwarf regime benchmark")
          print('-------------------------------------------------------')

          return False

     
