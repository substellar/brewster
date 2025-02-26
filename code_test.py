# This script tests the forward model in the context of cloudless T8 dwarf
# parameters a drawn from a retrieval on G570D
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mgimg
import matplotlib.colors as colors
import scipy as sp
import numpy as np
import emcee
import corner
import pickle as pickle
import forwardmodel
import ciamod
import TPmod
import brewtools
from astropy.convolution import convolve, convolve_fft
from astropy.convolution import Gaussian1DKernel
from scipy.io.idl import readsav
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from bensconv import prism_non_uniform
from bensconv import conv_uniform_R
from bensconv import conv_uniform_FWHM
from collections import namedtuple
import utils
import settings
import gas_nonuniform
import test_module





def NoCloud_Tdwarf(xpath,xlist):
    
     fwhm=0
     wavelength_range=[1,2.5]
     ndata=1

     ##gas
     chemeq=0
     gaslist = ['h2o','ch4','co','co2','nh3','h2s','K','Na']
     gastype_list=['U','U','U','U','U','U','U','U']
     ptype=9

     ## clouds

     do_clouds=0
     npatches=1

     cloudname = ['clear']
     cloudpacth_index=[[1]]
     particle_dis=[]

     # ModelConfig:
     do_fudge=0
     samplemode='mcmc'

     instrument_instance = utils.Instrument(fwhm,wavelength_range,ndata)
     re_params = utils.Retrieval_params(samplemode,chemeq,gaslist,gastype_list,fwhm,do_fudge,ptype,do_clouds,npatches,cloudname,cloudpacth_index,particle_dis)
     model_config_instance = utils.ModelConfig(samplemode,do_fudge)
     io_config_instance = utils.IOConfig()



     model_config_instance.do_bff=0
     model_config_instance.malk=0
     model_config_instance.pfile='data/test_data/G570D_model_benchmark_PROFILE.dat'
     model_config_instance.xlist=xlist #'gaslistR10K.dat'
     model_config_instance.xpath=xpath
     model_config_instance.update_dictionary()


     obspec = []#np.asfortranarray(np.loadtxt("LSR1835_data_realcalib_new_trimmed.dat",dtype='d',unpack='true')) # obs is not actually using

     args_instance = utils.ArgsGen(re_params,model_config_instance,instrument_instance,obspec)
     settings.init(args_instance)
     args_instance=settings.runargs

     all_params,all_params_values =utils.get_all_parametres(re_params.dictionary) 
     params_master = namedtuple('params',all_params)
     print(all_params)
     print(model_config_instance.do_fudge)
     theta=[-3.27,-3.36,-7.27,-8.28,-4.73,-8.71,-5.36]+[4.89]+[1.50901046e-19]+[0.00258329]
     params_instance = params_master(*theta)

     gnostics=0
     shiftspec, cloud_phot_press,other_phot_press,cfunc=test_module.modelspec(params_instance,re_params,args_instance,gnostics)

     modspec = np.array([shiftspec[0,::-1],shiftspec[1,::-1]])
     benchspec = np.loadtxt('data/test_data/No_cloud_800K_model_benchmark_SPEC.dat',skiprows=3,unpack=True)
     outspec = prism_non_uniform(benchspec,modspec,3.3)

     difference_spectrum = outspec / benchspec[1,:]
     print('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-')
     print('------------------------------------------------------------')
     print(' Test for T dwarf case, forward model only, No clouds')
     print('------------------------------------------------------')
     print('mean value of modelspectrum / T benchmark = '+str(np.mean(difference_spectrum)))
     print('std deviation of modelspectrum / T benchmark = '+str(np.std(difference_spectrum)))

     percent_change = np.mean(abs(outspec - benchspec[1,:])/benchspec[1,:])
     
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


def MieClouds_Ldwarf(xpath,xlist):
     fwhm=-3
     wavelength_range=[1,15]
     ndata=1

     #retrieval_params
     ##gas
     chemeq=0
     gaslist = ['h2o','co','co2','ch4','tio','vo','crh','feh','k','na']
     gastype_list=['U','U','U','U','U','U','U','U','U','U']
     ptype=2

     ## clouds
     do_clouds=1
     npatches=1
     cloudname = ['Mie scattering cloud slab--MgSiO3Cry','Mie scattering cloud deck--Fe2O3_WS15']
     cloudpacth_index=[[1],[1]]
     particle_dis=['hansan','hansan']
     # ModelConfig:

     do_fudge = 1
     samplemode='mcmc'
     instrument_instance = utils.Instrument(fwhm,wavelength_range,ndata)
     re_params = utils.Retrieval_params(samplemode,chemeq,gaslist,gastype_list,fwhm,do_fudge,ptype,do_clouds,npatches,cloudname,cloudpacth_index,particle_dis)
     model_config_instance = utils.ModelConfig(samplemode,do_fudge)
     io_config_instance = utils.IOConfig()



     model_config_instance.do_bff=1
     model_config_instance.malk=0
     model_config_instance.pfile="t1700g1000f3.dat"
     model_config_instance.xlist=xlist #'gaslistR10K.dat'
     model_config_instance.xpath=xpath
     model_config_instance.dist=11.35
     model_config_instance.update_dictionary()

     # obspec = np.asfortranarray(np.loadtxt("LSR1835_data_realcalib_new_trimmed.dat",dtype='d',unpack='true'))

     obspec= np.loadtxt('data/test_data/Mie_cloud_1800K_model_benchmark_SPEC.dat',skiprows=3,unpack=True)
     args_instance = utils.ArgsGen(re_params,model_config_instance,instrument_instance,obspec)
     settings.init(args_instance)
     args_instance=settings.runargs

     all_params,all_params_values =utils.get_all_parametres(re_params.dictionary) 
     params_master = namedtuple('params',all_params)
     theta =[-3.55278369e+00,-2.83012757e+00,-4.31062021e+00,-4.95190596e+00,-9.77059307e+00,-8.85409603e+00,-8.41539430e+00,-7.85745521e+00,-6.57250890e+00,5.46814691e+00,2.68655361e-20,1.07209671e+00,1.11922607e+00,2.83604013e-03,-3.16119701e+01,-3.32775232e+01,-3.46762823e+01,5.42024548e+00,-2.76574938e+00,4.38059949e-01,-5.73919866e-01,8.58329576e-02,8.72374998e-01,4.39392990e+00,-1.96757779e+00,6.24967679e-02,3.45025551e-01,6.78307874e-02,7.56891116e-02,1.71616709e+00,4.88646433e+03]
     theta_master=theta[0:17]+theta[26:]+theta[17:26]
     params_instance = params_master(*theta_master)
     gnostics=0
     shiftspec, cloud_phot_press,other_phot_press,cfunc=test_module.modelspec(params_instance,re_params,args_instance,gnostics)

     obspec=args_instance.obspec
     modspec = np.array([shiftspec[0,::-1],shiftspec[1,::-1]])
     mr1 = np.where(modspec[0,:] < 2.5)
     or1  = np.where(obspec[0,:] < 2.5)
     spec1 = prism_non_uniform(obspec[:,or1],modspec,3.3)

     # Mike Cushing supplied L band R = 425 
     # dispersion constant across order 0.0097um
     # R = 425
     R = 425
     mr2 = np.where(np.logical_and(modspec[0,:] > 2.5,modspec[0,:] < 5.0))
     or2 = np.where(np.logical_and(obspec[0,:] > 2.5,obspec[0,:] < 5.0))
     spec2 = params_instance.scale1 * conv_uniform_R(obspec[:,or2],modspec,R)

     # Spitzer IRS
     # R roughly constant within orders, and orders both appear to
     # have R ~ 100
     R = 100.0
     mr3 = np.where(modspec[0,:] > 5.0)
     or3 = np.where(obspec[0,:] > 5.0)
     spec3 = params_instance .scale2 * conv_uniform_R(obspec[:,or3],modspec,R)
     outspec =  np.array(np.concatenate((spec1,spec2,spec3),axis=0))


     difference_spectrum = outspec / obspec[1,:]

     print('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-')
     print('------------------------------------------------------------')
     print(' Test for L dwarf case, plotting via theta and testmodule')
     print(' theta taken from 2M2224 case with Mie clouds')
     print(' crystalline enstatite slab + rust deck')
     print(' optical data as at February 2020 in benchmark spectrum')
     print('------------------------------------------------------')
     print('mean value of modelspectrum / L benchmark = '+str(np.mean(difference_spectrum)))
     print('std deviation of modelspectrum / L benchmark = '+str(np.std(difference_spectrum)))
     percent_change = np.mean(abs(outspec - obspec[1,:])/obspec[1,:])
     
     if (percent_change < 0.01):
          print("less than 1percent difference with L dwarf regime benchmark")
          print('-------------------------------------------------------')
          return True
     else:
          print("greater than 1 percent difference with L dwarf regime benchmark")
          print('-------------------------------------------------------')

          return False

     
