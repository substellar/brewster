#!/usr/bin/env python


""" Module of bits to define Instrument, ModelConfig, IOConfig and Retrieval_params class"""

import numpy as np
import ciamod
import os
import gc
import sys
import pickle
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
import TPmod
from collections import namedtuple

__author__ = "Fei Wang"
__copyright__ = "Copyright 2024 - Fei Wang"
__credits__ = ["Fei Wang", "Ben Burningham"]
__license__ = "GPL"
__version__ = "0.2"  
__maintainer__ = ""
__email__ = ""
__status__ = "Development"


class Instrument:
    """
    A class to represent an instrument object.
    
    Parameters
    ----------
    data : str
        Telescope observing modes: 'IFS', 'Imaging'
    R_file : float
        Spectral resolution  file
    wavelength_range : float
        Maximum wavelength (um)
    ndata : float
        number of instruments
    
    Methods
    -------
    instrument_dic_gen()
    """
    
    def __init__(self, fwhm=None, wavelength_range=None, ndata=None,R_file=None,wavpoints=None):
        self.fwhm  = fwhm 
        self.wavelength_range = wavelength_range
        self.ndata = ndata
        self.wavpoints = wavpoints
        self.R_file = R_file

        self.dictionary = self.instrument_dic_gen()

#     @classmethod
    def instrument_dic_gen(self):
        """
        Initialize telescope object using current parameters.
        """
        return {
            'params': {
                'fwhm': self.fwhm,
                'wavelength_range': self.wavelength_range,
                'ndata': self.ndata,
                'wavpoints': self.wavpoints,
                'R_file': self.R_file
            }
        }

    def __str__(self):
        string = 'Instrument: \n------------\n' +\
            '- fwhm : ' + "%s" % (self.fwhm) + '\n' +\
            '- wavelength_range : ' + "%s" % (self.wavelength_range) + '\n' +\
            '- ndata : ' + "%s" % (self.ndata) + ' \n'  +\
            '- wavpoints : ' + "%s" % (self.wavpoints) + ' \n' +\
            '- R_file : ' + "%s" % (self.R_file) + ' \n'
        return string



class ModelConfig:
    """
    A class to represent a model configuration object.

    Parameters
    ----------
    samplemode: str
        The samplemode to use for parameter generation ('mcmc' or 'pymultinest')
    use_disort : int, optional
        Use DISORT flag (default: 0)
    do_fudge : int, optional
        Do fudge factor flag (default: 0)
    malk : int, optional
        Alkali metals flag (default: 0)
    mch4 : int, optional
        Methane flag (default: 0)
    do_bff : int, optional
        Do BFF flag (default: 1)
    fresh : int, optional
        Fresh start flag (default: 0)
    xpath : str, optional
        Path to line lists (default: "../Linelists/")
    xlist : str, optional
        Line list file (default: "data/gaslistRox.dat")
    dist : float, optional
        Distance parameter (default: None)
    pfile : str, optional
        Pt file (default: "LSR1835_eqpt.dat")

    Methods
    -------
    model_normal_dic_gen()
        Initialize normal model configuration dictionary.
    model_sampling_dic_gen()
        Initialize model configuration dictionary based on the specified method.
    update_dictionary()
        Update the model configuration dictionary with the current attributes.
    """

    def __init__(self, samplemode, do_fudge, use_disort=0, malk=0, mch4=0, do_bff=1, fresh=0, xpath="../Linelists/", xlist="data/gaslistRox.dat", dist=None, pfile="data/LSR1835_eqpt.dat"):
        self.samplemode = samplemode
        self.use_disort = use_disort
        self.do_fudge = do_fudge
        self.malk = malk
        self.mch4 = mch4
        self.do_bff = do_bff
        self.fresh = fresh
        self.xpath = xpath
        self.xlist = xlist
        
        self.dist = dist
        self.dist_err = 0
        self.pfile = pfile
        
        # Default MCMC parameters
        self.nwalkers = None
        self.nburn = 10000
        self.niter = 30000

        # Default MultiNest parameters
        self.LogLikelihood = None
        self.Prior = None
        self.ndim = None
        self.nparam= None
        self.n_clustering_params = None
        self.wrapped_params = None
        self.importance_nested_sampling = True
        self.multimodal = True
        self.const_efficiency_mode = False
        self.n_live_points = 400
        self.evidence_tolerance = 0.5
        self.sampling_efficiency = 0.8
        self.n_iter_before_update = 100
        self.null_log_evidence = -1.e90
        self.max_modes = 100
        self.mode_tolerance = -1.e90
        self.outputfiles_basename = ''
        self.seed = -1
        self.verbose = True
        self.resume = True
        self.context = 0
        self.write_output = True
        self.log_zero = -1.e100
        self.max_iter = 0
        self.init_MPI = False
        self.dump_callback = None

        self.dictionary_normal = self.model_normal_dic_gen()
        self.dictionary_sampling = self.model_sampling_dic_gen()

    def model_normal_dic_gen(self):
        """
        Initialize normal model configuration dictionary.
        """
        return {
            'params': {
                'use_disort': self.use_disort,
                'do_fudge': self.do_fudge,
                'malk': self.malk,
                'mch4': self.mch4,
                'do_bff': self.do_bff,
                'fresh': self.fresh,
                'xpath': self.xpath,
                'xlist': self.xlist,
                'dist': self.dist,
                'pfile': self.pfile
            }
        }

    def model_sampling_dic_gen(self):
        """
        Initialize model configuration dictionary based on the specified method.
        """
        if self.samplemode.lower() == 'mcmc':
            return {
                'params': {
                    'ndim': self.ndim,
                    'nwalkers': self.nwalkers,
                    'nburn': self.nburn,
                    'niter': self.niter
                }
            }
        elif self.samplemode.lower() == 'multinest':
            return {
                'LogLikelihood': self.LogLikelihood,
                'Prior': self.Prior,
                'ndim': self.ndim,
                'nparam': self.nparam,
                'n_clustering_params': self.n_clustering_params,
                'wrapped_params': self.wrapped_params,
                'importance_nested_sampling': self.importance_nested_sampling,
                'multimodal': self.multimodal,
                'const_efficiency_mode': self.const_efficiency_mode,
                'n_live_points': self.n_live_points,
                'evidence_tolerance': self.evidence_tolerance,
                'sampling_efficiency': self.sampling_efficiency,
                'n_iter_before_update': self.n_iter_before_update,
                'null_log_evidence': self.null_log_evidence,
                'max_modes': self.max_modes,
                'mode_tolerance': self.mode_tolerance,
                'outputfiles_basename': self.outputfiles_basename,
                'seed': self.seed,
                'verbose': self.verbose,
                'resume': self.resume,
                'context': self.context,
                'write_output': self.write_output,
                'log_zero': self.log_zero,
                'max_iter': self.max_iter,
                'init_MPI': self.init_MPI,
                'dump_callback': self.dump_callback
            }
        else:
            raise ValueError("Unsupported samplemode. Please choose 'mcmc' or 'multinest'.")

    def update_dictionary(self):
        """Update the model configuration dictionary with the current attributes."""
        self.dictionary_normal = self.model_normal_dic_gen()
        self.dictionary_sampling = self.model_sampling_dic_gen()

    def __str__(self):
        base_info = (
            f"ModelConfig: \n------------\n"
            f"- use_disort : {self.use_disort}\n"
            f"- do_fudge : {self.do_fudge}\n"
            f"- malk : {self.malk}\n"
            f"- mch4 : {self.mch4}\n"
            f"- do_bff : {self.do_bff}\n"
            f"- fresh : {self.fresh}\n"
            f"- xpath : {self.xpath}\n"
            f"- xlist : {self.xlist}\n"
            f"- dist : {self.dist}\n"
            f"- pfile : {self.pfile}\n"
            f"\n"
        )
        if self.samplemode.lower() == 'mcmc':
            samplemode_info = (
                f"(MCMC): \n------------\n"
                f"- ndim : {self.ndim}\n"
                f"- nwalkers : {self.nwalkers}\n"
                f"- nburn : {self.nburn}\n"
                f"- niter : {self.niter}\n"
            )
        elif self.samplemode.lower() == 'multinest':
            samplemode_info = (
                f"PyMultiNest: \n----------------------\n"
                f"- LogLikelihood: {self.LogLikelihood} \n"
                f"- Prior: {self.Prior} \n"
                f"- ndim: {self.ndim} \n"
                f"- nparam: {self.nparam}\n"
                f"- n_clustering_params: {self.n_clustering_params}\n"
                f"- wrapped_params: {self.wrapped_params}\n"
                f"- importance_nested_sampling: {self.importance_nested_sampling}\n"
                f"- multimodal: {self.multimodal}\n"
                f"- const_efficiency_mode: {self.const_efficiency_mode}\n"
                f"- n_live_points: {self.n_live_points}\n"
                f"- evidence_tolerance: {self.evidence_tolerance}\n"
                f"- sampling_efficiency: {self.sampling_efficiency}\n"
                f"- n_iter_before_update: {self.n_iter_before_update}\n"
                f"- null_log_evidence: {self.null_log_evidence}\n"
                f"- max_modes: {self.max_modes}\n"
                f"- mode_tolerance: {self.mode_tolerance}\n"
                f"- outputfiles_basename: {self.outputfiles_basename}\n"
                f"- seed: {self.seed}\n"
                f"- verbose: {self.verbose}\n"
                f"- resume: {self.resume}\n"
                f"- context: {self.context}\n"
                f"- write_output: {self.write_output}\n"
                f"- log_zero: {self.log_zero}\n"
                f"- max_iter: {self.max_iter}\n"
                f"- init_MPI: {self.init_MPI}\n"
                f"- dump_callback: {self.dump_callback}\n"
            )
        return base_info + samplemode_info


    
class IOConfig:
    """
    A class to represent an IO configuration object.
    
    Parameters
    ----------
    runname : str, optional
        Name of the run (default: 'retrieval')
    outdir : str, optional
        Output directory (default: None)
    finalout : str, optional
        Final output file (default: None)
    chaindump : int, optional
        Chain dump flag (default: 0)
    picdump : int, optional
        Picdump flag (default: 0)
    statfile : str, optional
        status file (default: None)
    rfile : str, optional
        runtimes file (default: None)
    runtest : int, optional
        Run test flag (default: 1)
    make_arg_pickle : int, optional
        Make argument pickle flag (default: 2)
    
    Methods
    -------
    IO_config_dic_gen()
        Initialize IO configuration dictionary.
    update_dictionary()
        Update the model configuration dictionary with the current attributes.
    """

    def __init__(self, runname='retrieval', outdir=None, finalout=None, chaindump=0, picdump=0, statfile=None, rfile=None, runtest=1, make_arg_pickle=2):
        self.runname = runname
        self.outdir = outdir
        self.finalout = finalout
        self.chaindump = chaindump
        self.picdump = picdump
        self.statfile = statfile
        self.rfile = rfile
        self.runtest = runtest
        self.make_arg_pickle = make_arg_pickle
        self.dictionary = self.IO_config_dic_gen()
        self.update_dictionary()

    def IO_config_dic_gen(self):
        """
        Initialize IO configuration dictionary.
        """
        return {
            'params': {
                'runname': self.runname,
                'outdir': self.outdir,
                'finalout': self.runname + ".pk1",
                'chaindump': self.runname + "_last_nc.pic",
                'picdump': self.runname + "_snapshot.pic",
                'statfile': "status_ball_" + self.runname + ".txt",
                'rfile':  "runtimes_" + self.runname + ".dat",
                'runtest': self.runtest,
                'make_arg_pickle': self.make_arg_pickle
            }
        }
    
    def update_dictionary(self):
        """
        Update the model configuration dictionary with the current attributes.
        """
        config_dict = self.IO_config_dic_gen()
        self.runname = config_dict['params']['runname']
        self.finalout = config_dict['params']['finalout']
        self.chaindump = config_dict['params']['chaindump']
        self.picdump = config_dict['params']['picdump']
        self.statfile = config_dict['params']['statfile']
        self.rfile = config_dict['params']['rfile']
        self.dictionary = config_dict

    def __str__(self):
        params = self.dictionary['params']
        return (
            f"IOConfig: \n------------\n"
            f"- runname : {params['runname']}\n"
            f"- outdir : {params['outdir'] if params['outdir'] is not None else 'N/A'}\n"
            f"- finalout : {params['finalout']}\n"
            f"- chaindump : {params['chaindump']}\n"
            f"- picdump : {params['picdump']}\n"
            f"- statfile : {params['statfile']}\n"
            f"- rfile : {params['rfile']}\n"
            f"- runtest : {params['runtest']}\n"
            f"- make_arg_pickle : {params['make_arg_pickle']}\n"
        )




cloud_dic = {
    "Al2O3": 1,
    "Fe": 2,
    "H2O": 3,
    "Mg2SiO4": 4,
    "MgSiO3": 5,
    "Mg2SiO4rich": 6,
    "NH3": 7,
    "MgSiO3Cry": 8,
    "Fe2O3_WS15": 9,
    "Mg2SiO4Cry": 10,
    "MgSiO3Cry_s2": 11,
    "SiO2_WS15": 12,
    "MgAl2O4_WS15": 13,
    "TiO2_WS15": 14,
    "Na2S": 15,
    "CaTiO3": 16,
    "MnS-newIR": 17,
    "KCl": 18,
    "ZnS": 19,
    "SiO": 20,
    "mixto": 99
}



def dp_customized_distribution(x):
    return np.abs(0.1 * np.random.randn(x))


def hansan_b_customized_distribution(x):
    return np.abs(0.2+0.05 * np.random.randn(x))


# def tolerance_parameter_customized_distribution(x):

#     args_instance=settings.runargs

#     return np.log10((np.random.rand(x)* (max(args_instance.obspec[2,:]**2)*(0.1 - 0.01))) + (0.01*min(args_instance.obspec[2,10::3]**2))) 




class Retrieval_params:
    """
    A class to represent retrieval parameters for atmospheric retrieval.

    Parameters
    ----------
    chemeq : int
        A flag to determine if chemical equilibrium is considered (1) or not (0).
    gaslist : list of str
        List of gas names.
    gastype_list : list of str
        List of gas types, corresponding to the gas names in gaslist.
    fwhm : float, optional
        Full width at half maximum of the spectral lines. 
    do_fudge : int, optional
        Flag indicating whether to apply tolerance_parameter to the data. 
    ndata : int
        number of instruments
    ptype : int
        Type of pressure-temperature profile.
    do_clouds : int, optional
        Flag indicating whether to include cloud parameters in the retrieval. 
        - 1: Include clouds.
    npatches : int, optional
        Number of patches for cloud distribution. 
    cloudname : list of str
        List of cloud names.
    cloudpatch_index : list of int, optional
        Indexes for `cloudname` corresponding to which cloud patches 
    particle_dis : str, optional
        Distribution type for particles in the cloud. Default is None.
        E.g., 'log_normal', 'hansan', etc.
        only used when include Mie cloud.
    
    Methods
    -------
    gas_dic_gen(gasname, gastype):
        Generates dictionary for gas parameters.
    pt_dic_gen(ptype):
        Generates dictionary for pressure-temperature profile parameters.
    cloud_dic_gen(cloudname, cloudspecies=None, docloud=None):
        Generates dictionary for cloud parameters.
    refinement_params_dic_gen(ndata):
        Generates dictionary for additional parameters.
    gas_allparams_gen(chemeq, gaslist, gastype_list):
        Generates dictionary for all gas parameters.
    cloud_allparams_gen(cloudname):
        Generates dictionary for all cloud parameters.
    retrieval_para_dic_gen(chemeq, gaslist, gastype_list, ndata, ptype, cloudname):
        Generates the complete dictionary for retrieval parameters.
    __str__():
        String representation of the class instance.
    """
    
    def __init__(self, samplemode,chemeq=None, gaslist=None, gastype_list=None,fwhm=None,do_fudge=1,ptype=None,do_clouds=1,npatches=None,cloudname=None,cloudpatch_index=None,particle_dis=None):
        self.samplemode = samplemode
        self.chemeq = chemeq
        self.gaslist = gaslist
        self.gastype_list = gastype_list
        self.fwhm = fwhm
        self.do_fudge = do_fudge
        self.ptype = ptype
        self.do_clouds = do_clouds
        self.cloudname = cloudname
        self.cloudpatch_index=cloudpatch_index
        self.particle_dis=particle_dis
        
        self.dictionary = self.retrieval_para_dic_gen(chemeq, gaslist, gastype_list,fwhm,do_fudge, ptype,do_clouds,npatches,cloudname,cloudpatch_index,particle_dis)
        
        
        
    def gas_dic_gen(self,gasname,gastype):
        
        dictionary = {}
        if gastype=='U':

            dictionary[gasname]={
                'gastype':gastype,
                'params':{'log_abund':
                           {'initialization':None,
                            'distribution':['normal',-4.0,0.5],
                            'range':None,
                            'prior':None}
                         }}

        elif gastype=='N':

            dictionary[gasname]={
                'gastype':gastype,
                'params':{'log_abund':
                           {'initialization':None,
                            'distribution':['normal',-4.0,0.5],
                            'range':None,
                            'prior': None},

                           "p_ref": 
                            {'initialization':None,
                              'distribution':['normal',-1,0.2],
                              'range':None,
                              'prior': None},

                           "alpha":
                            {'initialization':None,
                             'distribution':['uniform',0,1],
                             'range':None,
                             'prior': None}    
                           }}
        elif gastype=='H':
            dictionary[gasname]={
                'gastype':gastype,
                'params':{'log_abund':
                           {'initialization':None,
                            'distribution':['normal',-4.0,0.5],
                            'range':None,
                            'prior': None},
                            
                           "p_ref": 
                            {'initialization':None,
                              'distribution':['normal',-1,0.2],
                              'range':None,
                              'prior': None}

                           }}

        return dictionary
            
    
    
    
    def pt_dic_gen(self,ptype):
        dictionary = {}

        if ptype==1:
            dictionary={
                'ptype':ptype,
                'params':{'gamma':
                           {'initialization':None,
                            'distribution':['normal',50,1],
                             'range':None,
                            'prior':None}}}
            
            for i in range(13):
                dictionary['params']["T_%d" % (i+1)] = {
                    'initialization': None,
                    'distribution': ['normal',500, 50],
                    'range':[0,5000],
                    'prior': None
                }
            


        elif ptype==2:

            dictionary={
                'ptype':ptype,
                'params':{'alpha1':
                           {'initialization':None,
                            'distribution':['normal',0.2,0.1],
                            'range':None,
                            'prior':None},

                          'alpha2':
                           {'initialization':None,
                            'distribution':['normal',0.18,0.05],
                            'range':None,
                            'prior':None},

                          'logP1':
                           {'initialization':None,
                            'distribution':['normal',-1,0.2],
                            'range':None,
                            'prior':None},

                          'logP3':
                           {'initialization':None,
                            'distribution':['normal',2,0.2],
                            'range':None,
                            'prior':None},

                          'T3':
                           {'initialization':None,
                            'distribution':['normal',3600,500],
                            'range':None,
                            'prior':None}
                         }}
        elif ptype==3:

            dictionary={
                'ptype':ptype,
                'params':{'alpha1':
                           {'initialization':None,
                            'distribution':['normal',0.2,0.1],
                            'range':None,
                            'prior':None},

                          'alpha2':
                           {'initialization':None,
                            'distribution':['normal',0.18,0.05],
                            'range':None,
                            'prior':None},

                          'logP1':
                           {'initialization':None,
                            'distribution':['normal',-1,0.2],
                            'range':None,
                            'prior':None},

                          'logP2':
                           {'initialization':None,
                            'distribution':['normal',0.1,0.2],
                            'range':None,
                            'prior':None},

                          'logP3':
                           {'initialization':None,
                            'distribution':['normal',2,0.2],
                            'range':None,
                            'prior':None},

                          'T3':
                           {'initialization':None,
                            'distribution':['normal',3600,500],
                            'range':None,
                            'prior':None}
                         }}

        elif ptype==7 :

            dictionary={
                'ptype':ptype,
                'params':{'Tint':
                           {'initialization':None,
                            'distribution':['normal',1200,200],
                            'range':None,
                            'prior':None},

                          'alpha':
                           {'initialization':None,
                            'distribution':['uniform',1,2],
                            'range':None,
                            'prior':None},

                          'lndelta':
                           {'initialization':None,
                            'distribution':['normal',0,1],
                            'range':None,
                            'prior':None},

                          'T1':
                           {'initialization':None,
                            'distribution':['normal',1200,200],
                            'range':None,
                            'prior':None},

                          'T2':
                           {'initialization':None,
                            'distribution':['normal',1200,200],
                            'range':None,
                            'prior':None},

                          'T3':
                           {'initialization':None,
                            'distribution':['normal',1200,200],
                            'range':None,
                            'prior':None}
                         }}

        elif ptype==77:

            dictionary={
                'ptype':ptype,
                'params':{'gamma':
                           {'initialization':None,
                            'distribution':['normal',50,1],
                            'range':None,
                            'prior':None},

                          'Tint':
                           {'initialization':None,
                            'distribution':['normal',1200,200],
                            'range':None,
                            'prior':None},

                          'alpha':
                           {'initialization':None,
                            'distribution':['uniform',1,2],
                            'range':None,
                            'prior':None},
                            
                          'lndelta':
                           {'initialization':None,
                            'distribution':['normal',0,1],
                            'range':None,
                            'prior':None},

                          'T1':
                           {'initialization':None,
                            'distribution':['normal',1200,200],
                            'range':None,
                            'prior':None},

                          'T2':
                           {'initialization':None,
                            'distribution':['normal',1200,200],
                            'range':None,
                            'prior':None},

                          'T3':
                           {'initialization':None,
                            'distribution':['normal',1200,200],
                            'range':None,
                            'prior':None}
                         }}

        elif ptype==9:
                dictionary={'ptype':ptype,
                'params':{}}


        return dictionary



    def cloud_dic_gen(self,do_clouds,cloudname,particle_dis=None):
        dictionary = {}
        if do_clouds==0:
            return {}
        
        else:

            if cloudname=='grey cloud deck':

                dictionary["patch"]={
                    'cloudnum': 99,
                    'cloudtype':2,
                    'params':{'logp_gcd':
                               {'initialization':None,
                                'distribution':['normal',1,0.1],
                                'range':None,
                                'prior':None},
                              'dp_gcd':
                               {'initialization':None,
                                'distribution':['customized',dp_customized_distribution],  #lambda x: np.abs(0.1 * np.random.randn(x))
                                'range':None,
                                'prior':None},
                              'omega_gcd':
                              {'initialization':None,
                               'distribution':['uniform',0,1],
                               'range':None,
                                'prior':None}
                             }}

            elif cloudname=='grey cloud slab':

                dictionary["patch"]={
                    'cloudnum': 99,
                    'cloudtype':1,

                    'params':{'tau_gcs':
                               {'initialization':None,
                                'distribution':['normal',10,1],
                                'range':None,
                                'prior':None},
                              'logp_gcs':
                               {'initialization':None,
                                'distribution':['normal',-0.2,0.1],
                                'range':None,
                                'prior':None},
                              'dp_gcs':
                              {'initialization':None,
                               'distribution':['customized',dp_customized_distribution],  #lambda x: np.abs(0.5+0.01* np.random.randn(x))
                               'range':None,
                                'prior':None},
                              'omega_gcs':
                              {'initialization':None,
                               'distribution':['uniform',0,1],
                               'range':None,
                                'prior':None}
                             }}

            elif cloudname=='powerlaw cloud deck':

                dictionary["patch"]={
                    'cloudnum': 89,
                    'cloudtype':2,

                    'params':{'logp_pcd':
                               {'initialization':None,
                                'distribution':['normal',-0.2,0.1],
                                'range':None,
                                'prior':None},
                              'dp_pcd':
                               {'initialization':None,
                                'distribution':['customized',dp_customized_distribution], #lambda x: np.abs(0.1 * np.random.randn(x))
                                'range':None,
                                'prior':None},
                              'omega_pcd':
                              {'initialization':None,
                               'distribution':['uniform',0,1],
                               'range':None,
                                'prior':None},
                              'alpha_pcd':
                              {'initialization':None,
                               'distribution':['normal',0,1],
                               'range':None,
                                'prior':None}
                             }}

            elif 'Mie scattering cloud deck' in cloudname:
                
                cloudspecies=cloudname.split('--')[1].strip()
                cloudnum=cloud_dic.get(cloudspecies,None)


                if particle_dis=="hansan":
                    dictionary["patch"]={
                        'cloudnum': cloudnum,
                        'cloudtype':2,
                        "particle_dis":"hansan",
                        'params':{'logp_mcd_%s'%cloudspecies:
                                    {'initialization':None,
                                    'distribution':['normal',1,0.1],
                                    'range':None,
                                    'prior':None},
                                   'dp_mcd_%s'%cloudspecies:
                                    {'initialization':None,
                                    'distribution':['customized',dp_customized_distribution], #lambda x: np.abs(0.1 * np.random.randn(x))
                                    'range':None, 
                                    'prior':None},
                                    'hansan_a_mcd_%s'%cloudspecies:
                                    {'initialization':None,
                                    'distribution':['normal',-1.4,0.1],
                                    'range':None,
                                    'prior':None},
                                    'hansan_b_mcd_%s'%cloudspecies:
                                    {'initialization':None,
                                    'distribution':['customized',hansan_b_customized_distribution], #lambda x: np.abs(0.2+0.05 * np.random.randn(x))
                                    'range':None,
                                    'prior':None}
                                        }}
                                
                if particle_dis=="log_normal":
                    dictionary["patch"]={
                        'cloudnum': cloudnum,
                        'cloudtype':1,
                        "particle_dis":"log_normal",
                        'params':{'logp_mcd_%s'%cloudspecies:
                                    {'initialization':None,
                                    'distribution':['normal',1,0.1],
                                    'range':None,
                                    'prior':None},
                                    'dp_mcd_%s'%cloudspecies:
                                    {'initialization':None,
                                    'distribution':['customized',dp_customized_distribution], #lambda x: np.abs(0.1 * np.random.randn(x))
                                    'range':None,
                                    'prior':None},
                                    'mu_mcd_%s'%cloudspecies:
                                    {'initialization':None,
                                    'distribution':['normal',0,1],
                                    'range':None,
                                    'prior':None},
                                    'sigma_mcd_%s'%cloudspecies:
                                    {'initialization':None,
                                    'distribution':['normal',0,1],
                                    'range':None,
                                    'prior':None}
                                        }}

            elif cloudname=='power law cloud slab':

                dictionary["patch"]={
                    'cloudnum': 89,
                    'cloudtype':1,

                    'params':{'tau_pcs':
                               {'initialization':None,
                                'distribution':['normal',10,1],
                                'range':None,
                                'prior':None},
                              'logp_pcs':
                               {'initialization':None,
                                'distribution':['normal',-0.2,0.5],
                                'range':None,
                                'prior':None},
                              'dp_pcs':
                              {'initialization':None,
                               'distribution':['customized', dp_customized_distribution], #lambda x: np.abs(0.1 * np.random.randn(x))
                               'range':None,
                                'prior':None},
                              'omega_pcs':
                              {'initialization':None,
                               'distribution':['uniform',0,1],
                               'range':None,
                                'prior':None},
                              'alpha_pcs':
                              {'initialization':None,
                               'distribution':['normal',0,1],
                               'range':None,
                                'prior':None}
                             }}


            elif 'Mie scattering cloud slab' in cloudname:
                cloudspecies=cloudname.split('--')[1].strip()
                cloudnum=cloud_dic.get(cloudspecies,None)

                if particle_dis=="hansan":
                    dictionary["patch"]={
                        'cloudnum': cloudnum,
                        'cloudtype':1,
                        'particle_dis':"hansan",
                        'params':{'tau_mcs_%s'%cloudspecies:
                                   {'initialization':None,
                                    'distribution':['normal',10,1],
                                    'range':None,
                                    'prior':None},
                                  'logp_mcs_%s'%cloudspecies:
                                   {'initialization':None,
                                    'distribution':['normal',-0.2,0.5],
                                    'range':None,
                                    'prior':None},
                                  'dp_mcs_%s'%cloudspecies:
                                  {'initialization':None,
                                   'distribution':['customized',dp_customized_distribution], # lambda x: np.abs(0.1 * np.random.randn(x))
                                   'range':None,
                                    'prior':None},
                                  'hansan_a_mcs_%s'%cloudspecies:
                                  {'initialization':None,
                                   'distribution':['normal',-1.4,0.1],
                                   'range':None,
                                    'prior':None},
                                  'hansan_b_mcs_%s'%cloudspecies:
                                  {'initialization':None,
                                   'distribution':['customized',hansan_b_customized_distribution], #lambda x: np.abs(0.2+0.05 * np.random.randn(x))
                                   'range':None,
                                   'prior':None}
                                     }}

                if particle_dis=="log_normal":
                    dictionary["patch"]={
                        'cloudnum': cloudnum,
                        'cloudtype':1,
                        'particle_dis':"log_normal",
                        'params':{'tau_mcs_%s'%cloudspecies:
                                   {'initialization':None,
                                    'distribution':['normal',10,1],
                                    'range':None,
                                    'prior':None},
                                  'logp_mcs_%s'%cloudspecies:
                                   {'initialization':None,
                                    'distribution':['normal',-0.2,0.5],
                                    'range':None,
                                    'prior':None},
                                  'dp_mcs_%s'%cloudspecies:
                                  {'initialization':None,
                                   'distribution':['customized',dp_customized_distribution], #lambda x: np.abs(0.1 * np.random.randn(x))
                                   'range':None,
                                    'prior':None},
                                  'mu_mcs_%s'%cloudspecies:
                                  {'initialization':None,
                                   'distribution':['normal',0,1],
                                   'range':None,
                                    'prior':None},
                                  'sigma_mcs_%s'%cloudspecies:
                                  {'initialization':None,
                                   'distribution':['normal',0,1],
                                   'range':None,
                                    'prior':None}
                                     }}

            elif cloudname=='clear':
                dictionary["patch"]={'params':{}}


            return dictionary
    
    
    def refinement_params_dic_gen(self):
        """
        Initialize model configuration dictionary based on the specified method.
    
        """
        
        
        dictionary = {}


        if self.samplemode.lower() == 'mcmc':
               dictionary['params'] = {
                    'logg': {
                        'initialization': None,
                        'distribution': ['normal', 4.5, 0.1],
                        'prior': None
                    },
                    'r2d2': {
                        'initialization': None,
                        'distribution': ['normal', 0, 1],
                        'prior': None
                    },
                    'scale1': {
                        'initialization': None,
                        'distribution': ['normal', 0, 0.001],
                        'prior': None
                    },
                    'scale2': {
                        'initialization': None,
                        'distribution': ['normal', 0, 0.001],
                        'prior': None
                    },
                    'dlambda': {
                        'initialization': None,
                        'distribution': ['normal', 0, 0.001],
                        'prior': None
                    }
                }
        elif self.samplemode.lower() == 'multinest':
                    dictionary['params'] = {
                    'M': {
                        'initialization': None,
                        'distribution': ['normal', 4.5, 0.1],
                        'range':None,
                        'prior': None
                    },
                    'R': {
                        'initialization': None,
                        'distribution': ['normal', 0, 1],
                        'range':None,
                        'prior': None
                    },
                    'scale1': {
                        'initialization': None,
                        'distribution': ['normal', 0, 0.001],
                        'range':None,
                        'prior': None
                    },
                    'scale2': {
                        'initialization': None,
                        'distribution': ['normal', 0, 0.001],
                        'range':None,
                        'prior': None
                    },
                    'dlambda': {
                        'initialization': None,
                        'distribution': ['normal', 0, 0.001],
                        'range':None,
                        'prior': None
                    }
                }
        else:
            raise ValueError("Unsupported samplemode. Please choose 'mcmc' or 'multinest'.")
            
        # Remove 'scale1' and 'scale2' if fwhm condition is not met
    
        if self.fwhm >=0 or self.fwhm in [-5,-6] and self.do_fudge==1:

            del dictionary['params']['scale1']
            del dictionary['params']['scale2']
            ndata=1

        if self.fwhm in [-2] and self.do_fudge==1:
            del dictionary['params']['scale1']
            ndata=2

        if self.fwhm in [-1, -3, -4] and self.do_fudge==1:
            ndata=3

        # Add tolerance parameters after 'dlambda'
        if self.do_fudge==1:
            for i in range(ndata):
                dictionary['params']["tolerance_parameter_%d" % (i+1)] = {
                    'initialization': None,
                    'distribution': ['customized', 0],
                    'range':None,
                    'prior': None
                }
        return dictionary


    
    
    def gas_allparams_gen(self,chemeq,gaslist,gastype_list):
        gas_dic= {}
        if chemeq==1:
             gas_dic= {'params':
                       {'mh':
                           {'initialization':None,
                            'distribution':['normal',-1,3],
                            'range':None,
                            'prior':None},
                       'co':
                           {'initialization':None,
                            'distribution':['normal',0.25,2.25],
                            'range':None,
                            'prior':None}
                      }}

        elif chemeq==0:
            gasnum=len(gaslist)
            gaslist_lower = [gas.lower() for gas in gaslist]

            if gaslist_lower[-2:] == ['k', 'na']:
                gasnum=gasnum-1
                gaslist=list(gaslist[0:-2])+["K_Na"]
                
            elif gaslist_lower[-3:] == ['k','na','cs']:
                gasnum=gasnum-2
                gaslist=list(gaslist[0:-3])+['K_Na_Cs']
                
                
            for i in range(0,gasnum):
                dic=self.gas_dic_gen(gaslist[i],gastype_list[i])
                gas_dic.update(dic)  

        return gas_dic
    
    
    
    def cloud_allparams_gen(self,do_clouds,npatches,cloudname,cloudpatch_index,particle_dis):
        cloud_dic = {}

        if do_clouds > 0:
            if npatches > 1:
                cloud_dic["fcld"] = {
                    'initialization': None,
                    'distribution': ['uniform', 0, 1],
                    'range':None,
                    'prior': None
                }

            for i in range(npatches):
                for j in range(len(cloudpatch_index)):
                    if i + 1 in cloudpatch_index[j]:
                        dic = self.cloud_dic_gen(do_clouds, cloudname[j], particle_dis[j])
                        patch_key = f"patch {i + 1}"

                        if patch_key not in cloud_dic:
                            cloud_dic[patch_key] = {}

                        cloud_dic[patch_key][cloudname[j]] = dic["patch"]

                # if i+1 > len(cloudpatch_index):
                #     patch_key = f"patch {i + 1}"
                #     cloud_dic[patch_key] = {}

            return cloud_dic
        else:
            return {}


    def added_params(self,param):
        added_dic = {
            f"{param}": {
                'initialization': None,
                'distribution': ['normal', 0, 1],
                'range': None,
                'prior': None}}
        self.dictionary["added_params"]=added_dic

    
    
    def retrieval_para_dic_gen(self,chemeq,gaslist,gastype_list,fwhm,do_fudge,ptype,do_clouds,npatches,cloudname,cloudpacth_index,particle_dis):
        retrieval_param={}
        gas_dic=self.gas_allparams_gen(chemeq,gaslist,gastype_list)
        refinement_dic=self.refinement_params_dic_gen()
        pt_dic=self.pt_dic_gen(ptype)
        cloud_dic=self.cloud_allparams_gen(do_clouds,npatches,cloudname,cloudpacth_index,particle_dis) 
        retrieval_param["gas"]=gas_dic
        retrieval_param["refinement_params"]=refinement_dic
        retrieval_param["pt"]=pt_dic
        retrieval_param["cloud"]=cloud_dic

        return retrieval_param
    
    
    def update_dictionary(self):
        """Update the model configuration dictionary with the current attributes."""
        self.dictionary.update(self.dictionary)

            
    def __str__(self):
        s = []
        
        # Retrieve cloud keys
        cloud_keys = self.dictionary.get('cloud', {}).keys()
        ncloud = len(cloud_keys) - 1 if len(cloud_keys) > 1 else 1

        # Loop through cloud patches and retrieve their parameter keys
        for i in range(ncloud):
            patch_key = f'patch {i + 1}'
            cloud_types = self.dictionary['cloud'].get(patch_key, {}).keys()
            patch_params = [param for cloud_type in cloud_types for param in self.dictionary['cloud'][patch_key][cloud_type]['params'].keys()]
            s.append(f'  -- {patch_key}: {patch_params}\n')
        s = ''.join(s)

        # Retrieve gas, refinement_params, and pt parameter keys
        gas_keys = self.dictionary.get('gas', {}).keys()
        refinement_keys = self.dictionary.get('refinement_params', {}).get('params', {}).keys()
        pt_keys = self.dictionary.get('pt', {}).get('params', {}).keys()

        add_params_keys  = self.dictionary.get("added_params", {}).keys()

        # Format and return the complete string
        return (
            'retrieval_param: \n------------\n'
            f'- gas : {list(gas_keys)}\n'
            f'- refinement_params : {list(refinement_keys)}\n'
            f'- pt : {list(pt_keys)}\n'
            f'- cloud : {list(cloud_keys)}\n'
            f'{s}'
            f'- added_params: {list(add_params_keys)}\n'
        )
            

def get_all_parametres(dic):

    gaslist = list(dic['gas'].keys())
    gastype_values = [info['gastype'] for key, info in dic['gas'].items() if 'gastype' in info]
    
    gas=[]
    gas_values=[]
    
    for i in range(len(gaslist)):
        gas.append(gaslist[i])
        gas_values.append(dic['gas'][gaslist[i]]['params']['log_abund']['initialization'])
        if  gastype_values[i]=='N':
            gas.append("p_ref_%s"%gaslist[i])
            gas.append("alpha_%s"%gaslist[i])
            gas_values.append(dic['gas'][gaslist[i]]['params']['p_ref']['initialization'])
            gas_values.append(dic['gas'][gaslist[i]]['params']['alpha']['initialization'])    

        elif gastype_values[i]=='H':
            gas.append("p_ref_%s"%gaslist[i])
            gas_values.append(dic['gas'][gaslist[i]]['params']['p_ref']['initialization'])
       

            

    refinement_params = list(dic['refinement_params']['params'].keys())
    refinement_params_values = [dic['refinement_params']['params'][p]['initialization'] for p in refinement_params]

    pt = list(dic['pt']['params'].keys())
    pt_values = [dic['pt']['params'][p]['initialization'] for p in pt]
    
    
    
    cloud = []
    cloud_values = []
    if 'cloud' in dic:
        if 'fcld' in dic['cloud']:
            cloud.append('fcld')
            cloud_values.append(dic['cloud']['fcld']['initialization'])
        for patch_key in dic['cloud']:
            if patch_key.startswith('patch'):
                for cloud_key in dic['cloud'][patch_key]:
                    for param_key in dic['cloud'][patch_key][cloud_key]['params']:
                        cloud.append(f"{param_key}")
                        cloud_values.append(dic['cloud'][patch_key][cloud_key]['params'][param_key]['initialization'])

    all_params = gas + refinement_params + pt + cloud
    all_params_values = gas_values + refinement_params_values + pt_values + cloud_values

    if 'added_params' in dic.keys():
        added_dic = list(dic['added_params'].keys())
        added_values = [dic['added_params'][p]['initialization'] for p in added_dic]

        all_params += added_dic
        all_params_values += added_values


    unique_params = {}
    for param, value in zip(all_params, all_params_values):
        if param not in unique_params:
            unique_params[param] = value
            
    return list(unique_params.keys()), list(unique_params.values())



def update_dictionary(dic, params_instance):
    # Update gas parameters
    
    gastype_values = [info['gastype'] for key, info in dic['gas'].items() if 'gastype' in info]
    gaslist=list(dic['gas'].keys())
    for i in range(len(gaslist)):
        dic['gas'][gaslist[i]]['params']['log_abund']['initialization'] = getattr(params_instance, gaslist[i])
        if gastype_values[i]=='N':
            dic['gas'][gaslist[i]]['params']['p_ref']['initialization'] = getattr(params_instance, "p_ref_%s"%gaslist[i])
            dic['gas'][gaslist[i]]['params']['alpha']['initialization'] = getattr(params_instance, "alpha_%s"%gaslist[i])

        if gastype_values[i]=='H':
            dic['gas'][gaslist[i]]['params']['p_ref']['initialization'] = getattr(params_instance, "p_ref_%s"%gaslist[i])

    
    # Update refinement parameters
    for param in dic['refinement_params']['params'].keys():
        dic['refinement_params']['params'][param]['initialization'] = getattr(params_instance, param)
    
    # Update pt parameters
    for param in dic['pt']['params'].keys():
        dic['pt']['params'][param]['initialization'] = getattr(params_instance, param)


    # Update added parameters 
    if 'added_params' in dic.keys():
        for param in dic['added_params'].keys():
            dic['added_params'][param]['initialization'] = getattr(params_instance, param)
    
    # Update cloud parameters
    
    if 'cloud' in dic and dic['cloud']:
        if 'patch 1' in dic['cloud']:
            for cloud_type in dic['cloud']['patch 1'].keys():
                for param in dic['cloud']['patch 1'][cloud_type]['params'].keys():
                    if hasattr(params_instance, param):
                        dic['cloud']['patch 1'][cloud_type]['params'][param]['initialization'] = getattr(params_instance, param)

        if 'fcld' in dic['cloud'] and hasattr(params_instance, 'fcld'):
            dic['cloud']['fcld']['initialization'] = getattr(params_instance, 'fcld')

        if 'patch 2' in dic['cloud']:
            for cloud_type in dic['cloud']['patch 2'].keys():
                for param in dic['cloud']['patch 2'][cloud_type]['params'].keys():
                    if hasattr(params_instance, param):
                        dic['cloud']['patch 2'][cloud_type]['params'][param]['initialization'] = getattr(params_instance, param)
    
    return dic





def get_distribution_values(dic):
    distribution_values = []

    def recurse(d):
        if isinstance(d, dict):
            for key, value in d.items():
                if key == 'distribution':
                    distribution_values.append(value)
                else:
                    recurse(value)
        elif isinstance(d, list):
            for item in d:
                recurse(item)

    recurse(dic)
    return distribution_values


def MC_P0_gen(updated_dic,model_config_instance,args_instance):

    nwalkers=model_config_instance.nwalkers
    ndim=model_config_instance.ndim
    p0 = np.empty([nwalkers,ndim])
    
    
    all_distributions=get_distribution_values(updated_dic)
    
    for i in range(model_config_instance.ndim):
        if all_distributions[i][0]=='normal':
            mu,sigma=all_distributions[i][1:]
            p0[:,i]=mu+sigma*np.random.randn(nwalkers).reshape(nwalkers)
            
        elif  all_distributions[i][0]=='uniform':
            pmin,pmax=all_distributions[i][1:]
            p0[:,i]= np.random.uniform(pmin, pmax, nwalkers).reshape(nwalkers)

        elif  all_distributions[i][0]=='customized':
            f = all_distributions[i][1]
            p0[:, i] = f(nwalkers).reshape(nwalkers)

    if args_instance.proftype==1:

        all_params,all_params_values =get_all_parametres(updated_dic) 
        params_master = namedtuple('params',all_params)
        params_instance = params_master(*all_params_values)
        T_1_index=params_instance._fields.index('T_1')
        T_13_index=params_instance._fields.index('T_13')
        BTprof = np.loadtxt("data/BTtemp800_45_13.dat")

        for i in range(0, 13):  # 13 layer points ====> Total: 13 + 13 (gases+) +no cloud = 26
            p0[:,T_1_index+i] = (BTprof[i] - 200.) + (150. * np.random.randn(nwalkers).reshape(nwalkers))

        for i in range(0, nwalkers):
            while True:
                Tcheck = TPmod.set_prof(args_instance.proftype,args_instance.coarsePress,args_instance.press,p0[i, T_1_index:T_13_index+1])
                if min(Tcheck) > 1.0:
                    break
                else:
                    for i in range(0,13):
                        p0[:,T_1_index+i] = BTprof[i] + (50. * np.random.randn(nwalkers).reshape(nwalkers))

    return p0

    
    
def cloud_para_gen(dic):


    # Determine number of patches and clouds
    npatches = 1
    nclouds = 1


    if 'cloud' not in dic or not dic['cloud']:
        return np.zeros([npatches, nclouds], dtype='i'), np.zeros([npatches, nclouds], dtype='i'), np.zeros([npatches], dtype='i')

    if 'patch 1' in dic['cloud']:
        # Determine if there is more than one patch
        if 'patch 2' in dic['cloud'] and dic['cloud']['patch 2']:
            npatches = 2
            
        cloudnum_list=[]
        for cloud_type, cloud_info in dic['cloud']['patch 1'].items():
            if 'cloudnum' in cloud_info:
                cloudnum_list.append(cloud_info['cloudnum'])
            
        # if len(cloudnum_list) ==2: #and len(set(cloudnum_list)) == len(cloudnum_list):
        #     nclouds = 2
        nclouds= len(cloudnum_list)

    # Initialize arrays
    do_clouds = np.zeros([npatches], dtype='i')
    cloudnum = np.zeros([npatches, nclouds], dtype='i')
    cloudtype = np.zeros([npatches, nclouds], dtype='i')

    # Populate arrays
    patch_index = 0
    for patch_key in dic['cloud'].keys():
        if patch_key.startswith('patch 1'):            
            for i in range(nclouds):
                cloudkey=list(dic['cloud'][patch_key].keys())
                cloud_info = dic['cloud'][patch_key][cloudkey[i]]
                if 'cloudnum' in cloud_info and 'cloudtype' in cloud_info:
                    cloudnum[patch_index, i] = cloud_info['cloudnum']
                    cloudtype[patch_index, i] = cloud_info['cloudtype']
                    do_clouds[patch_index] = 1

            
    if npatches==2:
        cloudkey=list(dic['cloud']['patch 2'].keys())
        cloud_info = dic['cloud']['patch 2'][cloudkey[0]]
        if 'cloudnum' in cloud_info and 'cloudtype' in cloud_info:
            cloudnum[1, 0] = cloud_info['cloudnum']
            cloudnum[1, 1] = cloud_info['cloudnum']
            cloudtype[1, 1] = cloud_info['cloudtype']
            do_clouds[1] = 1
        
                    
    return cloudnum, cloudtype, do_clouds



# def args_gen(re_params,model,instrument,obspec):

#     # set up pressure grids in log(bar) cos its intuitive
#     logcoarsePress = np.arange(-4.0, 2.5, 0.53)
#     logfinePress = np.arange(-4.0, 2.4, 0.1)

#     # but forward model wants pressure in bar
#     coarsePress = pow(10,logcoarsePress)
#     press = pow(10,logfinePress)        
    
#     dist=model.dist
#     use_disort=model.use_disort
#     xpath=model.xpath
#     xlist=model.xlist
#     malk=model.malk
#     do_fudge=model.do_fudge
#     pfile=model.pfile
#     do_bff=model.do_bff
#     gaslist=list(re_params.dictionary['gas'].keys())
#     gaslist_lower = [gas.lower() for gas in gaslist]

#     chemeq=re_params.chemeq
    
#     if gaslist_lower[-1] == 'k_na':
#         gaslist=list(gaslist[0:-1])+['k', 'na']

#     elif gaslist_lower[-1] == 'k_na_cs':
#         gaslist=list(gaslist[0:-1])+['k','na','cs']

    
#     fwhm=instrument.fwhm
#     w1,w2=instrument.wavelength_range[:]
#     proftype=re_params.ptype

#     cloudnum, cloudtype, do_clouds=cloud_para_gen(re_params.dictionary)
    
#     prof = np.full(13, 100.)
#     if proftype == 9:
#         modP, modT = np.loadtxt(pfile, skiprows=1, usecols=(1, 2), unpack=True)
#         tfit = InterpolatedUnivariateSpline(np.log10(modP), modT, k=1)
#         prof = tfit(logcoarsePress)
        
        
#     inlinetemps,inwavenum,linelist,gasnum,nwave = get_opacities(gaslist,w1,w2,press,xpath,xlist,malk)
#     tmpcia, ciatemps = ciamod.read_cia("CIA_DS_aug_2015.dat",inwavenum)
#     cia = np.asfortranarray(np.empty((4,ciatemps.size,nwave)),dtype='float32')
#     cia[:,:,:] = tmpcia[:,:,:nwave] 
#     ciatemps = np.asfortranarray(ciatemps, dtype='float32')
    
    
#     # grab BFF and Chemical grids
#     bff_raw,ceTgrid,metscale,coscale,gases_myP = sort_bff_and_CE(chemeq,"chem_eq_tables_P3K.pic",press,gaslist)
    
#     return gases_myP,chemeq,dist,cloudtype,do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge, prof,do_bff,bff_raw,ceTgrid,metscale,coscale




def get_opacities(gaslist,w1,w2,press,xpath='../Linelists',xlist='data/gaslistR10K.dat',malk=0):
    # Now we'll get the opacity files into an array
    ngas = len(gaslist)

    totgas = 0
    gasdata = []
    with open(xlist) as fa:
        for line_aa in fa.readlines():
            if len(line_aa) == 0:
                break
            totgas = totgas +1 
            line_aa = line_aa.strip()
            gasdata.append(line_aa.split())


    list1 = []
    for i in range(0,ngas):
        for j in range(0,totgas):
            if (gasdata[j][1].lower() == gaslist[i].lower()):
                list1.append(gasdata[j])

    if (malk == 1):
        for i in range (0,ngas):
            list1[i] = [w.replace('K_', 'K_Mike_') for w in list1[i]]
            list1[i] = [w.replace('Na_', 'Na_Mike_') for w in list1[i]]

    if (malk == 2):
        for i in range (0,ngas):
            list1[i] = [w.replace('K_', 'K_2021_') for w in list1[i]]
            list1[i] = [w.replace('Na_', 'Na_2021_') for w in list1[i]]



    lists = [xpath+i[3] for i in list1[0:ngas]]
    gasmass = np.asfortranarray(np.array([i[2] for i in list1[0:ngas]],dtype='float32'))


    # get the basic framework from water list
    rawwavenum, inpress, inlinetemps, inlinelist = pickle.load(open(lists[0], "rb"))

    wn1 = 10000. / w2
    wn2 = 10000. / w1
    inwavenum = np.asfortranarray(rawwavenum[np.where(np.logical_not(np.logical_or(rawwavenum[:] > wn2, rawwavenum[:] < wn1)))],dtype='float64')
    ntemps = inlinetemps.size
    npress= press.size
    nwave = inwavenum.size
    r1 = np.amin(np.where(np.logical_not(np.logical_or(rawwavenum[:] > wn2, rawwavenum[:] < wn1))))
    r2 = np.amax(np.where(np.logical_not(np.logical_or(rawwavenum[:] > wn2, rawwavenum[:] < wn1))))

    # Here we are interpolating the linelist onto our fine pressure scale.
    # pickles have linelist as 4th entry....
    linelist = (np.zeros([ngas,npress,ntemps,nwave],order='F')).astype('float64', order='F')
    for gas in range (0,ngas):
        inlinelist= pickle.load(open(lists[gas], "rb" ) )[3]
        for i in range (0,ntemps):
            for j in range (r1,r2+1):
                pfit = interp1d(np.log10(inpress),np.log10(inlinelist[:,i,j]))
                linelist[gas,:,i,(j-r1)] = np.asfortranarray(pfit(np.log10(press)))
    linelist[np.isnan(linelist)] = -50.0
    
    # convert gaslist into fortran array of ascii strings for fortran code 
    gasnames = np.empty((len(gaslist), 10), dtype='c')
    for i in range(0,len(gaslist)):
        gasnames[i,0:len(gaslist[i])] = gaslist[i]

    gasnames = np.asfortranarray(gasnames,dtype='c')

    return inlinetemps,inwavenum,linelist,gasnames,gasmass,nwave




def sort_bff_and_CE(chemeq,ce_table,press,gaslist):

    # Sort out the BFF opacity stuff and chemical equilibrium tables:
    metscale,coscale,Tgrid,Pgrid,gasnames,abunds = pickle.load( open(ce_table, "rb" ) )
    nabpress = Pgrid.size
    nabtemp = Tgrid.size
    nabgas = abunds.shape[4]
    nmet = metscale.size
    nco = coscale.size
    nlayers = press.size
    ngas = len(gaslist)


    bff_raw = np.zeros([nabtemp,nlayers,3])
    gases_myP = np.zeros([nmet,nco,nabtemp,nlayers,ngas+3])
    gases = np.zeros([nmet,nco,nabtemp,nabpress,ngas+3])

    if (chemeq == 0):
        # Just want the ion fractions for solar metallicity in this case
        ab_myP = np.empty([nabtemp,nlayers,nabgas])
        i1 = np.where(metscale == 0.0)  #high metallicity, solar metallicity=0.0
        i2 = np.where(coscale == 1.0)
        for gas in range (0,nabgas):
            for i in range (0,nabtemp):
                pfit = InterpolatedUnivariateSpline(Pgrid,np.log10(abunds[i1[0],i2[0],i,:,gas]),k=1)
                ab_myP[i,:,gas] = pfit(np.log10(press))
                
        bff_raw[:,:,0] = ab_myP[:,:,0]
        bff_raw[:,:,1] = ab_myP[:,:,2]
        bff_raw[:,:,2] = ab_myP[:,:,4]
        #bff_raw[:,:,2] = 10**-50

    else:

        # In this case we need the rows for the gases we're doing and ion fractions
        gases[:,:,:,:,0] = abunds[:,:,:,:,0] # 'e-'
        gases[:,:,:,:,1] = abunds[:,:,:,:,2] # 'H'
        gases[:,:,:,:,2] = abunds[:,:,:,:,4] # 'H-'
        nmatch = 0

        if (gaslist[len(gaslist)-1] == 'h_mins'):

            for i in range(0,ngas-1):
                for j in range(0,nabgas):
                    if (gasnames[j].lower() == gaslist[i].lower()):
                        gases[:,:,:,:,i+3] = abunds[:,:,:,:,j]
                        nmatch = nmatch + 1
            # print(nmatch)
            if (nmatch != ngas-1):
                print("you've requested a gas that isn't in the Vischer table. Please check and try again.")
                sys.exit()
            
            gases[:,:,:,:,-1] = abunds[:,:,:,:,0]*0+1e-90

            for i in range(0,nmet):
                for j in range(0,nco):
                    for k in range(0,ngas+3):
                        for l in range(0,nabtemp):
                            pfit = InterpolatedUnivariateSpline(Pgrid,np.log10(gases[i,j,l,:,k]),k=1)
                            gases_myP[i,j,l,:,k] = pfit(np.log10(press))

        else:
            # In this case we need the rows for the gases we're doing and ion fractions
            # gases[:,:,:,:,0] = abunds[:,:,:,:,0] # 'e-'
            # gases[:,:,:,:,1] = abunds[:,:,:,:,2] # 'H'
            # gases[:,:,:,:,2] = abunds[:,:,:,:,4] # 'H-'
            # nmatch = 0


            for i in range(0,ngas):
                for j in range(0,nabgas):
                    if (gasnames[j].lower() == gaslist[i].lower()):
                        gases[:,:,:,:,i+3] = abunds[:,:,:,:,j]
                        nmatch = nmatch + 1
            if (nmatch != ngas):
                print("you've requested a gas that isn't in the Vischer table. Please check and try again.")
                sys.exit()

            for i in range(0,nmet):
                for j in range(0,nco):
                    for k in range(0,ngas+3):
                        for l in range(0,nabtemp):
                            pfit = InterpolatedUnivariateSpline(Pgrid,np.log10(gases[i,j,l,:,k]),k=1)
                            gases_myP[i,j,l,:,k] = pfit(np.log10(press))


    return bff_raw,Tgrid,metscale,coscale,gases_myP






class ArgsGen:
    """
    A class to generate and manage forward model arguments.
    
    Parameters
    ----------
    re_params : object
        Object containing retrieval parameters.
    model : object
        Object containing the model configuration.
    instrument : Instrument
        Instrument class containing instrument specifications.
    obspec : object
        Observation spectra object.
    
    Attributes
    ----------
    gases_myP : list
        Processed gas list at pressures.
    chemeq : bool
        Chemical equilibrium flag.
    dist : float
        Distance to the object.
    gasnames: c
        List of gasnames in fortran compatible character arrays
    gasmass: float
        List of molecular masses for gases
    cloudtype : str
        Type of cloud parameterization used.
    do_clouds : bool
        Whether clouds are considered.
    cloudnum : int
        Number of cloud layers.
    inlinetemps : np.array
        Inline temperatures from opacities.
    coarsePress : np.array
        Coarse pressure grid.
    press : np.array
        Fine pressure grid.
    inwavenum : np.array
        Input wave numbers.
    linelist : list
        Line list for the opacities.
    cia : np.array
        CIA data array.
    ciatemps : np.array
        CIA temperature grid.
    use_disort : bool
        DISORT flag for radiative transfer.
    fwhm : float
        Full width at half maximum for the instrument.
    obspec : object
        Observation spectra.
    proftype : int
        Profile type.
    do_fudge : bool
        Fudge factor flag.
    prof : np.array
        Atmospheric profile.
    do_bff : bool
        BFF flag.
    bff_raw : np.array
        Raw BFF grid.
    ceTgrid : np.array
        Chemical equilibrium temperature grid.
    metscale : np.array
        Metallicity scaling grid.
    coscale : np.array
        Carbon-to-oxygen ratio scaling grid.
    
    Methods
    -------
    generate()
        Generate the required model arguments.
    """

    def __init__(self, re_params, model, instrument, obspec):
        self.re_params = re_params
        self.model = model
        self.instrument = instrument
        self.obspec = obspec

        # Generate all necessary model arguments on initialization
        self.generate()

    def generate(self):
        # Set up pressure grids in log(bar)
        logcoarsePress = np.arange(-4.0, 2.5, 0.53)
        logfinePress = np.arange(-4.0, 2.4, 0.1)
        
        # Pressure in bar
        self.coarsePress = pow(10, logcoarsePress)
        self.press = pow(10, logfinePress)
        
        # Retrieve model parameters
        self.dist = self.model.dist
        self.dist_err = self.model.dist_err
        self.use_disort = self.model.use_disort
        self.xpath = self.model.xpath
        self.xlist = self.model.xlist
        self.malk = self.model.malk
        self.do_fudge = self.model.do_fudge
        self.pfile = self.model.pfile
        self.do_bff = self.model.do_bff
        self.chemeq = self.re_params.chemeq
        
        # Process gas list
        self.gaslist = list(self.re_params.dictionary['gas'].keys())
        gaslist_lower = [gas.lower() for gas in self.gaslist]
        
        if gaslist_lower[-1] == 'k_na':
            self.gaslist = list(self.gaslist[:-1]) + ['k', 'na']
        elif gaslist_lower[-1] == 'k_na_cs':
            self.gaslist = list(self.gaslist[:-1]) + ['k', 'na', 'cs']
        
        # Retrieve instrument parameters
        self.fwhm = self.instrument.fwhm
        self.w1, self.w2 = self.instrument.wavelength_range
        
        # Profile type and cloud parameters
        self.proftype = self.re_params.ptype
        self.cloudnum, self.cloudtype, self.do_clouds = cloud_para_gen(self.re_params.dictionary)
        
        # Generate temperature profile
        self.prof = np.full(13, 100.0)
        if self.proftype == 9:
            logmodP,modT = np.loadtxt(self.pfile,skiprows=0,usecols=(0,1),unpack=True)
            tfit = InterpolatedUnivariateSpline(logmodP,modT,k=1)
            self.prof = tfit(np.log10(self.coarsePress))
        
        # Get opacities, CIA data
        self.inlinetemps, self.inwavenum, self.linelist,self.gasnames,self.gasmass, self.nwave = get_opacities(
            self.gaslist, self.w1, self.w2, self.press, self.xpath, self.xlist, self.malk)

        self.tmpcia, self.ciatemps = ciamod.read_cia("data/CIA_DS_aug_2015.dat", self.inwavenum)
        self.cia = np.asfortranarray(np.empty((4, self.ciatemps.size, self.nwave)), dtype='float32')
        self.cia[:, :, :] = self.tmpcia[:, :, :self.nwave]
        self.ciatemps = np.asfortranarray(self.ciatemps, dtype='float32')
        
        # BFF and Chemical grids
        self.bff_raw, self.ceTgrid, self.metscale, self.coscale, self.gases_myP = sort_bff_and_CE(
            self.chemeq, "data/chem_eq_tables_P3K.pic", self.press, self.gaslist)

        
    def __str__(self):
        return f"""
        ArgsGen Model Parameters:
        -------------------------
        Distance: {self.dist} (Error: {self.dist_err})
        Chemical Equilibrium: {self.chemeq}
        FWHM: {self.fwhm}
        Wavelength Range: {self.w1} - {self.w2}
        DISORT Flag: {self.use_disort}
        do_fudge {self.do_fudge}
        do_bff: {self.do_bff}
        Profile Type: {self.proftype}
        Gas List: {self.gaslist}
        do_clouds: {self.do_clouds}
        Number of Clouds: {self.cloudnum}
        Cloud Type: {self.cloudtype} 
        Metallicity Scale: {self.metscale}
        C/O Ratio Scale: {self.coscale}
        Coarse Pressure Grid: {self.coarsePress}
        """

        # Temperature Profile: {self.prof}
        # Inlined Temperatures: {self.inlinetemps}
        # Line List: {self.linelist}
        # CIA Data: {self.cia}
        # CIA Temperatures: {self.ciatemps}
        # BFF Raw Grid: {self.bff_raw}
        # Fine Pressure Grid: {self.press}
        # Chemical Grid T: {self.ceTgrid}
