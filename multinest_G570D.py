#!/usr/bin/env python
"""MCMC Retrieval Setup Template"""
import os 
import utils
import numpy as np
import retrieval_run
import settings


__author__ = "Fei Wang"
__copyright__ = "Copyright 2024 - Fei Wang"
__credits__ = ["Fei Wang", "Ben Burningham"]
__license__ = "GPL"
__version__ = "0.2"  
__maintainer__ = ""
__email__ = ""
__status__ = "Development"


fwhm=0
wavelength_range=[1.0,2.5]
ndata=1
# wavpoints=None

chemeq=0
gaslist =  ['h2o','co','co2','ch4','nh3','h2s','k','na']
gastype_list=['U','U','U','U','U','U','U','U']

ptype=1
do_clouds=0
npatches=1
cloudname = ['power law cloud slab']  
cloudpacth_index=[[1]] 
particle_dis=['hansan']


# cloudname = []  
# cloudpacth_index=[] 


# particle_dis=['hansan','log_normal']
# cloudname = ['power law cloud slab']  
# do_bff=1

do_fudge=1
# samplemode='mcmc'
samplemode='multinest'

instrument_instance = utils.Instrument(fwhm,wavelength_range,ndata)
re_params = utils.Retrieval_params(samplemode,chemeq,gaslist,gastype_list,fwhm,do_fudge,ptype,do_clouds,npatches,cloudname,cloudpacth_index,particle_dis)
model_config_instance = utils.ModelConfig(samplemode)
io_config_instance = utils.IOConfig()


io_config_instance.outdir="/beegfs/car/fei/lsr1835/test/"
io_config_instance.runname='multinest_test_G570D_clear'
io_config_instance.update_dictionary()


model_config_instance.dist= 5.84
model_config_instance.xlist ='gaslistRox.dat'
model_config_instance.do_bff=0
model_config_instance.malk=1
model_config_instance.ch4=0
model_config_instance.const_efficiency_mode=True
model_config_instance.sampling_efficiency=0.3
model_config_instance.multimodal = False
model_config_instance.log_zero= -1e90
model_config_instance.importance_nested_sampling= False
model_config_instance.evidence_tolerance=0.1
model_config_instance.update_dictionary()


obspec = np.asfortranarray(np.loadtxt("G570D_2MHcalib.dat",dtype='d',unpack='true')) # G570D_2MassJcalib.dat
args_instance = utils.ArgsGen(re_params,model_config_instance,instrument_instance,obspec)
settings.init(args_instance)
retrieval_run.brewster_reterieval_run(re_params,model_config_instance,io_config_instance)

