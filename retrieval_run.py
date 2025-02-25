#!/usr/bin/env python

""" Module of processes to set up the retrieval"""
from __future__ import print_function
import time
import numpy as np
import scipy as sp
import emcee
import os
import sys
import pickle
import cloud_dic
from builtins import str
from builtins import range
import utils
from schwimmbad import MPIPool
# import pymultinest as mn
import mpi4py
import test_module
import settings

__author__ = "Fei Wang"
__copyright__ = "Copyright 2024 - Fei Wang"
__credits__ = ["Fei Wang", "Ben Burningham"]
__license__ = "GPL"
__version__ = "0.2"  
__maintainer__ = ""
__email__ = ""
__status__ = "Development"


def tolerance_parameter_customized_distribution(x):
    args_instance=settings.runargs
    return np.log10((np.random.rand(x)* (max(args_instance.obspec[2,:]**2)*(0.1 - 0.01))) + (0.01*min(args_instance.obspec[2,10::3]**2))) 


def brewster_reterieval_run(re_params,model_config_instance,io_config_instance):

    args_instance=settings.runargs
    all_params,all_params_values =utils.get_all_parametres(re_params.dictionary)

    if not os.path.exists(f"{io_config_instance.outdir}"):
        os.makedirs(f"{io_config_instance.outdir}", exist_ok=True)

    
    # Write the arguments to a pickle if needed
    if (io_config_instance.make_arg_pickle > 0):
        pickle.dump(settings.runargs,open(io_config_instance.outdir+io_config_instance.runname+"_runargs.pic","wb"))

        with open(io_config_instance.outdir+io_config_instance.runname+'_configs.pkl', 'wb') as file:
            pickle.dump({
                're_params': re_params,
                'model_config': model_config_instance,
                'io_config': io_config_instance,
            }, file)

        if(io_config_instance.make_arg_pickle == 1):
            sys.exit()

        

    if re_params.samplemode=='mcmc':

        model_config_instance.ndim=len(all_params)
        model_config_instance.nwalkers=len(all_params)*16

        r2d2 = (71492e3)**2. / (model_config_instance.dist * 3.086e+16)**2.
        re_params.dictionary['refinement_params']['params']['r2d2']['distribution']=['normal', r2d2, 0.1*r2d2]

        for i in range(len(all_params)):
            if all_params[i].startswith('tolerance_parameter'):
                re_params.dictionary['refinement_params']['params'][all_params[i]]['distribution']=['customized',tolerance_parameter_customized_distribution]

        if model_config_instance.fresh == 0:
            p0=utils.MC_P0_gen(re_params.dictionary,model_config_instance,args_instance)

        elif (model_config_instance.fresh!= 0):
            fname=io_config_instance.chaindump
            pic=pickle.load(open(fname,'rb'))
            p0=pic
            if (model_config_instance.fresh == 2):
                gaslist = list(re_params.dictionary['gas'].keys())
                gastype_values = [info['gastype'] for key, info in re_params.dictionary['gas'].items() if 'gastype' in info]
                gas=[]
                for i in range(len(gaslist)):
                    gas.append(gaslist[i])
                    if  gastype_values[i]=='N':
                        gas.append("p_ref_%s"%gaslist[i])
                        gas.append("alpha_%s"%gaslist[i])
                for i in range(len(gas)):
                    p0[:,i] = (np.random.rand(model_config_instance.nwalkers).reshape(model_config_instance.nwalkers)*0.5) + p0[:,i]
                    
                    
        # Now we set up the MPI bits
        pool = MPIPool()
        if not pool.is_master():
            pool.wait()
            sys.exit()

        sampler = emcee.EnsembleSampler(model_config_instance.nwalkers, model_config_instance.ndim,test_module.lnprob,args=(re_params,),pool=pool)
        # '''
        # run the sampler
        print("running the sampler")
        clock = np.empty(80000)
        k = 0
        times = open(io_config_instance.rfile, "w")
        times.close()
        if io_config_instance.runtest == 0 and model_config_instance.fresh == 0:
            pos, prob, state = sampler.run_mcmc(p0,model_config_instance.nburn)
            sampler.reset()
            p0 = pos
        for result in sampler.sample(p0, iterations=model_config_instance.niter):
            clock[k] = time.perf_counter()
            if k > 1:
                tcycle = clock[k] - clock[k-1]
                times = open(io_config_instance.rfile, "a")
                times.write("*****TIME FOR CYCLE*****\n")
                times.write(str(tcycle))
                times.close()
            k = k+1
            position = result.coords
            f = open(io_config_instance.statfile, "w")
            f.write("****Iteration*****")
            f.write(str(k))
            f.write("****Reduced Chi2*****")
            f.write(str(result.log_prob * -2.0/(args_instance.obspec.shape[1] / 3.0)))
            f.write("****Accept Fraction*****")
            f.write(str(sampler.acceptance_fraction))
            f.write("*****Values****")
            f.write(str(result.coords))
            f.close()

            if (k==10 or k==1000 or k==1500 or k==2000 or k==2500 or k==3000 or k==3500 or k==4000 or k==4500 or k==5000 or k==6000 or k==7000 or k==8000 or k==9000 or k==10000 or k==11000 or k==12000 or k==15000 or k==18000 or k==19000 or k==20000 or k==21000 or k==22000 or k==23000 or k==24000 or k==25000 or k==26000 or k==27000 or k==28000 or k==29000 or k == 30000 or k == 35000 or k == 40000 or k == 45000 or k == 50000 or k == 55000 or k == 60000 or k == 65000):
                chain=sampler.chain
                lnprob=sampler.lnprobability
                output=[chain,lnprob]
                pickle.dump(output,open(io_config_instance.outdir+io_config_instance.picdump,"wb"))
                pickle.dump(chain[:,k-1,:], open(io_config_instance.chaindump,'wb'))



        # get rid of problematic bit of sampler object
        del sampler.__dict__['pool']


        def save_object(obj, filename):
            with open(filename, "wb") as output:
                pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

        pool.close()
        save_object(sampler, io_config_instance.outdir+io_config_instance.finalout)
    
        
    if re_params.samplemode=='multinest':

        model_config_instance.ndim=len(all_params)
        model_config_instance.nparam=len(all_params)
        
        
        def initialize_mnest():
            mnest_args = {'LogLikelihood':'',
                            'Prior':'',
                            'n_dims':0,
                            'n_params':0,
                            'n_clustering_params':None,
                            'wrapped_params':None,
                            'importance_nested_sampling':True,
                            'multimodal':True,
                            'const_efficiency_mode':False,
                            'n_live_points':400,
                            'evidence_tolerance':0.5,
                            'sampling_efficiency':0.8,
                            'n_iter_before_update':100,
                            'null_log_evidence':-1.e90,
                            'max_modes':100,
                            'mode_tolerance':-1.e90,
                            'outputfiles_basename':'',
                            'seed':-1,
                            'verbose':True,
                            'resume':True,
                            'context':0,
                            'write_output':True,
                            'log_zero':-1.e100,
                            'max_iter':0,
                            'init_MPI':False,
                            'dump_callback':None}
            return mnest_args


        # Create the LogLikelihood function with additional parameters
        def log_likelihood_call(re_params):
            def log_likelihood(cube, ndim, nparams):
                theta=cube[:ndim]
                lnLik=test_module.lnlike(theta,re_params)
                return lnLik
            return log_likelihood
        
        # # Create the Prior function with additional parameters
        # def prior_call(re_params):
        #     def prior(cube, ndim, nparams):
        #         theta=cube[:ndim]
        #         phi=test_module.priormap_dic(theta,re_params)
        #         return phi
        #     return prior

        def prior_call(re_params):
            def prior(cube, ndim, nparams):
                theta = cube[:ndim]
                phi = test_module.priormap_dic(theta, re_params)
                for i in range(ndim):
                    cube[i] = phi[i]  # Update cube with transformed values
            return prior



        mnest_args = initialize_mnest()
        mnest_args['n_params'] = model_config_instance.nparam
        mnest_args['n_dims'] = model_config_instance.ndim
        mnest_args['n_live_points'] = model_config_instance.n_live_points
        mnest_args['outputfiles_basename'] = io_config_instance.outdir+io_config_instance.runname
        
        mnest_args['LogLikelihood']= log_likelihood_call(re_params)
        mnest_args['Prior'] = prior_call(re_params)

        mnest_args['multimodal'] = model_config_instance.multimodal
        mnest_args['log_zero'] = model_config_instance.log_zero
        mnest_args['importance_nested_sampling']= model_config_instance.importance_nested_sampling
        mnest_args['sampling_efficiency']= model_config_instance.sampling_efficiency
        mnest_args['const_efficiency_mode']= model_config_instance.const_efficiency_mode
        mnest_args['evidence_tolerance']= model_config_instance.evidence_tolerance
        
        result = mn.run(**mnest_args)
        # result = mn.solve(**mnest_args)


        print()
        print('evidence: %(logZ).1f +- %(logZerr).1f' % result)
        print()
        print('parameter values:')
        for col in zip(result['samples'].transpose()):
            print('%.3f +- %.3f' % (col.mean(), col.std()))


