"""
Retrieval Class

author: Evert Nasedkin
contact: nasedkinevert@gmail.com
"""
import numpy as np
from matplotlib import pyplot as plt
import sys, os
import json
from scipy.interpolate import interp1d

from config import *
from .util import show, Surf_To_Meas
from .master_retrieval_model import retrieval_model_plain, calc_MMW
from .rebin_give_width import rebin_give_width

from petitRADTRANS import Radtrans
from petitRADTRANS import nat_cst as nc

# Declare diagnostics
function_calls = 0
computed_spectra = 0
NaN_spectra = 0
write_threshold =  delta_wt = WRITE_THRESHOLD

class Retrieval:
    def __init__(self,
                 data_obj,
                 rt_obj,
                 prior_obj,
                 name_in = "",
                 data_path = "",
                 output_path= "",
                 live = 1000,
                 plotting = False,
                 diagnostics = False,
                 resume = False,
                 verbose = True,
                 sampling_efficiency = 0.3):
        """
        Retrieval Class
        
        This class implements the Prior and Loglikelihood functions
        as required by pyMultinest.
        
        The Prior function takes in samples from a unit hypercube of ndim 
        dimensions, and transforms them to physical parameter space. The actual 
        prior functions used for this are located in util.
        
        The Loglikelihood function takes a vector of parameters sampled from the 
        unit hypercube and transformed by the prior function, and uses these 
        parameters to compute a model spectrum using petitRadTrans.This model is 
        then compared to data to compute the loglikelihood.
        
        Parameters
        ------------------------
        data_obj: Data
            An instance of the Data class implemented in data.py. This class stores
            the measured flux in each wavelength bin.
        rt_obj: Radtrans
            An instance of the Radtrans class implemented by petitRadTrans. This object
            is used to compute a model atmospheric spectrum given parameters
        prior_obj: Prior
            An instance of the Prior class implemented in priors.py. It contains the
            dictionary of prior functions used to transform the unit hypercube sample
            space.
        name_in: str
            Name for outputting files.
        data_path: str
            Location of data files. Not needed.
        output_path: str
            Location of output files.
        live: int
            Number of live points used for sampling in Multinest
        plotting: bool
            Plot figures while running. Do not use if running full retrieval
        diagnostics: bool
            Print information about each sample. Do not use for full retrieval
        resume: bool
            Pymultinest parameter
        verbose: bool
            Pymultinest parameter
        sampling_efficiency: float [0,1]
            Pymultinest parameter
        Returns
        ------------------
        Prior and loglikelihood functions to be used by pymultinest
        """
        self.data = data_obj
        self.rt_obj = rt_obj
        self.prior_obj = prior_obj
        self.data_directory = data_path
        self.output_directory = output_path
        self.name_in = name_in
        
        self.live = live
        self.plotting = plotting
        self.diagnostics = diagnostics
        self.ndim = self.prior_obj.n_params
        self.analyzer = None

        self.resume = resume
        self.verbose = verbose
        self.efficiency = sampling_efficiency
        return
    
    def Prior(self,cube,ndim,nparam):
        """
        Prior

        This function transforms samples from the unit hypercube to samples
        from an arbitrary distribution specified by the prior_obj.
        
        For example, to generate a uniform distribution between hi and lo
        cube[0] = cube[0] * (hi - lo) + lo
        
        Most priors are based on 
        https://github.com/JohannesBuchner/MultiNest/blob/master/src/priors.f90

        Parameters
        --------------------
        cube: ndarray?
             Actually a ctype pointer or something. Use like a np array. A vector of ndim samples 
             from the unit hypercube, with each entry in [0,1]
        ndim:
             Number of dimensions, typically equal to nparam
        nparam:
             Total number of parameters, can be different from ndim if not all are being sampled
        """
        ind = 0
        for key in self.prior_obj.log_priors.keys():
            # This is a dictionary of lambda functions that will output the correct transformed sample
            cube[ind] = self.prior_obj.priors[key](cube[ind])
            ind += 1
            
    def LogLikelihood(self,cube,ndim,nparam):
        return self.computeLogLikelihood(cube,ndim)
    
    def computeLogLikelihood(self,cube,ndim):
        """
        computeLogLikelihood

        This function implements the log-likelihood function for an emission spectra
        using a petitRadTrans model. It is more fully described in:
        https://arxiv.org/abs/1904.11504,
        https://petitradtrans.readthedocs.io/en/latest/content/notebooks/ret_emission_master.html

        Parameters
        -------------------
        cube: ndarray?
            The transformed samples from the unit hypercube, transformed by Prior
        ndim: int
            Number of dimensions
        """
        params = []
        for i in range(ndim):
            params.append(cube[i])
        params = np.array(params)
        nlines = len(LINE_SPECIES)
        nclouds = len(CLOUD_SPECIES)

        # Make dictionary for modified Guillot parameters
        temp_params = {}
        for i,param in enumerate(ATMOSPHERE):
            temp_params[param] = params[i] 
        for key,value in FIXED_PARAMS.items():
            temp_params[key] = value
            
        # Make dictionary for log 'metal' abundances
        ab_metals = {}
        for i,line in enumerate(LINE_SPECIES):
            ab_metals[line] = params[-nlines:][i]

        if self.diagnostics:
            global function_calls
            global computed_spectra
            global NaN_spectra
            global write_threshold   
            function_calls += 1
        
        # Prior calculation of all input parameters
        log_prior = 0.
        
        # Alpha should not be smaller than -1, this
        # would lead to negative temperatures!
        if temp_params['alpha'] < -1:
            if self.diagnostics:
                print("Alpha too small!")
            return -np.inf
        
        for key in temp_params.keys():
            # Fixed parameters aren't in the prior dictionary
            if key in self.prior_obj.log_priors:
                log_prior += self.prior_obj.log_priors[key](temp_params[key])         
                log_prior += self.prior_obj.log_priors['log_g'](temp_params['log_g'])
                log_prior += self.prior_obj.log_priors['log_P0'](temp_params['log_P0'])
                if log_prior == -np.inf:
                    if self.diagnostics:
                        print(self.prior_obj.log_priors[key](temp_params[key]),\
                              self.prior_obj.log_priors['log_g'](temp_params['log_g']),\
                              self.prior_obj.log_priors['log_P0'](temp_params['log_P0']))
                    return -np.inf
        # Metal abundances: check that their
        # summed mass fraction is below 1.
        metal_sum = 0.
        for name in ab_metals.keys():
            log_prior += self.prior_obj.log_priors[name](ab_metals[name])
            metal_sum += 1e1**ab_metals[name]
                
        if metal_sum > 1.:
            if self.diagnostics:
                print("Metal sum > 1")
            log_prior += -np.inf
                    
        # Return -inf if parameters fall outside prior distribution
        if (log_prior == -np.inf):
            return -np.inf
    
        # Calculate the log-likelihood
        log_likelihood = 0.
        
        # Calculate the forward model, this
        # retur<ns the wavelengths in cm and the flux F_nu
        # in erg/cm^2/s/Hz
        wlen, flux_nu = retrieval_model_plain(self.rt_obj, temp_params, ab_metals)

        # Just to make sure that a long chain does not die
        # unexpectedly:
        # Return -inf if forward model returns NaN values
        if np.sum(np.isnan(flux_nu)) > 0:
            print("NaN spectrum encountered")
            if self.diagnostics:
                NaN_spectra += 1
            return -np.inf
            self.n_params = len(self.log_priors.keys())
        # Convert to observation for emission case
        if IS_COMPANION:
            flux_star = fstar(wlen)
            flux_sq   = flux_nu/flux_star*(temp_params['R_pl']*nc.r_jup_mean/temp_params['R_star'])**2 
        else:
            # TODO - Check to see if this is the correct way to retrieve a field BD
            # should correct for distance (reduce model by d**-2 or normalize data)
            # distance should be fixed or floating param?
            # - probably depends on what data is available (GAIA?)
            flux_sq = Surf_To_Meas(flux_nu,temp_params['R_pl']*nc.r_jup_mean,temp_params['D_pl'])
        # Calculate log-likelihood
        for instrument in self.data.data_wlen.keys():
            # Rebin model to observation
            flux_rebinned = rebin_give_width(wlen, flux_sq, \
                                             self.data.data_wlen[instrument],\
                                             self.data.data_wlen_bins[instrument])
            if self.plotting:
                plt.errorbar(self.data.data_wlen[instrument], \
                             self.data.data_flux_nu[instrument], \
                             self.data.data_flux_nu_error[instrument], \
                             fmt = 'o', \
                             zorder = -20, \
                             color = 'red')
                
                plt.plot(self.data.data_wlen[instrument], \
                         flux_rebinned, \
                         's', \
                         zorder = -20, \
                         color = 'blue')
                
            # Calculate log-likelihood
            log_likelihood += -np.sum(((flux_rebinned - self.data.data_flux_nu[instrument])/ \
                                       self.data.data_flux_nu_error[instrument])**2.)/2.

        if self.plotting:
            plt.plot(wlen, flux_sq, color = 'black')
            plt.xscale('log')
            plt.show()
            
        if self.diagnostics:
            computed_spectra += 1
            if (function_calls >= write_threshold):
                
                write_threshold += delta_wt
                #hours = (time.time() - start_time)/3600.0
                info_list = [function_calls, computed_spectra, NaN_spectra, \
                             log_prior + log_likelihood, hours] 
                
                file_object = open(self.output_directory + 'diag_' + self.name_in + '.dat', 'a')

                for i in np.arange(len(info_list)):
                    if (i == len(info_list) - 1):
                        file_object.write(str(info_list[i]).ljust(15) + "\n")
                    else:
                        file_object.write(str(info_list[i]).ljust(15) + " ")
                file_object.close()
            print(log_prior + log_likelihood)
            print("--> ", function_calls, " --> ", computed_spectra)
            
        if np.isnan(log_prior + log_likelihood):
            return -np.inf
        else:
            return log_likelihood # + log_prior
