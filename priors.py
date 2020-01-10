import numpy as np
import sys
import os

from petitRADTRANS import nat_cst as nc
class Prior:
    def __init__(self,
                 log_prior_dict = None,
                 prior_dict = None):
        """
        Prior class

        This class stores the dictionaries of prior functions and 
        sanity checks. 

        Parameters
        ----------------
        log_prior_dict: dict
            This dictionary stores sanity checks on parameter ranges
            for uniform priors, and computes the exponent of a gaussian
            distribution for gaussian distributed parameters
        prior_dict: dict
            This dictionary stores the actual prior functions which transform
            samples from a unit hypercube into physical parameter space.
        """
        self.log_priors = {}
        self.priors = {}
        self.n_params = 0
        if log_prior_dict is not None:
            self.log_priors = log_prior_dict
        if prior_dict is not None:
            self.priors = prior_dict
            self.n_params = len(self.priors.keys())

    def getLogPriors(self):
        return self.log_priors
    
    def getPriors(self):
        return np.array(self.priors)
    
    def addLogPrior(self,param, value):
        self.log_priors[param] = value
        self.n_params += 1
        
    def setLogDict(self,logdict):
        self.log_priors = logdict
        
    def addPrior(self,param,value):
        self.priors[param] = value
        
    def setPriors(self,prior_dict):
        self.priors = prior_dict
        self.n_params = len(self.priors.keys())


        
