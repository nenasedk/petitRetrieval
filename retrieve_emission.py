"""
petitRetrieval

Retreive Emission

This script uses pymultinest to sample the prior parameter space in order to compute
posterior distributions and Bayesian evidence based on the the log-likelihood
function implemented in retrieval.Retrieval

Author: Evert Nasedkin

Changelog:
Jan 10 2020 - Initial construction

TODO:
- Compare having the pymultinest.run as part of the script vs
  part of the retrieval class.

- Find a better way of dealing with global variables
- Automated generation of config files + config parser
- Notification option
- Improve corner plotting

Credits:
petitRadTrans - https://arxiv.org/abs/1904.11504
HELIOS-t - https://arxiv.org/abs/1610.03216
Multinest - https://arxiv.org/abs/0809.3437
pyMultinest - https://arxiv.org/abs/1402.0004
"""
from __future__ import absolute_import, unicode_literals, print_function

from config import *
from util import Surf_To_Meas
from priors import Prior
from data import Data
from retrieval import Retrieval
from corner_plots import CornerPlot

import pymultinest
import numpy as np
from astropy.io import ascii,fits
from astropy.table import Table,Column
import sys, os
import time
import json
from matplotlib import pyplot as plt

from petitRADTRANS import Radtrans
import master_retrieval_model as rm
from petitRADTRANS import nat_cst as nc
import rebin_give_width as rgw
from scipy.interpolate import interp1d

sys.stdout.flush()
start = time.time()

### Create and setup radiative transfer object
# Create random P-T profile to create RT arrays of the Radtrans object.
temp_params = {}
temp_params['log_delta'] = -6.
temp_params['log_gamma'] = np.log10(0.4)
temp_params['t_int'] = 750.
temp_params['t_equ'] = 0.
temp_params['log_p_trans'] = -3.
temp_params['alpha'] = 0.
p, t = nc.make_press_temp(temp_params)

# Create the Radtrans object here
rt_object = Radtrans(line_species=LINE_SPECIES, \
                    rayleigh_species=RAYLEIGH_SPECIES, \
                    continuum_opacities = CONTINUUM_OPACITIES, \
                    mode='c-k', \
                    wlen_bords_micron=WLEN)
rt_object.setup_opa_structure(p)

# Read in data and convert to CGS
data = Data(observation_files,
            RETRIEVAL_NAME,
            INPUT_DIR,
            OUTPUT_DIR + RETRIEVAL_NAME + '/')

# Ensure output directories exist, and change to output dir to shorten path name
# Multinest has a 100 char limit for paths (recently updated to 1000?)
if not os.path.exists(OUTPUT_DIR + RETRIEVAL_NAME + "/"):
    os.mkdir(OUTPUT_DIR + RETRIEVAL_NAME + "/")
if not os.path.isfile(OUTPUT_DIR + RETRIEVAL_NAME + "/" + 'ev.dat'):
    f = open(OUTPUT_DIR  + RETRIEVAL_NAME + "/"+ 'ev.dat','w+')
    f.close()
os.chdir(OUTPUT_DIR)

# Setup Priors
prior = Prior(LOG_PRIORS,PRIORS)

# Create object to compute priors and log likelihood
retrieval = Retrieval(data_obj = data,
                      rt_obj = rt_object,
                      prior_obj = prior,
                      planet_radius = R_pl,
                      star_radius = R_star,
                      planet_distance = D_pl,
                      name_in = RETRIEVAL_NAME,
                      data_path = INPUT_DIR,
                      output_path = RETRIEVAL_NAME + '/',
                      plotting = PLOTTING,
                      live = LIVE,
                      diagnostics = False)

# Run retrieval using pymultinest
# outputfiles_basename must be <100 chars
# https://johannesbuchner.github.io/PyMultiNest/pymultinest_run.html
# init_MPI should be False even when using MPI for parallelization!

pymultinest.run(retrieval.LogLikelihood,
                retrieval.Prior,
                retrieval.ndim,
                outputfiles_basename=retrieval.output_directory + retrieval.name_in+'_',
                resume=retrieval.resume,
                verbose=retrieval.verbose,
                n_live_points=retrieval.live,
                sampling_efficiency = retrieval.efficiency, # 0.3 for parameter est, 0.8 for model comp
                importance_nested_sampling=True, # True has higher memory requirements
                write_output=True)    

# Read output files for analysis
a = pymultinest.Analyzer(outputfiles_basename=retrieval.output_directory \
                                + retrieval.name_in + '_', n_params=retrieval.ndim)
# Save in human readable format
with open('%sparams.json' % a.outputfiles_basename, 'w') as f:
    json.dump(parameters, f, indent=2)
stats = a.get_stats()
with open('%sstats.json' % a.outputfiles_basename, mode='w') as f:
    json.dump(stats, f, indent=2)
bestfit_params = a.get_best_fit()
with open('%sbestfit.json' % a.outputfiles_basename, mode='w') as f:
    json.dump(bestfit_params, f, indent=2)    
samples = a.get_equal_weighted_posterior()[:, :-1]

################################################################################
################################################################################
### Analysis (From HELIOS)
################################################################################
################################################################################

retrieved_results = list(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                             zip(*np.percentile(samples, [16, 50, 84], axis=0))))
plot_percentiles = []
new_param_dict = {}
retrieved_parameters_full = {}
retrieved_parameters_list = []
for i in range(prior.n_params):
    retrieved_parameters_full[parameters[i]] = retrieved_results[i]
    new_param_dict[parameters[i]] = bestfit_params['parameters'][i]
    retrieved_parameters_list.append(retrieved_results[i][0])
    plot_percentiles.append(retrieved_results[i])


print("""PyMultinest result:""")
for param in parameters:
    print(param,' = ',retrieved_parameters_full[param][0],' +',
          retrieved_parameters_full[param][1],' -',retrieved_parameters_full[param][2])

a_lnZ = a.get_stats()['global evidence']
log_evidence = a_lnZ/np.log(10)

print
print ('************************')
print ('MAIN RESULT: Evidence Z ')
print ('************************')
print ('  log Z for model = %.1f' % (log_evidence))
print

################################################################################
################################################################################
### Plotting and model output
################################################################################
################################################################################

retrieved_parameters_list = np.array(retrieved_parameters_list)
log_delta, log_gamma, t_int, t_equ, log_p_trans, alpha, log_g, log_P0 = retrieved_parameters_list[:-1 * len(LINE_SPECIES)]

# Make dictionary for modified Guillot parameters
temp_params = {}
temp_params['log_delta'] = log_delta
temp_params['log_gamma'] = log_gamma
temp_params['t_int'] = t_int
temp_params['t_equ'] = t_equ
temp_params['log_p_trans'] = log_p_trans
temp_params['alpha'] = alpha

# Make dictionary for log 'metal' abundances
ab_metals = {}
for i,line in enumerate(LINE_SPECIES):
    ab_metals[line] = retrieved_parameters_list[-1*len(LINE_SPECIES):][i]


## compute model for retrieved results ##
wlen, flux_nu = rm.retrieval_model_plain(rt_object, temp_params, log_g, \
                                         log_P0, R_pl, ab_metals)
flux_nu = Surf_To_Meas(fluxes,R_pl,D_pl.value)
output = Table([wlen,flux_nu],names = ['wavelength','flux_nu'])
ascii.write(output, RETRIEVAL_NAME + '/' + RETRIEVAL_NAME + "_BestFitModel.dat", overwrite=True)

# Corner plots
plots = CornerPlot(analyzer = a,
                   name_in = RETRIEVAL_NAME + "_corner_plots",
                   input_directory = RETRIEVAL_NAME + '/',
                   output_directory = RETRIEVAL_NAME + '/',
                   ndim = prior.n_params)
plots.basic_plot(RETRIEVAL_NAME + "_marginals")

end = time.time()
print("time = ",end-start)
