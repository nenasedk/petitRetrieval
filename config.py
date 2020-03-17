from petitRADTRANS import nat_cst as nc
from scipy.interpolate import interp1d
from petitRetrieval.util import *

# Data Setup
RETRIEVAL_NAME = 'VHS_nonoise_nocorr_CH1'
INPUT_DIR = '/home/ipa/quanz/user_accounts/evertn/SimulatedObs/VHS1256b_lowered_nonoise_rfc/jwst_spectra_nocorr/' # end with forward slash!
OUTPUT_DIR = '/home/ipa/quanz/user_accounts/evertn/AtmoRets/'

# Observation file format:
# wlen (micron), F_pl/F_s, F_pl/F_s (error)
# Space separated
observation_files = {}
observation_files['MIRIMRS'] = 'combined_cal_spectra.dat'

# Pymultinest and Verbosity Parameters
PLOTTING = False # (not used, seems to break multinest)
LIVE = 500
WRITE_THRESHOLD = 200

# Wavelength range of observations, fixed parameters that will not be retrieved
# Any of these parameters can be floating, but must be included in the PRIORS dictionary
OBJ_NAME = 'VHS1256b' #R_pl =   1.17*nc.r_jup_mean
IS_COMPANION = False # False if field object
WLEN = [4.5, 18.5]
LOG_G = 3.20 # Obsolete, included as parameter

# Fixed Parameters
FIXED_PARAMS = {}
FIXED_PARAMS['t_equ'] = 2.7
FIXED_PARAMS['R_star'] = 1.81*nc.r_sun
FIXED_PARAMS['D_pl'] = 12.7 * nc.pc

# Get host star spectrum to calculate F_pl / F_star later.
# Only used if using a contrast measurement.
T_star = 6295.
x = nc.get_PHOENIX_spec(T_star)
fstar = interp1d(x[:,0], x[:,1])


# =======================================
# TEMPERATURE-PRESSURE PROFILE PARAMETERS
# =======================================
# One set of P-T profile parameters must be included in
# the prior dictionary.
#
# Modified Guillot Parameters:
# 'log_delta', 'log_gamma', 't_int', 't_equ', 'log_p_trans', 'alpha'
#
# Guillot Global Parameters
# 'log_delta', 'log_gamma', 't_int', 't_equ' OR
# 'log_g', 'kappa_IR', 'log_gamma', 't_int', 't_equ',
# ========================================


# ========================================
# PRIOR FUNCTIONS
# ========================================
# Parameters to be retrieved
# Should use BIC to determine how many parameters to use automatically.
RAYLEIGH_SPECIES = ['H2','He']
CONTINUUM_OPACITIES = ['H2-H2','H2-He']
CLOUD_SPECIES = [] #['Mg2SiO4(c)_cd']
LINE_SPECIES =  ['CO_all_iso','H2O','CH4','NH3','CO2','H2S','C2H2','TiO']

PRIORS = {}
# use uniform prior rather than log prior because everything loglike function uses exp(y)
#PRIORS['log_delta']      = lambda x: gaussian_prior(x,-5.5,2.5)                          
PRIORS['log_gamma']      = lambda x: gaussian_prior(x,0.0,2.0)
PRIORS['t_int']          = lambda x: uniform_prior(x, 0., 3500.)
PRIORS['R_pl']           = lambda x: uniform_prior(x, 0.01, 2.5)
PRIORS['kappa_IR']       = lambda x: uniform_prior(x,0.,1.)
#PRIORS['log_p_trans']    = lambda x: gaussian_prior(x,-3.,3.)
#PRIORS['alpha']          = lambda x: gaussian_prior(x,0.25,0.4)
PRIORS['log_g']          = lambda x: uniform_prior(x, 2.0, 4.5) 
PRIORS['log_P0']         = lambda x: uniform_prior(x, -6, 2.)

ATMOSPHERE = PRIORS.keys()
parameters = []
parameters.extend(ATMOSPHERE)
parameters.extend(CLOUD_SPECIES)
parameters.extend(LINE_SPECIES)

# Boundary sanity checks on parameter values
# Speeds up computation when non-physical values are tried.
RANGE = {}
#RANGE['log_delta']      = lambda x: -((x-(-5.5))/2.5)**2./2.                           
RANGE['log_gamma']      = lambda x: -((x-(-0.0))/2.)**2./2. 
RANGE['t_int']          = lambda x: a_b_range(x, 0., 3500.) 
RANGE['R_pl']           = lambda x: a_b_range(x, 0.01, 2.5)
#RANGE['log_p_trans']    = lambda x: -((x-(-3))/3.)**2./2.
#RANGE['alpha']          = lambda x: -((x-0.25)/0.4)**2./2.
RANGE['kappa_IR']       = lambda x: a_b_range(x,0.,1.)
RANGE['log_g']          = lambda x: a_b_range(x, 2.0, 4.5) 
RANGE['log_P0']         = lambda x: a_b_range(x, -6, 2.)

for line in LINE_SPECIES:
    RANGE[line] = lambda x: a_b_range(x, -12., 0.)
    PRIORS[line] = lambda x: uniform_prior(x, -12., 0.)
if len(CLOUD_SPECIES) > 0:
    for cloud in CLOUD_SPECIES:
        RANGE[cloud] = lambda x: a_b_range(x, -12., 0.)
        PRIORS[cloud] = lambda x: uniform_prior(x, -12., 0.)
