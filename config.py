from petitRADTRANS import nat_cst as nc
from scipy.interpolate import interp1d
from util import *

# Data Setup
RETRIEVAL_NAME = 'JWST_emission_g395M_v5'
INPUT_DIR = 'data/' # end with forward slash!
OUTPUT_DIR = '/home/ipa/quanz/user_accounts/evertn/AtmoRets/'
#OUTPUT_DIR = '/home/evert/Documents/ThesisNotebooks/AtmoRets/' 

# Observation file format:
# wlen (micron), F_pl/F_s, F_pl/F_s (error)
# Space separated
observation_files = {}
observation_files['NIRISS SOSS'] = 'NIRISS_SOSS_flux.dat'
observation_files['NIRSpec G395M'] = 'NIRSpec_G395M_flux.dat'
observation_files['MIRI LRS'] = 'MIRI_LRS_flux.dat'


# Pymultinest and Verbosity Parameters
PLOTTING = False # (not used, seems to break multinest)
LIVE = 500
WRITE_THRESHOLD = 200

# Wavelength range of observations, fixed parameters that will not be retrieved
OBJ_NAME = 'g395M'
IS_COMPANION = True # False if brown dwarf
WLEN = [0.8, 14.0]
LOG_G =  2.58
R_pl =   1.84*nc.r_jup_mean
R_star = 1.81*nc.r_sun

# Get host star spectrum to calculate F_pl / F_star later.
T_star = 6295.
x = nc.get_PHOENIX_spec(T_star)
fstar = interp1d(x[:,0], x[:,1])

# Parameters to be retrieved
# I really hate this. Should use BIC to determine how many parameters to use.
parameters = ['log_delta','log_gamma','t_int','t_equ','log_p_trans','alpha','log_g','log_P0','CO_all_iso','H2O','CH4','NH3','CO2','H2S','Na','K']

# Boundary sanity checks on parameter values
LOG_PRIORS = {}
LOG_PRIORS['log_delta']      = lambda x: -((x-(-5.5))/2.5)**2./2.                           
LOG_PRIORS['log_gamma']      = lambda x: -((x-(-0.0))/2.)**2./2. 
LOG_PRIORS['t_int']          = lambda x: a_b_range(x, 0., 1500.)
LOG_PRIORS['t_equ']          = lambda x: a_b_range(x, 0., 4000.)
LOG_PRIORS['log_p_trans']    = lambda x: -((x-(-3))/3.)**2./2.
LOG_PRIORS['alpha']          = lambda x: -((x-0.25)/0.4)**2./2.
LOG_PRIORS['log_g']          = lambda x: a_b_range(x, 2.0, 3.7) 
LOG_PRIORS['log_P0']         = lambda x: a_b_range(x, -4, 2.)
LOG_PRIORS['CO_all_iso']     = lambda x: a_b_range(x, -10., 0.)
LOG_PRIORS['H2O']            = lambda x: a_b_range(x, -10., 0.)
LOG_PRIORS['CH4']            = lambda x: a_b_range(x, -10., 0.)
LOG_PRIORS['NH3']            = lambda x: a_b_range(x, -10., 0.)
LOG_PRIORS['CO2']            = lambda x: a_b_range(x, -10., 0.)
LOG_PRIORS['H2S']            = lambda x: a_b_range(x, -10., 0.)
LOG_PRIORS['Na']             = lambda x: a_b_range(x, -10., 0.)
LOG_PRIORS['K']              = lambda x: a_b_range(x, -10., 0.)

PRIORS = {}
# Prior functions on parameter values
# use uniform prior rather than log prior because everything loglike function uses 10**y
PRIORS['log_delta']      = lambda x: gaussian_prior(x,-5.5,2.5)                          
PRIORS['log_gamma']      = lambda x: gaussian_prior(x,0.0,2.0)
PRIORS['t_int']          = lambda x: uniform_prior(x, 0., 1500.)
PRIORS['t_equ']          = lambda x: uniform_prior(x, 0., 4000.)
PRIORS['log_p_trans']    = lambda x: gaussian_prior(x,-3.,3.)
PRIORS['alpha']          = lambda x: gaussian_prior(x,0.25,0.4)
PRIORS['log_g']          = lambda x: uniform_prior(x, 2.0, 3.7) 
PRIORS['log_P0']         = lambda x: uniform_prior(x, -4, 2.)

# Priors for log mass fractions
PRIORS['CO_all_iso']     = lambda x: uniform_prior(x, -10., 0.)
PRIORS['H2O']            = lambda x: uniform_prior(x, -10., 0.)
PRIORS['CH4']            = lambda x: uniform_prior(x, -10., 0.)
PRIORS['NH3']            = lambda x: uniform_prior(x, -10., 0.)
PRIORS['CO2']            = lambda x: uniform_prior(x, -10., 0.)
PRIORS['H2S']            = lambda x: uniform_prior(x, -10., 0.)
PRIORS['Na']             = lambda x: uniform_prior(x, -10., 0.)
PRIORS['K']              = lambda x: uniform_prior(x, -10., 0.)

