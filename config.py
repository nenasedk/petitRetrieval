from petitRADTRANS import nat_cst as nc
from scipy.interpolate import interp1d
from util import *

# Data Setup
RETRIEVAL_NAME = 'WISE0855_fringe'
INPUT_DIR = '/home/ipa/quanz/user_accounts/evertn/SimulatedObs/WISE0855/jwst_spectra/' # end with forward slash!
OUTPUT_DIR = '/home/ipa/quanz/user_accounts/evertn/AtmoRets/'
#OUTPUT_DIR = '/home/evert/Documents/ThesisNotebooks/AtmoRets/' 

# Observation file format:
# wlen (micron), F_pl/F_s, F_pl/F_s (error)
# Space separated
observation_files = {}
#observation_files['NIRISS SOSS'] = 'NIRISS_SOSS_flux.dat'
#observation_files['NIRSpec G395M'] = 'NIRSpec_G395M_flux.dat'
#observation_files['MIRI LRS'] = 'MIRI_LRS_flux.dat'
observation_files['MIRIMRS'] = 'combined_cal_spectra.dat'

# Pymultinest and Verbosity Parameters
PLOTTING = False # (not used, seems to break multinest)
LIVE = 1600
WRITE_THRESHOLD = 200

# Wavelength range of observations, fixed parameters that will not be retrieved
# Any of these parameters can be floating, but must be included in the PRIORS dictionary
OBJ_NAME = 'WISE0855'
IS_COMPANION = False # False if field object
WLEN = [4.5, 18.5]
LOG_G = 3.20
R_pl =   1.17*nc.r_jup_mean
R_star = 1.81*nc.r_sun
D_pl = 2.23 * nc.pc

# Get host star spectrum to calculate F_pl / F_star later.
T_star = 6295.
x = nc.get_PHOENIX_spec(T_star)
fstar = interp1d(x[:,0], x[:,1])

# Parameters to be retrieved
# I really hate this. Should use BIC to determine how many parameters to use.
RAYLEIGH_SPECIES = ['H2','He']
CONTINUUM_OPACITIES = ['H2-H2','H2-He']
CLOUD_SPECIES = []
#CLOUD_SPECIES = ['Mg2SiO4(c)_cd']
LINE_SPECIES =  ['CO_all_iso','H2O','CH4','NH3','CO2','HCN','H2S','C2H2']

# MUST INCLUDE: 'log_delta','log_gamma','t_int','t_equ','log_p_trans','alpha','log_g','log_P0'
ATMOSPHERE = ['log_delta','log_gamma','t_int','t_equ','log_p_trans','alpha','log_g','log_P0']#,'Kzz','fsed','sigma_lnorm']
parameters = []
parameters.extend(ATMOSPHERE)
parameters.extend(CLOUD_SPECIES)
parameters.extend(LINE_SPECIES)


# Prior functions on parameter values
PRIORS = {}
# use uniform prior rather than log prior because everything loglike function uses exp(y)
PRIORS['log_delta']      = lambda x: gaussian_prior(x,-5.5,2.5)                          
PRIORS['log_gamma']      = lambda x: gaussian_prior(x,0.0,2.0)
PRIORS['t_int']          = lambda x: uniform_prior(x, 0., 3500.)
PRIORS['t_equ']          = lambda x: delta_prior(x, 3.4, 3.4)
PRIORS['log_p_trans']    = lambda x: gaussian_prior(x,-3.,3.)
PRIORS['alpha']          = lambda x: gaussian_prior(x,0.25,0.4)
PRIORS['log_g']          = lambda x: uniform_prior(x, 2.0, 4.5) 
PRIORS['log_P0']         = lambda x: uniform_prior(x, -6, 2.)

# Boundary sanity checks on parameter values
# Speeds up computation when non-physical values are tried.
LOG_PRIORS = {}
LOG_PRIORS['log_delta']      = lambda x: -((x-(-5.5))/2.5)**2./2.                           
LOG_PRIORS['log_gamma']      = lambda x: -((x-(-0.0))/2.)**2./2. 
LOG_PRIORS['t_int']          = lambda x: a_b_range(x, 0., 3500.) 
LOG_PRIORS['t_equ']          = lambda x: a_b_range(x, 3., 4.)
LOG_PRIORS['log_p_trans']    = lambda x: -((x-(-3))/3.)**2./2.
LOG_PRIORS['alpha']          = lambda x: -((x-0.25)/0.4)**2./2.
LOG_PRIORS['log_g']          = lambda x: a_b_range(x, 2.0, 4.5) 
LOG_PRIORS['log_P0']         = lambda x: a_b_range(x, -6, 2.)

for line in LINE_SPECIES:
    LOG_PRIORS[line] = lambda x: a_b_range(x, -12., 0.)
    PRIORS[line] = lambda x: uniform_prior(x, -12., 0.)
if len(CLOUD_SPECIES) > 0:
    for cloud in CLOUD_SPECIES:
        LOG_PRIORS[cloud] = lambda x: a_b_range(x, -12., 0.)
        PRIORS[cloud] = lambda x: uniform_prior(x, -12., 0.)
