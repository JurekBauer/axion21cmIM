from units import *
import numpy as np

# Define fiducial cosmology and parameters (not needed if input file is specified)
cosmo_dictionary = {
    # Cosmological parameters from Hlozek et al. (2017) "Future CMB tests...".
    'omega_M_0': (0.1197 + 0.02222 + 0.06 / 93.14) / 0.69 ** 2.,
    'omega_lambda_0': 1.0 - (0.1197 + 0.02222 + 0.06 / 93.14) / 0.69 ** 2.,  # 0.72,
    'omega_cdm_0': 0.1197 / 0.69 ** 2. * (1.0 - 0.02),
    'omega_b_0': 0.02222 / 0.69 ** 2.,
    'omega_ax_0': 0.02 * 0.1197 / 0.69 ** 2.,
    'ma': -28,
    'h': 0.69,
    'ns': 0.9655,
    'As': 2.1955e-9,
    'mnu': 0.06,
    'omega_HI_0': 2.45e-4 / 0.69,  # 4.86e-4, #6.50e-4,
    'k_piv': 0.05,  # n_s
    # components for the calculation of the power spectrum, total matter power spectrum will be loaded anyways.
    'components_for_P': ('CDM', 'baryon'),  # ('total',), #
    # range and number of mesh points in between for the halo mass, in order to calculate the HI bias or omega_HI
    'M_min': 8,
    'M_max': 17,
    'N_mass_mesh': 2000,
    # which M_HI-M relation to use and additional arguments for the relation
    'HI_halo_formula': 'PRA2017',
    'HI_halo_formula_args': {'alpha': 0.09, 'beta': -0.58, 'vc0': 10 ** 1.56},
    'gamma_HI': 1.45,
    'c_HI0': 28.65,
    'sigma_8': 0.811,
}

SURVEY = {
    'ttot': 10e3 * HRS_MHZ,  # Total integration time [MHz^-1]
    'nu_line': 1420.406,  # Rest-frame freq. of emission line [MHz]
}

################################################################################
# Experiments
################################################################################

SKA1MIDfull = {
    'mode': 'dish',  # Interferometer or single dish
    'Ndish': 190,  # No. of dishes
    'Nbeam': 1,  # No. of beams (for multi-pixel detectors)
    'Ddish': 15,  # Single dish diameter [m]
    'Tinst': 28. * (1e3),  # System temp. [mK]
    'survey_dnutot': 1070.,  # Total bandwidth of *entire* survey [MHz]
    'survey_numax': 1420.,  # Max. freq. of survey
    'dnu': 0.1,  # Bandwidth of single channel [MHz]
    'Sarea': 25e3 * (D2RAD) ** 2.,  # Total survey area [radians^2]
    'n(x)': "array_config/nx_SKAM190_dec30.dat"  # Interferometer antenna density
}
SKA1MIDfull.update(SURVEY)

exptS = {
    'mode': 'dish',  # Interferometer or single dish
    'Ndish': 1,  # No. of dishes
    'Nbeam': 50,  # No. of beams (for multi-pixel detectors)
    'Ddish': 30.,  # Single dish diameter [m]
    'Tinst': 50. * (1e3),  # System temp. [mK]
    'survey_dnutot': 300.,  # Total bandwidth of *entire* survey [MHz]
    'survey_numax': 1100.,  # Max. freq. of survey
    'dnu': 0.1,  # Bandwidth of single channel [MHz]
    'Sarea': 5e3 * (D2RAD) ** 2.  # Total survey area [radians^2]
}
exptS.update(SURVEY)

exptM = {
    'mode': 'interferom',  # Interferometer or single dish
    'Ndish': 160,  # No. of dishes
    'Nbeam': 1,  # No. of beams (for multi-pixel detectors)
    'Ddish': 4.,  # Single dish diameter [m]
    'Tinst': 35. * (1e3),  # System temp. [mK]
    'nu_crit': 1000.,  # Critical frequency [MHz]
    'survey_dnutot': 400.,  # Total bandwidth of *entire* survey [MHz]
    'survey_numax': 1000.,  # Max. freq. of survey
    'dnu': 0.1,  # Bandwidth of single channel [MHz]
    'Sarea': 2e3 * (D2RAD) ** 2.,  # Total survey area [radians^2]
    'Dmax': 53.,  # Max. interferom. baseline [m]
    'Dmin': 4.  # Min. interferom. baseline [m]
}
exptM.update(SURVEY)

exptL = {
    'mode': 'combined',  # Interferometer or single dish
    'Ndish': 250,  # No. of dishes
    'Nbeam': 1,  # No. of beams (for multi-pixel detectors)
    'Ddish': 15.,  # Single dish diameter [m]
    'Tinst': 20. * (1e3),  # System temp. [mK]
    'survey_dnutot': 700.,  # Total bandwidth of *entire* survey [MHz]
    'survey_numax': 1100.,  # Max. freq. of survey
    'dnu': 0.1,  # Bandwidth of single channel [MHz]
    'Sarea': 25e3 * (D2RAD) ** 2.,  # Total survey area [radians^2]
    'Dmax': 600.,  # Max. interferom. baseline [m]
    'Dmin': 15.  # Min. interferom. baseline [m]
}
exptL.update(SURVEY)

GBT = {
    'mode': 'dish',  # Interferometer or single dish
    'Ndish': 1,  # No. of dishes
    'Nbeam': 1,  # No. of beams (for multi-pixel detectors)
    'Ddish': 100.,  # Single dish diameter [m]
    'Tinst': 29. * (1e3),  # System temp. [mK]
    'survey_dnutot': 240.,  # Total bandwidth of *entire* survey [MHz]
    'survey_numax': 920.,  # Max. freq. of survey
    'dnu': 0.1,  # Bandwidth of single channel [MHz]
    'Sarea': 1e2 * (D2RAD) ** 2.,  # Total survey area [radians^2]
}
GBT.update(SURVEY)

Parkes = {
    'mode': 'dish',  # Interferometer or single dish
    'Ndish': 1,  # No. of dishes
    'Nbeam': 13,  # No. of beams (for multi-pixel detectors)
    'Ddish': 64.,  # Single dish diameter [m]
    'Tinst': 23. * (1e3),  # System temp. [mK]
    'survey_dnutot': 265.,  # Total bandwidth of *entire* survey [MHz]
    'survey_numax': 1420.,  # Max. freq. of survey
    'dnu': 0.1,  # Bandwidth of single channel [MHz]
    'Sarea': 1e2 * (D2RAD) ** 2.,  # Total survey area [radians^2]
}
Parkes.update(SURVEY)

GMRT = {
    'mode': 'dish',  # Interferometer or single dish
    'Ndish': 30,  # No. of dishes
    'Nbeam': 1,  # No. of beams (for multi-pixel detectors)
    'Ddish': 45.,  # Single dish diameter [m]
    'Tinst': 70. * (1e3),  # System temp. [mK]
    'survey_dnutot': 420.,  # Total bandwidth of *entire* survey [MHz]
    'survey_numax': 1420.,  # Max. freq. of survey
    'dnu': 0.1,  # Bandwidth of single channel [MHz]
    'Sarea': 1e2 * (D2RAD) ** 2.,  # Total survey area [radians^2]
}
GMRT.update(SURVEY)

BINGO = {
    'mode': 'dish',  # Interferometer or single dish
    'Ndish': 1,  # No. of dishes
    'Nbeam': 50,  # No. of beams (for multi-pixel detectors)
    'Ddish': 25.,  # Single dish diameter [m]
    # 'Ddish':            40.,               # Single dish diameter [m]
    'Tinst': 50. * (1e3),  # System temp. [mK]
    'survey_dnutot': 300.,  # Total bandwidth of *entire* survey [MHz]
    'survey_numax': 1260.,  # Max. freq. of survey
    'dnu': 0.1,  # Bandwidth of single channel [MHz]
    'Sarea': 2e3 * (D2RAD) ** 2.,  # Total survey area [radians^2]
    # 'Sarea':            5e3*(D2RAD)**2.,   # Total survey area [radians^2]
}
BINGO.update(SURVEY)

FAST = {
    'mode': 'dish',  # Interferometer or single dish
    'Ndish': 1,  # No. of dishes
    'Nbeam': 20,  # No. of beams (for multi-pixel detectors)
    'Ddish': 500.,  # Single dish diameter [m]
    'Tinst': 20. * (1e3),  # System temp. [mK]
    'survey_dnutot': 600.,  # Total bandwidth of *entire* survey [MHz]
    'survey_numax': 1000.,  # Max. freq. of survey
    'dnu': 0.1,  # Bandwidth of single channel [MHz]
    'Sarea': 2e3 * (D2RAD) ** 2.,  # Total survey area [radians^2]
}
FAST.update(SURVEY)

MeerKATb1 = {
    'mode': 'dish',  # Interferometer or single dish
    'Ndish': 64,  # No. of dishes
    'Nbeam': 1.,  # No. of beams (for multi-pixel detectors)
    'Ddish': 13.5,  # Single dish diameter [m]
    'Tinst': 29. * (1e3),  # System temp. [mK]
    'survey_dnutot': 435.,  # Total bandwidth of *entire* survey [MHz]
    'survey_numax': 1015.,  # Max. freq. of survey
    'dnu': 0.1,  # Bandwidth of single channel [MHz]
    # 'Sarea':            5e3*(D2RAD)**2.,  # Total survey area [radians^2]
    'Sarea': 25e3 * (D2RAD) ** 2.,  # Total survey area [radians^2]
    'n(x)': "array_config/nx_MKREF2_dec30.dat"  # Interferometer antenna density
}
MeerKATb1.update(SURVEY)

MeerKATb2 = {
    'mode': 'dish',  # Interferometer or single dish
    'Ndish': 64,  # No. of dishes
    'Nbeam': 1,  # No. of beams (for multi-pixel detectors)
    'Ddish': 13.5,  # Single dish diameter [m]
    'Tinst': 20. * (1e3),  # System temp. [mK]
    'survey_dnutot': 520.,  # Total bandwidth of *entire* survey [MHz]
    'survey_numax': 1420.,  # Max. freq. of survey
    'dnu': 0.1,  # Bandwidth of single channel [MHz]
    'Sarea': 5e3 * (D2RAD) ** 2.,  # Total survey area [radians^2]
    'n(x)': "array_config/nx_MKREF2_dec30.dat"  # Interferometer antenna density
}
MeerKATb2.update(SURVEY)

HIRAX = {
    'mode': 'interferom',  # Interferometer or single dish
    'Ndish': 1024,  # No. of dishes
    'Nbeam': 1,  # No. of beams (for multi-pixel detectors)
    'Ddish': 6.,  # Single dish diameter [m]
    'Tinst': 50. * (1e3),  # System temp. [mK]
    'Dmin': 6.,  # Min. interferom. baseline [m] from Ferreira et al. 2019
    'Dmax': 300.,   # Max. interferom. baseline [m] from Ferreira et al. 2019
    'survey_dnutot': 400.,
    'survey_numax': 800.,
    'dnu': 0.4,  # Bandwidth of single channel [MHz]
    'Sarea': 15e3 * (D2RAD) ** 2.,  # Total survey area [radians^2]
}
HIRAX.update(SURVEY)

CVlimited_z0to3 = {
    'mode':          'dish',       # Interferometer or single dish
    'Ndish':         1e10,         # No. of dishes (HUGE, to make sure CV liited!)
    'Nbeam':         1,            # No. of beams (for multi-pixel detectors)
    'Ddish':         100.,           # Single dish diameter [m]
    'Tinst':         50.*(1e3),    # System temp. [mK]
    'survey_dnutot': 1065.,         #1065., 
    'survey_numax':  1420.,        #1419.,  
    'dnu':           0.4,          # Bandwidth of single channel [MHz]
    'Sarea':         4.*np.pi,     # Total survey area [radians^2]
}
CVlimited_z0to3.update(SURVEY)

CVlimited_z0to5 = {
    'mode':          'dish',       # Interferometer or single dish
    'Ndish':         1e10,         # No. of dishes (HUGE, to make sure CV limited!)
    'Nbeam':         1,            # No. of beams (for multi-pixel detectors)
    'Ddish':         100.,           # Single dish diameter [m]
    'Tinst':         50.*(1e3),    # System temp. [mK]
    'survey_dnutot': 1184.,
    'survey_numax':  1420,
    'dnu':           0.4,          # Bandwidth of single channel [MHz]
    'Sarea':         4.*np.pi,     # Total survey area [radians^2]
}
CVlimited_z0to5.update(SURVEY)




