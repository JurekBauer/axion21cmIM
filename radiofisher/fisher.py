import copy
import hashlib
import os
import uuid
import scipy.interpolate
from tempfile import gettempdir
import camb_wrapper as camb
import halo_model as hm
from cosmological_functions import *
from load_powerspec_data import *
from units import *

verbos = 0

EXP_OVERFLOW_VAL = 250.  # Max. value of exponent for np.exp() before assuming overflow
INF_NOISE = 1e200 # Very large finite no. used to denote infinite noise

# Location of CAMB executable
CAMB_EXEC = '/home/jurek/Dropbox/axionCAMB'
# CAMB_EXEC = '/usr/users/jbauer5/bin/axionCAMB'


def window_function(z, z_min, z_max):
    assert z_max > z_min, \
        '''
        z_max is not bigger than z_min. Somethings wrong here.
        z_max = %.5f, z_min %.5f
        ''' % (z_max, z_min)
    return np.where(z_min <= z <= z_max, 1.0 / (z_max - z_min), 0.0)


def Cl_signal(ell, z, z_min, z_max, cosmo, cachefile, analysis_specifications):
    '''
    Computes C_\ell, given an array of ell values and redshift bin [z_min, z_max].

    ell: array of ell values; default from 1 to 1000
    z_min: minimum redshift of redshift bin
    z_max: maximum redshift of redshift bin
    cosmo: dictionary with all important cosmological and astrophysical parameters
    cachefile: name of the output from axionCAMB
    analysis_specifications: e.g. noise_expression, ...
    '''
    assert z == cosmo['z'], '''z in cosmo dictionary does not agree with z passed to this function.'''
    # Hubble value H(z) in (km * h)/(s * Mpc)
    HH = hubble_value(z=z, H0=cosmo['h'] * 100, om_m=get_omega_M(cosmo=cosmo)) / cosmo['h']
    # Comoving distance r(z) in Mpc/h
    rr = comoving_distance(z=z, om_m=get_omega_M(cosmo=cosmo))
    # Effective k_arr in h/Mpc
    k = (ell + 0.5) / rr
    # Load power spectrum dictionary from axionCAMB
    powerspec_dic = load_power_spectrum(cosmo=cosmo, analysis_specifications=analysis_specifications,
                                        cachefile=cachefile,
                                        force=False, force_load=False)
    # Get hydrogen power spectrum P_HI = P_1-halo + P_2-halo in (Mpc/h)^3
    P_HI = hm.hydrogen_powerspectrum(k_arr=k, z=z, powerspec_dic=powerspec_dic, cosmo=cosmo,
                                     analysis_specifications=analysis_specifications, output_opt='full')
    # Prefactor in (h/Mpc)^3, C is speed of light in km/s
    prefactor = (window_function(z=z, z_min=z_min, z_max=z_max) / rr) ** 2 * HH * (z_max - z_min) / C
    if verbos >= 2:
        plt.plot(ell, prefactor * P_HI)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$C_\ell$')
        plt.show()
        k_overview = np.logspace(-3, 3, 100)
        P_1h_overview, P_2h_overview, P_overview = hm.hydrogen_powerspectrum(k_arr=k_overview, z=z,
                                                                             powerspec_dic=powerspec_dic,
                                                                             analysis_specifications=analysis_specifications,
                                                                             cosmo=cosmo,
                                                                             output_opt='full_plot')
        plt.plot(k_overview, P_overview)
        plt.plot(k_overview, P_1h_overview, label='1-halo', ls='--')
        plt.plot(k_overview, P_2h_overview, label='2-halo', ls='--')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$k$ [1/Mpc]')
        plt.ylabel(r'$P_{\rm{HI}}(k)$ [Mpc$^3$]')
        plt.title(r'$z = %.4f$' % z)
        plt.legend()
        plt.show()
    return prefactor * P_HI


def dish_response(ell, z, cosmo, expt):
    """
    Dish beam factor, B, for the noise
    covariance of a single-dish mode instrument.
    """
    # Define parallel/perp. beam scales
    l = 3e8 * (1. + z) / (1e6 * expt['nu_line'])
    theta_fwhm = l / expt['Ddish']
    # Sanity check: Require that Sarea > Nbeam * (beam)^2
    if (expt['Sarea'] < expt['Nbeam'] * theta_fwhm ** 2 / (16.0 * np.log(2))):
        raise ValueError("Sarea is less than (Nbeam * beam^2)")

    # Single-dish experiment has Gaussian beams in perp. and par. directions
    B_tot = (ell * theta_fwhm) ** 2 / (16.0 * np.log(2))
    B_tot[np.where(B_tot > EXP_OVERFLOW_VAL)] = EXP_OVERFLOW_VAL
    invbeam2 = np.exp(B_tot)
    return invbeam2


def load_interferom_file(fname):
    """
    Load n(u) file for interferometer and return linear interpolation function.
    """
    x, _nx = np.genfromtxt(fname).T
    interp_nx = scipy.interpolate.interp1d(x, _nx, kind='linear', bounds_error=False, fill_value=1./INF_NOISE)
    return interp_nx


def Nl_noise(ell, z_min, z_max, expt, cosmo, analysis_specifications):
    '''
    Computes N_\ell (dimensionless) for the given survey specifications in dictionary expt.

    ell:            Array of ell
    z_min:          Minimum redshift of redshift bin
    z_max:          Maximum redshift of redshift bin
    expt:           Dictionary with necessary survey specifications.
    cosmo:          cosmo dictionary, necessary for T_b
    analysis_specifications: e.g. noise_expression, ...
    '''
    z = 0.5 * (z_min + z_max)  # Center of redshift bin
    Tb = mean_brightness_temperature(z=z, cosmo_dic=cosmo)
    print('\tmean brightness temperature T_b = %.3e mK' % Tb)
    nu = expt['nu_line'] / (1. + z)  # Frequency of 21cm line at center of redshift bin z (MHz)
    wavelength = 3e8 / (nu * 1e6)  # Wavelength (m)
    Tsky = 60e3 * (350. / nu) ** 2.5  # Temp. of sky (mK)
    Tsys = expt['Tinst'] + Tsky  # Temp. of system (mK); temp. of instrument (mK)
    numin = expt['nu_line'] / (1. + z_max)
    numax = expt['nu_line'] / (1. + z_min)
    dnu_bin = (numax - numin)  # dnu of redshift bin
    # Number of receiver polarisation channels (default is two) and dish effic.
    npol = expt['npol'] if 'npol' in expt.keys() else 2.
    effic = expt['effic'] if 'effic' in expt.keys() else 0.7
    if expt['mode'] == 'dish':
        if analysis_specifications['noise'] == 'bull':
            # Calculate base noise properties
            Vsurvey = expt['Sarea']  # * dnu_bin / expt['nu_line']
            noise = Tsys ** 2. * Vsurvey / (npol * expt['ttot'] * dnu_bin)
            Aeff = effic * 0.25 * np.pi * expt['Ddish'] ** 2. \
                if 'Aeff' not in expt.keys() else expt['Aeff']
            theta_b = wavelength / expt['Ddish']
            noise *= wavelength ** 4. / (Aeff ** 2. * theta_b ** 4.)
            noise *= 1. / (expt['Ndish'] * expt['Nbeam'])
            noise *= dish_response(ell=ell, z=z, cosmo=cosmo, expt=expt)
            noise /= Tb ** 2
        elif analysis_specifications['noise'] == 'padmanabhan':
            # dnu_bin = expt['dnu']                     # dnu of redshift bin
            theta_b = wavelength / (expt['Ddish'] * expt['Ndish'])
            sigma2_b = theta_b ** 2 / (8 * np.log(2))  # \sigma^2_{\rm{beam}}
            # Set together all the relevant parts
            sigma2_pix = Tsys ** 2. / (npol * expt['ttot'] * dnu_bin)  # sigma_^2_pix
            Wl = np.exp(- ell ** 2 * sigma2_b)
            omega_pix = theta_b ** 2  # Omega_pix
            noise = (sigma2_pix / Tb ** 2) * omega_pix / Wl
        else:
            theta_b = wavelength / (expt['Ddish'])
            sigma2_b = theta_b ** 2 / (8.0 * np.log(2))  # sigma^2_beam
            # Set together all the relevant parts
            sigma2_pix = Tsys ** 2. / (npol * expt['ttot'] * dnu_bin * expt['Ndish'] * expt['Nbeam'])  # sigma_^2_pix
            B_tot = ell ** 2 * sigma2_b
            B_tot[np.where(B_tot > EXP_OVERFLOW_VAL)] = EXP_OVERFLOW_VAL
            Wl = np.exp(B_tot)
            if analysis_specifications['noise'] == 'chen':
                noise = (sigma2_pix / Tb ** 2) * Wl * 4. * np.pi
            elif analysis_specifications['noise'] == 'knox':
                noise = (sigma2_pix / Tb ** 2) * Wl * expt['Sarea']
            elif analysis_specifications['noise'] == 'CVlimited':
                print('\tNl_noise(): Cosmic variance limited survey. Setting noise to 0.')
                noise = np.array([0.] * len(ell))
            else:
                noise = np.array([0.] * len(ell))
        if analysis_specifications['nonlinear_cutoff']:
            print('\tNl_noise(): Cut-off non-linear scales.')
            k_nl = 0.14 / cosmo['h'] * (1.0 + z) ** (2.0/(2.0 + cosmo['ns']))
            ell_nl = k_nl * comoving_distance(z=z, om_m=get_omega_M(cosmo=cosmo))
            noise[np.where(ell > ell_nl)] = INF_NOISE
        return noise
    elif expt['mode'] == 'interferometric':
        if "n(x)" in expt.keys():
            # Rescale n(x) with freq.-dependence
            print("\tUsing user-specified baseline density, n(u)")
            x = (ell / (2. * np.pi))/ nu  # x = u / (freq [MHz]) = 2 pi ell / (freq [MHz])
            n_u = expt['n(x)'](x) / nu ** 2.  # n(x) = n(u) * nu^2
            n_u[np.where(n_u == 0.)] = 1. / INF_NOISE
        else:
            # Approximate expression for n(u), assuming uniform density in UV plane
            print("\tUsing uniform baseline density, n(u) ~ const.")
            u_min = expt['Dmin'] / wavelength
            u_max = expt['Dmax'] / wavelength

            # Sanity check: Physical area of array must be > combined area of dishes
            ff = expt['Ndish'] * (expt['Ddish'] / expt['Dmax']) ** 2.  # Filling factor
            print("\tArray filling factor: %3.3f" % ff)
            if ff > 1.:
                raise ValueError(("Filling factor is > 1; dishes are too big to "
                                  "fit in specified area (out to Dmax)."))

            # Uniform density n(u)
            n_u = expt['Ndish'] * (expt['Ndish'] - 1.) * wavelength ** 2. * np.ones(ell.shape) \
                  / (2. * np.pi * (expt['Dmax'] ** 2. - expt['Dmin'] ** 2.))
            n_u[np.where(ell < u_min * 2.0 * np.pi)] = 1. / INF_NOISE
            n_u[np.where(ell > u_max * 2.0 * np.pi)] = 1. / INF_NOISE
        FOV = wavelength ** 2 / expt['Ddish']
        noise = Tsys ** 2 * FOV ** 2 / (npol * dnu_bin * expt['ttot'] * n_u * Tb ** 2)
        if analysis_specifications['noise'] == 'CVlimited':
            print('\tNl_noise(): Cosmic variance limited survey. Setting noise to 0.')
            noise = np.array([0.] * len(ell))
        if analysis_specifications['nonlinear_cutoff']:
            print('\tNl_noise(): Cut-off non-linear scales.')
            k_nl = 0.14 / cosmo['h'] * (1.0 + z) ** (2.0/(2.0 + cosmo['ns']))
            ell_nl = k_nl * comoving_distance(z=z, om_m=get_omega_M(cosmo=cosmo))
            noise[np.where(ell > ell_nl)] = INF_NOISE
        return noise
    else:
        print('\tNl_noise(): No other mode than \'dish\' and \'interferometric\' implemented, yet. Abort.')
        raise IOError


def deriv_wrapper(ell_arr, zmin, zmax, cosmo, analysis_specifications, deriv_opt, cachefile_fid,
                  Delta_om_b=0.004,
                  Delta_om_d=0.004,
                  Delta_h=0.01,
                  Delta_axfrac=0.005,
                  Delta_As=1.0e-13,
                  Delta_ns=0.0005,
                  Delta_vc0=0.01,
                  Delta_beta=0.003):
    '''
    Returns the numerically calculated derivative of C_ell wrt to some model parameter
    (units depending on model parameter).
    Could probably modularize and functionalize this function to shorten up the code.

    ell_arr:            Array of ell's
    zmin:               Minimum redshift in redshift bin
    zmax:               Maximum redshift in redshift bin
    cosmo:              Dictionary of the fiducial cosmology
    analysis_specifications: e.g. noise_expression, ...
    deriv_opt:          Option: Specifies the parameter the derivative is taken wrt.
    cachefile_fid:      Path to the cachefile of the fiducial cosmology.
    '''

    def return_derivative(cosmo_minus, cosmo_plus, cachefile_minus, cachefile_plus, Delta):
        Cl_minus = Cl_signal(ell_arr, z=z, z_min=zmin, z_max=zmax,
                             cosmo=cosmo_minus, cachefile=cachefile_minus,
                             analysis_specifications=analysis_specifications)
        Cl_plus = Cl_signal(ell_arr, z=z, z_min=zmin, z_max=zmax,
                            cosmo=cosmo_plus, cachefile=cachefile_plus,
                            analysis_specifications=analysis_specifications)
        return (Cl_plus - Cl_minus) / (2.0 * Delta), Delta

    z = 0.5 * (zmin + zmax)
    file_base = 'PS_dictionaries/cache_powerspec_ma{}_axfrac_fid{}_'.format(cosmo['ma'], cosmo['axion_fraction'])
    if deriv_opt == 'omega_b':
        cosmo_Bplus = copy.deepcopy(cosmo)
        cosmo_Bminus = copy.deepcopy(cosmo)
        cosmo_Bplus['omega_b_0'] += Delta_om_b
        cosmo_Bminus['omega_b_0'] -= Delta_om_b
        cachefile_plus = file_base + 'Bplus_Delta{}_z{}.npy'.format(Delta_om_b, z)
        cachefile_minus = file_base + 'Bminus_Delta{}_z{}.npy'.format(Delta_om_b, z)
        return return_derivative(cosmo_minus=cosmo_Bminus, cosmo_plus=cosmo_Bplus,
                                 cachefile_minus=cachefile_minus, cachefile_plus=cachefile_plus, Delta=Delta_om_b)
    elif deriv_opt == 'omega_d':
        cosmo_Dplus = copy.deepcopy(cosmo)
        cosmo_Dminus = copy.deepcopy(cosmo)
        cosmo_Dplus['omega_d_0'] += Delta_om_d
        cosmo_Dminus['omega_d_0'] -= Delta_om_d
        cachefile_plus = file_base + 'Dplus_Delta{}_z{}.npy'.format(Delta_om_d, z)
        cachefile_minus = file_base + 'Dminus_Delta{}_z{}.npy'.format(Delta_om_d, z)
        return return_derivative(cosmo_minus=cosmo_Dminus, cosmo_plus=cosmo_Dplus,
                                 cachefile_minus=cachefile_minus, cachefile_plus=cachefile_plus, Delta=Delta_om_d)
    elif deriv_opt == 'hubble':
        cosmo_Hplus = copy.deepcopy(cosmo)
        cosmo_Hminus = copy.deepcopy(cosmo)
        cosmo_Hplus['h'] += Delta_h
        cosmo_Hminus['h'] -= Delta_h
        cachefile_plus = file_base + 'Hplus_Delta{}_z{}.npy'.format(Delta_h, z)
        cachefile_minus = file_base + 'Hminus_Delta{}_z{}.npy'.format(Delta_h, z)
        return return_derivative(cosmo_minus=cosmo_Hminus, cosmo_plus=cosmo_Hplus,
                                 cachefile_minus=cachefile_minus, cachefile_plus=cachefile_plus, Delta=Delta_h)
    elif deriv_opt == 'axfrac':
        cosmo_Aplus = copy.deepcopy(cosmo)
        cosmo_Aminus = copy.deepcopy(cosmo)
        cosmo_Aplus['axion_fraction'] += Delta_axfrac
        cosmo_Aminus['axion_fraction'] -= Delta_axfrac
        cachefile_plus = file_base + 'Aplus_Delta{}_z{}.npy'.format(Delta_axfrac, z)
        cachefile_minus = file_base + 'Aminus_Delta{}_z{}.npy'.format(Delta_axfrac, z)
        return return_derivative(cosmo_minus=cosmo_Aminus, cosmo_plus=cosmo_Aplus,
                                 cachefile_minus=cachefile_minus, cachefile_plus=cachefile_plus, Delta=Delta_axfrac)
    elif deriv_opt == 'As':
        cosmo_Asplus = copy.deepcopy(cosmo)
        cosmo_Asminus = copy.deepcopy(cosmo)
        cosmo_Asplus['As'] += Delta_As
        cosmo_Asminus['As'] -= Delta_As
        cachefile_plus = file_base + 'Asplus_Delta{}_z{}.npy'.format(Delta_As, z)
        cachefile_minus = file_base + 'Asminus_Delta{}_z{}.npy'.format(Delta_As, z)
        return return_derivative(cosmo_minus=cosmo_Asminus, cosmo_plus=cosmo_Asplus,
                                 cachefile_minus=cachefile_minus, cachefile_plus=cachefile_plus, Delta=Delta_As)
    elif deriv_opt == 'ns':
        cosmo_nsplus = copy.deepcopy(cosmo)
        cosmo_nsminus = copy.deepcopy(cosmo)
        cosmo_nsplus['ns'] += Delta_ns
        cosmo_nsminus['ns'] -= Delta_ns
        cachefile_plus = file_base + 'nsplus_Delta{}_z{}.npy'.format(Delta_ns, z)
        cachefile_minus = file_base + 'nsminus_Delta{}_z{}.npy'.format(Delta_ns, z)
        return return_derivative(cosmo_minus=cosmo_nsminus, cosmo_plus=cosmo_nsplus,
                                 cachefile_minus=cachefile_minus, cachefile_plus=cachefile_plus, Delta=Delta_ns)
    elif deriv_opt == 'vc0':
        cosmo_vc0plus = copy.deepcopy(cosmo)
        cosmo_vc0minus = copy.deepcopy(cosmo)
        cosmo_vc0plus['HI_halo_formula_args']['vc0'] += Delta_vc0
        cosmo_vc0minus['HI_halo_formula_args']['vc0'] -= Delta_vc0
        return return_derivative(cosmo_minus=cosmo_vc0minus, cosmo_plus=cosmo_vc0plus,
                                 cachefile_minus=cachefile_fid, cachefile_plus=cachefile_fid, Delta=Delta_vc0)
    elif deriv_opt == 'beta':
        cosmo_betaplus = copy.deepcopy(cosmo)
        cosmo_betaminus = copy.deepcopy(cosmo)
        cosmo_betaplus['HI_halo_formula_args']['beta'] += Delta_beta
        cosmo_betaminus['HI_halo_formula_args']['beta'] -= Delta_beta
        return return_derivative(cosmo_minus=cosmo_betaminus, cosmo_plus=cosmo_betaplus,
                                 cachefile_minus=cachefile_fid, cachefile_plus=cachefile_fid, Delta=Delta_beta)
    else:
        print(
            'deriv_wrapper_numerically(): You specified a derivative option %r, which is not yet implemented (or it is a typo?).' % deriv_opt)
        raise IOError


def fisher_integrands(ell_arr, zmin, zmax, cosmo, expt, cachefile, analysis_specifications):
    '''
    Calculate derivatives of all parameters for the Fisher matrix and for all ell.
    Return list of arrays of these derivatives.
    Order: ln(As), ns, omega_b, omega_d, axfrac, h, vc0, beta

    ell_arr:        Array of ell
    zmin:           Minimum redshift for redshift bin
    zmax:           Maximum redshift for redshift bin
    cosmo:          cosmo dictionary
    analysis_specifications: e.g. noise_expression, ...
    expt:
    cachefile:
    Returns:
        DeltaCl_squared:     Array of \Delta C_\ell^2
        Cl:                  Array of C_\ell
        deriv_list:          List of arrays of the derivatives
        paramnames           List of parameter names to which the derivatives are taken (same order as deriv_list, obviously).f
    '''
    print('fisher_integrands(): Calculating C_ell and N_ell...')
    Cl = Cl_signal(ell_arr, z=0.5 * (zmin + zmax), z_min=zmin, z_max=zmax, cosmo=cosmo, cachefile=cachefile,
                   analysis_specifications=analysis_specifications)

    Nl = Nl_noise(ell_arr, z_min=zmin, z_max=zmax, expt=expt, cosmo=cosmo,
                  analysis_specifications=analysis_specifications)

    # Calculate sky fraction
    f_sky = expt['Sarea'] / (4 * np.pi)
    DeltaCl_squared = 2.0 / ((2.0 * ell_arr + 1.0) * f_sky) * (Cl + Nl) ** 2

    print('fisher_integrands(): Calculating derivative wrt to **As** numerically...')
    deriv_As_num, delta_As = deriv_wrapper(ell_arr=ell_arr, zmin=zmin, zmax=zmax,
                                           analysis_specifications=analysis_specifications, cosmo=cosmo,
                                           deriv_opt='As', cachefile_fid=cachefile)
    print('fisher_integrands(): Calculating derivative wrt to **ns** numerically...')
    deriv_ns_num, delta_ns = deriv_wrapper(ell_arr=ell_arr, zmin=zmin, zmax=zmax,
                                           analysis_specifications=analysis_specifications, cosmo=cosmo,
                                           deriv_opt='ns', cachefile_fid=cachefile)
    print('fisher_integrands(): Calculating derivative wrt to **omega_d** numerically...')
    deriv_d, delta_d = deriv_wrapper(ell_arr=ell_arr, zmin=zmin, zmax=zmax,
                                     analysis_specifications=analysis_specifications, cosmo=cosmo,
                                     deriv_opt='omega_d', cachefile_fid=cachefile)
    print('fisher_integrands(): Calculating derivative wrt to **axion fraction** numerically...')
    deriv_axfrac, delta_axfrac = deriv_wrapper(ell_arr=ell_arr, zmin=zmin, zmax=zmax,
                                               analysis_specifications=analysis_specifications, cosmo=cosmo,
                                               deriv_opt='axfrac', cachefile_fid=cachefile)
    print('fisher_integrands(): Calculating derivative wrt to **omega_b** numerically...')
    deriv_b, delta_b = deriv_wrapper(ell_arr=ell_arr, zmin=zmin, zmax=zmax,
                                     analysis_specifications=analysis_specifications, cosmo=cosmo,
                                     deriv_opt='omega_b', cachefile_fid=cachefile)
    print('fisher_integrands(): Calculating derivative wrt to **h** numerically...')
    deriv_h, delta_h = deriv_wrapper(ell_arr=ell_arr, zmin=zmin, zmax=zmax,
                                     analysis_specifications=analysis_specifications, cosmo=cosmo,
                                     deriv_opt='hubble', cachefile_fid=cachefile)
    if cosmo['HI_halo_formula'] == 'PRA2017':
        print('fisher_integrands(): Calculating derivative wrt to **vc0** numerically...')
        deriv_vc0, delta_vc0 = deriv_wrapper(ell_arr=ell_arr, zmin=zmin, zmax=zmax,
                                             analysis_specifications=analysis_specifications,
                                             cosmo=cosmo, deriv_opt='vc0', cachefile_fid=cachefile)
        print('fisher_integrands(): Calculating derivative wrt to **beta** numerically...')
        deriv_beta, delta_beta = deriv_wrapper(ell_arr=ell_arr, zmin=zmin, zmax=zmax,
                                               analysis_specifications=analysis_specifications,
                                               cosmo=cosmo, deriv_opt='beta', cachefile_fid=cachefile)
        deriv_list = [deriv_As_num * cosmo['As'],
                      deriv_ns_num,
                      deriv_b,
                      deriv_d,
                      deriv_axfrac,
                      deriv_h,
                      deriv_vc0,
                      deriv_beta]
        paramnames = ['ln_As', 'ns', 'omega_b_0', 'omega_d_0', 'axfrac', 'h', 'vc0', 'beta']
    else:
        print('fisher_integrands(): You specified a HI halo mass relation which isn\'t implemendet for \
        calculating the derivatives yet. Will assume fixed HI halo mass relation.')
        deriv_list = [deriv_As_num * cosmo['As'],
                      deriv_ns_num,
                      deriv_b,
                      deriv_d,
                      deriv_axfrac,
                      deriv_h]
        paramnames = ['ln_As', 'ns', 'omega_b_0', 'omega_d_0', 'axfrac', 'h']
    return DeltaCl_squared, Cl, deriv_list, paramnames


def integrate_fisher_elements(derivs, DeltaCl_squared_arr):
    """
    Construct a Fisher matrix by performing the integrals for each element and
    populating a matrix with the results.
    """
    Nparams = len(derivs)

    # Loop through matrix elements, performing the appropriate integral for each
    F = np.zeros((Nparams, Nparams))
    for i in range(Nparams):
        for j in range(Nparams):
            if j >= i:
                term = (derivs[i] * derivs[j]) / DeltaCl_squared_arr
                F[i, j] = np.sum(term)
                F[j, i] = F[i, j]
                if verbos >= 2:
                    plt.plot(term)
                    plt.title('Fisher matrix')
                    plt.show()
    return F


def fisher_Cl(zmin, zmax, cosmo, expt, cachefile, analysis_specifications, survey_name=''):
    """
    Return Fisher matrix (and C_ell with errors) for given fiducial cosmology and experimental settings.

    Parameters
    ----------

    zmin, zmax : float
        Redshift window of survey

    cosmo : dict
        Dictionary of fiducial cosmological parameters

    expt : dict
        Dictionary of experimental parameters
    cachefile : str
        String where the cachefile for the power spectrum is saved for the fiducial cosmology
    analysis_specifications: dict
        Dictionary with analysis specifications necessary for calculating the Fisher matrix.
    survey_name : str
        Name of the survey. Passed to this function if one wants to plot something with that info on it.

    Returns
    -------

    F : array_like (2D)
        Fisher matrix for the parameters.

    """
    print('Survey: %s' % survey_name)
    # Copy, to make sure we don't modify input expt or cosmo
    cosmo = copy.deepcopy(cosmo)
    expt = copy.deepcopy(expt)

    # Get derivative terms for Fisher matrix integrands, then perform the integrals and populate the matrix
    print('Calculate derivative terms for Fisher matrix integrands...')

    ell_arr = np.arange(analysis_specifications['lmin'], analysis_specifications['lmax'] + 1, 1)
    # Consistency check: Is kmax in axionCAMB set higher than the maximum range which will be tested?
    rr = comoving_distance(z=0.5 * (zmax + zmin), om_m=get_omega_M(cosmo=cosmo))
    if verbos >= 1:
        print(
            "\tSanity check:\n\t\tCalculated P(k,z) from axionCAMB goes up to %3.2f h/Mpc, maximum l/r(z) is %3.2f h/Mpc." %
            (analysis_specifications['CAMB_kmax'], np.max(ell_arr) / rr))
    if analysis_specifications['CAMB_kmax'] < np.max(ell_arr) / rr:
        print("\tWARNING: Calculated P(k,z) from axionCAMB goes up to %3.2f h/Mpc, but maximum l/r(z) is %3.2f h/Mpc." %
              (analysis_specifications['CAMB_kmax'], np.max(ell_arr) / rr))
    # Calculate Delta C_ell and the derivatives wrt to the parameters
    deltacl_arr, Cl, derivs, paramnames = fisher_integrands(ell_arr=ell_arr, zmin=zmin, zmax=zmax, cosmo=cosmo,
                                                            expt=expt,
                                                            cachefile=cachefile,
                                                            analysis_specifications=analysis_specifications)
    print('Done.')
    F = integrate_fisher_elements(derivs=derivs, DeltaCl_squared_arr=deltacl_arr)
    # Return results
    return F, paramnames, deltacl_arr, Cl, ell_arr


def load_power_spectrum(cosmo, analysis_specifications, cachefile,
                        force=False, force_load=False):
    """
    Load power spectrum from cachefile and return dictionary with powerspectrum, different transfer functions and other
    axionCAMB output information.

    Parameters
    ----------

    cosmo : dict
        Dictionary of cosmological parameters
    analysis_specifications: dict
        Dictionary of analysis specifications
    cachefile : string
        Path to a cached CAMB matter powerspectrum output file.

    Returns
    -------
    powerspec_dic: Dictionary of important axionCAMB output
    """
    # Set-up axionCAMB parameters
    if verbos >= 1: print('\tSet up axionCAMB parameters.')
    p = convert_to_camb(cosmo)
    p['transfer_kmax'] = analysis_specifications['CAMB_kmax'] # CAMB_KMAX / cosmo['h'] if kmax <= 14. else 14.
    p['transfer_high_precision'] = 'T'
    p['transfer_k_per_logint'] = analysis_specifications['transfer_k_perlogint']

    # Check for massive neutrinos
    if verbos >= 1: print('\tCheck for massive neutrinos and parameters to axionCAMB dictionary.')
    if cosmo['mnu'] > 0.001:
        p['omnuh2'] = cosmo['mnu'] / 93.14
        p['massless_neutrinos'] = 2.046
        p['massive_neutrinos'] = '1'  # "2 1"
        p['nu_mass_eigenstates'] = 1

    print("\tload_power_spectrum(): Loading matter P(k)...")
    powerspec_dic = cached_camb_output(p=p, cachefile=cachefile, cosmo=cosmo,
                                       force=force, force_load=force_load)
    sigma8_computed = np.sqrt(hm.sigma_R(R=8.0, k_arr=powerspec_dic['k'], PS_arr=powerspec_dic['PS_total'],
                                         scipy_opt=False))
    print('\t\tsigma8 with PS_total from transfer function:\t%f' % sigma8_computed)
    return powerspec_dic


def mean_brightness_temperature(z, cosmo_dic, formula=''):
    """
    Brightness temperature Tb(z), in mK. Several different expressions for the
    21cm line brightness temperature are available:

     * 'santos': obtained using a simple power-law fit to Mario's data.
       (Best-fit Tb: 0.1376)
     * 'powerlaw': simple power-law fit to Mario's updated data (powerlaw M_HI
       function with alpha=0.6) (Default)
     * 'chang': from Chang et al.
    """
    omegaHI = cosmo_dic['omega_HI_0']
    if formula == 'santos':
        Tb = 0.1376 * (1. + 1.44 * z - 0.277 * z ** 2.)
    elif formula == 'powerlaw':
        Tb = 5.5919e-02 + 2.3242e-01 * z - 2.4136e-02 * z ** 2.
    elif formula == 'chang':
        Tb = 0.3 * (omegaHI / 1e-3) * np.sqrt(0.29 * (1. + z) ** 3.)
        Tb *= np.sqrt((1. + z) / 2.5)
        Tb /= np.sqrt(get_omega_M(cosmo=cosmo_dic) * (1. + z) ** 3. + get_omega_lambda(cosmo=cosmo_dic))
    else:
        om = get_omega_M(cosmo=cosmo_dic)
        ol = get_omega_lambda(cosmo=cosmo_dic)
        E = hm.efunc(z=z, om_m=om, om_v=ol, om_k=0, om_r=0)
        Tb = 190. * cosmo_dic['h'] * omegaHI * (1. + z) ** 2. / E
    return Tb


def cached_camb_output(p, cachefile, cosmo,
                       force=False, force_load=False):
    """
    Load P(k) or T(k) from cache, or else use CAMB to recompute it.
    """
    if verbos >= 1: print('\tCreate hash of cosmo parameters passed to axionCAMB and the components of the \n '
                          '\ttransfer function used for calculating the power spectrum and halo mass function, etc..\n'
                          '\tThese parameters are:')
    # Create hash of input cosmo parameters
    m = hashlib.md5()
    keys = list(p.keys())
    keys.sort()
    for key in keys:
        if verbos >= 1:
            print('\t\t%s: %s' % (key, p[key]))
        if key is not "output_root":
            m.update(("%s:%s" % (key, p[key])).encode('utf-8'))
    in_hash = m.hexdigest()

    # Check if cached version exists; if so, make sure its hash matches
    if verbos >= 1: print('\tCheck if cached version exits and if so, make sure its hash matches.')
    try:
        # Get hash from header of file
        powerspec_dic = np.load(cachefile, allow_pickle=True).item()
        # powerspec_dic = np.load(cachefile.item()
        out_hash = powerspec_dic['hash']

        # Compare with input hash; quit if hash doesn't match (unless force=True)
        if in_hash != out_hash and not force and not force_load:
            print("\tcached_camb_output(): Hashes do not match; delete the cache file and run again to update.")
            print("\t\tFile: %s" % cachefile)
            print("\t\tHash in file:  %s" % out_hash)
            print("\t\tHash of input: %s" % in_hash)
            raise ValueError()

        # Check if recalculation has been forced; throw fake IOError if it has
        if force:
            print("\t\tcached_camb_output(): Forcing recalculation of P(k).")
            raise IOError
        # If everything was OK, try to load from cache file, then return
        print('\tcached_camb_output(): Could reload powerspec_dic and hashes matched.')
        dict_key = get_dict_key(cosmo['components_for_P'])
        PS_comp_string = 'PS_' + dict_key
        omega_dict_key_string = 'omega_%s_0' % get_dict_key(comps=cosmo['components_for_P'])
        if not PS_comp_string in powerspec_dic.keys():
            print('\tcached_camb_output(): Power spectrum for wished components not found. Try to recompute...')
            powerspec_dic = recompute_PS_from_transfer_dictionary(transfer_dic=powerspec_dic,
                                                                  comps=cosmo['components_for_P'])
            powerspec_dic[omega_dict_key_string] = get_omega_dict_key(cosmo)
        if not omega_dict_key_string in powerspec_dic.keys():
            print('\tcached_camb_output(): Found power spectrum for wished components. Update omega_comps...')
            powerspec_dic[omega_dict_key_string] = get_omega_dict_key(cosmo)
        return powerspec_dic
    except IOError:
        if verbos >= 1: print('\tNo cached file found which fits. Need to recompute.')
        pass  # Need to recompute
    except:
        raise

    # Set output directory to /tmp and check that paramfiles and PS_HMF_dictionaries directory exists
    if verbos >= 1:
        print('\tSet output directory to ~/tmp and check that paramfiles and output directory exists.')
    root = gettempdir() + "/"
    if not os.path.exists("paramfiles/"):
        os.makedirs("paramfiles/", exist_ok=True)
        print("\tcached_camb_output(): Created paramfiles/ directory.")
    if not os.path.exists("PS_dictionaries/"):
        os.makedirs("PS_dictionaries/", exist_ok=True)
        print("\tcached_camb_output(): Created PS_dictionaries/ directory.")

    # Generate unique filename, create parameter file, and run axionCAMB
    if verbos >= 1: print('\tGenerate unique filename, create parameter file, and run axionCAMB')
    fname = str(uuid.uuid4())
    if verbos >= 1: print('\t\tInput filename is %s.ini in paramfiles/' % fname)
    p['output_root'] = root + fname
    camb.camb_params("%s.ini" % fname, **p)
    output = camb.run_camb("%s.ini" % fname, camb_exec_dir=CAMB_EXEC)

    # Load requested datafile (matterpower or transfer) and calculate halo mass function, etc.
    if verbos >= 1: print('\tLoad axionCAMB output and calculate power spectrum (with wished components) \n'
                          '\thalo mass function, etc.')
    path_to_files = 'paramfiles/'
    ini_filename = '%s.ini' % fname

    powerspec_dic = calc_PS_from_camb_output(path_to_files=path_to_files,
                                             ini_filename=ini_filename,
                                             comps=cosmo['components_for_P'],
                                             normalization=1,
                                             test_PS_plotting=False)
    # Save data to file, adding hash to library
    if verbos >= 1: print(
        "\tcached_camb_output(): Saving powerspec_dictionary to %s including the hash for identification.\n"
        "\tpowerspec_dictionary includes all important quantities for calculating hmf and so on \n"
        "\t(e.g. transfer functions, power spectras, growth rate)" % cachefile)
    omega_dict_key_string = 'omega_%s_0' % get_dict_key(comps=cosmo['components_for_P'])
    powerspec_dic[omega_dict_key_string] = get_omega_dict_key(cosmo)
    powerspec_dic['hash'] = in_hash
    np.save(cachefile, powerspec_dic)
    return powerspec_dic


def convert_to_camb(cosmo):
    """
    Convert cosmological parameters to CAMB parameters.
    (N.B. CAMB derives Omega_Lambda from other density parameters)
    """
    p = {}
    p['transfer_redshift__1___'] = cosmo['z']
    p['hubble'] = 100. * cosmo['h']
    p['omch2'] = get_omega_cdm(cosmo=cosmo) * cosmo['h'] ** 2.
    p['ombh2'] = cosmo['omega_b_0'] * cosmo['h'] ** 2.
    p['omk'] = 0.0  # 1. - (omega_M + omega_Lambda)
    p['scalar_spectral_index__1___'] = cosmo['ns']
    p['scalar_amp__1___'] = cosmo['As']
    p['pivot_scalar'] = cosmo['k_piv']
    p['omaxh2'] = get_omega_ax(cosmo=cosmo) * cosmo['h'] ** 2. + 1e-10  # +1e-10 to make sure to have a non-zero value
    p['m_ax'] = 10.0 ** cosmo['ma']
    return p


def zbins_const_dnu(expt, cosmo, bins=None, dnu=None, initial_dz=None):
    """
    Return redshift bin edges and centroids for bins that are equally-spaced
    in frequency (nu). Will either split the full range into some number of
    bins, or else fill the range using bins of const. dnu set from an initial
    bin delta_z.
    """
    # Get redshift range
    zmin = expt['nu_line'] / expt['survey_numax'] - 1.
    zmax = expt['nu_line'] / (expt['survey_numax'] - expt['survey_dnutot']) - 1.
    numax = expt['nu_line'] / (1. + zmin)
    numin = expt['nu_line'] / (1. + zmax)

    # nu as a function of z
    nu_z = lambda zz: expt['nu_line'] / (1. + zz)
    z_nu = lambda f: expt['nu_line'] / f - 1.

    # Return bin edges and centroids
    if bins is not None:
        # Get certain no. of bins
        nubins = np.linspace(numax, numin, bins + 1)
        zbins = z_nu(nubins)
    else:
        # Get bins with const. dr from initial delta_z
        if dnu is None:
            dnu = nu_z(zmin + initial_dz) - nu_z(zmin)
        dnu = -1. * np.abs(dnu)  # dnu negative, to go from highest to lowest freq.

        # Loop through range, filling with const. dr bins as far as possible
        nu = numax
        zbins = []
        while nu > numin:
            zbins.append(z_nu(nu))
            nu += dnu
        zbins.append(z_nu(numin))

        # Check we haven't exactly hit the edge (to ensure no zero-width bins)
        if np.abs(zbins[-1] - zbins[-2]) < 1e-4:
            del zbins[-2]  # Remove interstitial bin edge
        zbins = np.array(zbins)

    zc = [0.5 * (zbins[i + 1] + zbins[i]) for i in range(zbins.size - 1)]
    return zbins, np.array(zc)


def zbins_const_dz(expt, dz=None):
    """
    Return redshift bin edges and centroids for bins that are equally-spaced
    in redshift (z). Will split the full range into some number of
    bins.
    """
    # Get redshift range
    zmin = expt['nu_line'] / expt['survey_numax'] - 1.
    zmax = expt['nu_line'] / (expt['survey_numax'] - expt['survey_dnutot']) - 1.

    dz = +1. * np.abs(dz)  # dz positive, to go from lowest to highest redshift
    z = zmin
    zbins = []
    while z < zmax:
        zbins.append(z)
        z += dz
    zbins.append(zmax)

    # Check we haven't exactly hit the edge (to ensure no zero-width bins)
    if np.abs(zbins[-1] - zbins[-2]) < 1e-4:
        del zbins[-2]  # Remove interstitial bin edge
    zbins = np.array(zbins)

    zc = [0.5 * (zbins[i + 1] + zbins[i]) for i in range(zbins.size - 1)]
    return zbins, np.array(zc)
