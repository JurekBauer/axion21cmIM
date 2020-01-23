import scipy.integrate
import scipy.interpolate
import radiofisher as rf
from HI_density_profile import *
from hydrogen_halo_mass_relation import *

verbos = 1


def mean_hydrogen_density(M_arr, dn_dlnM_arr, relation, plot_integrand=True, log_integration=True,
                          **relation_specific_args):
    '''
    Returns the mean hydrogen density \bar{\rho}_{HI} in M_sun/M_pc*h^2

    M_arr:                      Array of halo mass M in M_sun/h
    dn_dlnM_arr:                Array of halo mass function, dn_dlnM in h^3/Mpc^3, same dimension as M_arr
    relation:                   HI halo mass relation to use
    relation_specific_args:     Necessary additional arguments for the M_HI(M).
    '''
    if log_integration:
        assert (abs(np.log(M_arr[1] / M_arr[0]) - np.log(M_arr[-1] / M_arr[-2])) < 1e-6), \
            '''
            Caution: You applied an integration scheme, which works on an equidistant array of ln(k). 
            However, the k_arr you passed has first ln difference k_arr: %.5f and last ln difference in k_arr: %.5f
            ''' % (np.log(M_arr[1] / M_arr[0]), np.log(M_arr[-1] / M_arr[-2]))
    assert len(dn_dlnM_arr) == len(M_arr), \
        '''
        Length of dn_dlnM and M arrays are not equal. Something went wrong here. You need to get your problems solved! ;)
        '''
    M_HI_arr = hydrogen_halo_mass_relation(M=M_arr, equation=relation, **relation_specific_args)
    integrand_arr = M_HI_arr * dn_dlnM_arr
    if plot_integrand:
        plt.plot(M_arr, integrand_arr)
        plt.xscale('log')
        plt.xlabel('$M$ [$M_\\odot/h$]')
        plt.ylabel('Integrand $M_{\\rm{HI}}(M) n(M) M$ [$M_\\odot/h (h/$Mpc$)^3$]')
        plt.title('Integrand of $\\Omega_{\\rm{HI}}$, check for convergence')
        plt.show()
    if log_integration:
        dlnM = np.log(M_arr[1] / M_arr[0])
        mean_hydrogen_density = integrate.simps(y=integrand_arr, dx=dlnM)
    else:
        integrand_arr = integrand_arr / M_arr
        mean_hydrogen_density = integrate.simps(y=integrand_arr, x=M_arr)
    return mean_hydrogen_density


def PS_one_halo_term(k_arr, z, M_arr, dn_dlnM_arr, omega_0, omega_m, omega_lambda, hubble, c_HI0, gamma, relation,
                     log_integration=True, **relation_specific_args):
    '''
    Returns P_HI,1-halo in (Mpc/h)^3.

    k_arr:                      Array of Wavenumber k in h/Mpc
    z:                          Redshift z
    M_arr:                      Array of halo mass M in M_sun/h
    dn_dlnM_arr:                Array of halo mass function, dn_dlnM in h^3/Mpc^3, same dimension as M_arr
    omega_matter:               Omega_matter
    omega_lambda:               omega_lambda
    hubble:                     h
    c_HI0, gamma:               c_HI0 and gamma for the HI density profile
    relation:                   HI halo mass relation to use
    relation_specific_args:     necessary parameters for HI halo mass relation
    log_integration:            integrates over log-space
    '''
    if log_integration:
        assert (abs(np.log(M_arr[1] / M_arr[0]) - np.log(M_arr[-1] / M_arr[-2])) < 1e-6), \
            '''
            Caution: You applied an integration scheme, which works on an equidistant array of ln(k). 
            However, the k_arr you passed has first ln difference k_arr: %.5f and last ln difference in k_arr: %.5f
            ''' % (np.log(M_arr[1] / M_arr[0]), np.log(M_arr[-1] / M_arr[-2]))
    assert len(dn_dlnM_arr) == len(M_arr), \
        '''
        Length of dn_dlnM and M arrays are not equal. Something went wrong here. You need to get your problems solved! ;)
        '''
    M_HI_arr = hydrogen_halo_mass_relation(M=M_arr, equation=relation, **relation_specific_args)
    integrand_arr = M_HI_arr ** 2 * dn_dlnM_arr * \
                    np.abs(hydrogen_density_profile_kspace(k_arr=k_arr, M_arr=M_arr, z=z, om_0=omega_0, om_m=omega_m,
                                                           om_l=omega_lambda, c_HI0=c_HI0, gamma=gamma, h=hubble)) ** 2
    if log_integration:
        dlnM = np.log(M_arr[1] / M_arr[0])
        one_halo = integrate.simps(y=integrand_arr, dx=dlnM, axis=-1)
    else:
        integrand_arr = integrand_arr / M_arr
        one_halo = integrate.simps(y=integrand_arr, x=M_arr, axis=-1)
    mean_rho_HI = mean_hydrogen_density(M_arr=M_arr, dn_dlnM_arr=dn_dlnM_arr, relation=relation,
                                        plot_integrand=False, log_integration=True,
                                        **relation_specific_args)
    return one_halo / mean_rho_HI ** 2


def HI_bias(k_arr, z, M_arr, dn_dlnM_arr, sigma_arr, omega_0, omega_m, omega_lambda, hubble, c_HI0, gamma, relation,
            log_integration=True, **relation_specific_args):
    '''
    Returns b_HI(k) dimensionless. Exclude HI density profile info u_HI to calculate simply b_HI.

    k_arr:                      Array of Wavenumber k in h/Mpc
    z:                          Redshift z
    M_arr:                      Array of halo mass M in M_sun/h
    dn_dlnM_arr:                Array of halo mass function, dn_dlnM in h^3/Mpc^3, same dimension as M_arr
    omega_matter:               Omega_matter
    omega_lambda:               omega_lambda
    hubble:                     h
    c_HI0, gamma:               c_HI0 and gamma for the HI density profile
    relation:                   HI halo mass relation to use
    log_integration:
    relation_specific_args:
    '''
    if log_integration:
        assert (abs(np.log(M_arr[1] / M_arr[0]) - np.log(M_arr[-1] / M_arr[-2])) < 1e-6), \
            '''
            Caution: You applied an integration scheme, which works on an equidistant array of ln(k). 
            However, the k_arr you passed has first ln difference k_arr: %.5f and last ln difference in k_arr: %.5f
            ''' % (np.log(M_arr[1] / M_arr[0]), np.log(M_arr[-1] / M_arr[-2]))
    assert len(dn_dlnM_arr) == len(M_arr), \
        '''
        Length of dn_dlnM and M arrays are not equal. Something went wrong here. You need to get your problems solved! ;)
        '''
    M_HI_arr = hydrogen_halo_mass_relation(M=M_arr, equation=relation, **relation_specific_args)
    integrand_arr = M_HI_arr * dn_dlnM_arr * sheth_tormen_halo_bias(sigma=sigma_arr) * \
                    np.abs(hydrogen_density_profile_kspace(k_arr=k_arr, M_arr=M_arr, z=z, om_0=omega_0, om_m=omega_m,
                                                           om_l=omega_lambda, c_HI0=c_HI0,
                                                           gamma=gamma, h=hubble))
    mean_rho_HI = mean_hydrogen_density(M_arr=M_arr, dn_dlnM_arr=dn_dlnM_arr, relation=relation,
                                        plot_integrand=False, log_integration=True,
                                        **relation_specific_args)
    if log_integration:
        dlnM = np.log(M_arr[1] / M_arr[0])
        two_halo = integrate.simps(y=integrand_arr, dx=dlnM, axis=-1) / mean_rho_HI
    else:
        integrand_arr = integrand_arr / M_arr
        two_halo = integrate.simps(y=integrand_arr, x=M_arr, axis=-1) / mean_rho_HI
    return two_halo


def hydrogen_powerspectrum(k_arr, z, powerspec_dic, cosmo, analysis_specifications, output_opt=None):
    '''
    Calculates the hydrogen power spectrum.

    k_arr:          array of wavenumber k in h/Mpc
    z:              redshift z
    cachefile:      cachefile to load the power spectrum with which the derivative shall be calculated
    cosmo:          cosmo dictionary
    analysis_specifications
    output_opt:     option what the function shall return. Possible options: derivAs, derivns, vc0, beta
    '''
    # Sanity check. z needs to be a float, not an array e.g..
    assert isinstance(z, float), 'CAUTION: z is not a float! Type of z (redshift) is %s.' % type(z)

    # Create halo mass array (to be integrated over)
    M_arr = np.logspace(analysis_specifications['M_min'], analysis_specifications['M_max'],
                        analysis_specifications['N_mass_mesh'])
    # Load power spectrum with including the right components
    dict_key = rf.get_dict_key(cosmo['components_for_P'])  # powerspec_dic['dictionary_key']
    PS_arr = powerspec_dic['PS_' + dict_key]
    # Interpolate the power spectrum
    P_interpolated = scipy.interpolate.interp1d(powerspec_dic['k'], PS_arr,
                                                kind='quadratic',
                                                bounds_error=False,
                                                # fill_value='extrapolate')
                                                fill_value=0.)
    # Omega_0, depending on the components to be included
    omega_comps = powerspec_dic['omega_%s_0' % dict_key]
    # Get HMF dictionary, including the HMF and sigma
    hmf_dictionary = halo_mass_function_defining_equation(omega_0=omega_comps,
                                                          M_arr=M_arr,
                                                          k_arr=powerspec_dic['k'],
                                                          PS_arr=PS_arr,
                                                          plot_integrands=False, scipy_opt=False)
    # Calculate normalization of the M_HI(M) halo mass function, not necessary in principle!
    if cosmo['HI_halo_formula'] == 'Bull':
        if cosmo['HI_halo_formula_args']['f3'] == 1:
            pass
        else:
            print('\t\t\tf3 in cosmo dictionary was not set to 1 before assigning its correct value\n'
                  '\t\t\tin agreement with Omega_HI. Setting f3=1 and proceeding.')
            cosmo['HI_halo_formula_args']['f3'] = 1
        cosmo['HI_halo_formula_args']['redshift'] = z
        cosmo['HI_halo_formula_args']['h'] = cosmo['h']
        f3 = cosmo['omega_HI_0'] * rho_critical_at_redshift_zero() / \
             mean_hydrogen_density(M_arr=hmf_dictionary['M'], dn_dlnM_arr=hmf_dictionary['dn_dlnM'],
                                   relation=cosmo['HI_halo_formula'],
                                   plot_integrand=False, log_integration=True,
                                   **cosmo['HI_halo_formula_args'])
        cosmo['HI_halo_formula_args']['f3'] = f3
    elif cosmo['HI_halo_formula'] == 'PRA2017':
        cosmo['HI_halo_formula_args']['redshift'] = z
        cosmo['HI_halo_formula_args']['om_comps'] = omega_comps
        cosmo['HI_halo_formula_args']['om_b'] = cosmo['omega_b_0']
        cosmo['HI_halo_formula_args']['om_m'] = cosmo['omega_M_0']
        cosmo['HI_halo_formula_args']['om_lambda'] = cosmo['omega_lambda_0']
        cosmo['HI_halo_formula_args']['alpha'] = 1.0
        alpha = cosmo['omega_HI_0'] * rho_critical_at_redshift_zero() / \
                mean_hydrogen_density(M_arr=hmf_dictionary['M'], dn_dlnM_arr=hmf_dictionary['dn_dlnM'],
                                      relation=cosmo['HI_halo_formula'],
                                      plot_integrand=False, log_integration=True,
                                      **cosmo['HI_halo_formula_args'])
        print('alpha in HI halo mass relation = %f' % alpha)
        cosmo['HI_halo_formula_args']['alpha'] = alpha

    hydrogen_powerspec_one = PS_one_halo_term(k_arr=k_arr,
                                              z=z, M_arr=M_arr,
                                              dn_dlnM_arr=hmf_dictionary['dn_dlnM'],
                                              omega_lambda=cosmo['omega_lambda_0'],
                                              omega_m=cosmo['omega_M_0'],
                                              omega_0=omega_comps,
                                              hubble=cosmo['h'],
                                              c_HI0=cosmo['c_HI0'],
                                              gamma=cosmo['gamma_HI'],
                                              relation=cosmo['HI_halo_formula'],
                                              log_integration=True,
                                              **cosmo['HI_halo_formula_args'])

    bias_term_for_two = HI_bias(k_arr=k_arr,
                                z=z, M_arr=M_arr,
                                dn_dlnM_arr=hmf_dictionary['dn_dlnM'],
                                sigma_arr=np.sqrt(hmf_dictionary['sigma_squared']),
                                omega_lambda=cosmo['omega_lambda_0'],
                                omega_0=omega_comps,
                                omega_m=cosmo['omega_M_0'],
                                hubble=cosmo['h'],
                                c_HI0=cosmo['c_HI0'],
                                gamma=cosmo['gamma_HI'],
                                relation=cosmo['HI_halo_formula'],
                                log_integration=True,
                                **cosmo['HI_halo_formula_args'])
    hydrogen_powerspec_two = P_interpolated(k_arr) * bias_term_for_two ** 2
    if verbos >= 1:
        print('b_HI(k=1/r) = %f' % bias_term_for_two[0])
        print('b_HI(k=1000/r) = %f' % bias_term_for_two[-1])
    if output_opt == '1-halo':
        return hydrogen_powerspec_one
    elif output_opt == '2-halo':
        return hydrogen_powerspec_two
    elif output_opt == 'full':
        return hydrogen_powerspec_one + hydrogen_powerspec_two
    elif output_opt == 'full_plot':
        return hydrogen_powerspec_one, hydrogen_powerspec_two, hydrogen_powerspec_one + hydrogen_powerspec_two
    else:
        print('You specified an invalid output option for the function hydrogen_powerspectrum_derivatives().'
              'output_opt: %r' % output_opt)
        raise IOError
