import astropy.constants as const
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

EPSILON = 2e-5


def radius_from_mass(M, omega_0):
    '''
    Converts halo mass to radius. Assumed is a sphere of density rho_0.
    '''
    return (3.0 * M / (4.0 * np.pi * omega_0 * rho_critical_at_redshift_zero())) ** (1.0 / 3.0)


def rho_critical_at_redshift_zero():
    '''
    Returns rho_critical(z=0) given in units solar_mass/h*(h/Mpc)^3 = solar_mass/Mpc^3*h^2.
    '''
    gravitational_constant = const.G.to('km**2*Mpc/(Msun*s**2)').value
    rho_critical = 3.0 * 100.0 ** 2 / (8. * np.pi * gravitational_constant)
    return rho_critical


def spherical_tophat_window_function(k, R, conditional_return=True):
    '''
    k needs to be given in units of h/Mpc and R in Mpc/h.
    Returns a dimensionless spherical tophat window function.
    R should be float and k can be an array.
    '''
    x = np.outer(R, k)
    if conditional_return:
        # takes care of numerical issues for x -> 0 (did expand the term and checked)
        return np.where(x > 1.4e-6, (3.0 / x ** 3) * (np.sin(x) - x * np.cos(x)), 1.0)
    else:
        return 3.0 / x ** 3 * (np.sin(x) - x * np.cos(x))


def some_term_wrt_spherical_tophat(k_arr, R, conditional_return=True):
    '''
    This function returns a term which is necessary for computing the derivative term dsigma^2/dM for the halo mass
    function. It is basically the derivative of the spherical top-hat function squared with respect to M, ignoring the
    constants in front of it.
    '''
    x = np.outer(R, k_arr)
    term1 = np.sin(x) * (1.0 - (3.0 / x ** 2)) + 3.0 / x * np.cos(x)
    term2 = np.sin(x) - x * np.cos(x)
    if conditional_return:
        # takes care of numerical issues for x -> 0 (did expand the term and checked)
        return np.where(x > 1e-3, term1 * term2, 0.0)
    else:
        return term1 * term2


def sigma_R(R, k_arr, PS_arr, plotting=False, label='', scipy_opt=True):
    '''
    Computes the variance of a power spectrum when smoothed by a spherical top-hat window function of size R.
    The inital power spectrum needs to be a 1d array in units (Mpc/h)^3,
    the wavenumber k is 1d array of same length in units (h/Mpc) and the size R needs to be a float in dimension Mpc/h.
    Returns the dimensionless quantity \sigma^2(R).
    '''
    if not scipy_opt:
        assert (abs(np.log(k_arr[1] / k_arr[0]) - np.log(k_arr[-1] / k_arr[-2])) < EPSILON), \
            '''
            Caution: You applied an integration scheme, which works on an equidistant array of ln(k). 
            However, the k_arr you passed has first ln difference k_arr: %.5f and last ln difference in k_arr: %.5f
            and difference %.4e
            ''' % (np.log(k_arr[1] / k_arr[0]), np.log(k_arr[-1] / k_arr[-2]),
               (abs(np.log(k_arr[1] / k_arr[0]) - np.log(k_arr[-1] / k_arr[-2]))))
    assert len(PS_arr) == len(k_arr), \
        '''
        k array and power spectrum array are not of same dimension, 
        but length(PS_arr) = %i and length(k_arr) = %i
        ''' % (len(PS_arr), len(k_arr))
    integrand_arr = PS_arr * spherical_tophat_window_function(k=k_arr, R=R) ** 2 * k_arr ** 2
    if plotting:
        plt.plot(k_arr, integrand_arr[0], label='integrand')
        plt.xlabel('$k$ [$h$/Mpc]')
        plt.ylabel(' [some units]')
        plt.title('Ingredients for $\sigma^2(R)$ for %s' % label)
        plt.legend()
        plt.show()
    if scipy_opt:
        sigma_squared = integrate.simps(y=integrand_arr, x=k_arr, axis=-1) / (2 * np.pi ** 2)
    else:
        dlnk = np.log(k_arr[1] / k_arr[0])
        sigma_squared = integrate.simps(y=integrand_arr * k_arr, dx=dlnk, axis=-1) / (2 * np.pi ** 2)
    return sigma_squared


def sigma_M(M, k_arr, PS_arr, omega_0, plotting=False, scipy_opt=True):
    '''
    Computes \sigma^2(M(R)) (dimensionless).
    M needs to be given as a float in solar mass/h, k_arr in h/Mpc and rho_matter in solar_mass/h*(h/Mpc)^3.
    '''
    R = radius_from_mass(M=M, omega_0=omega_0)
    if plotting:
        label = 'R = %.4f Mpc/h from M = %.2e solar_mass/h' % (R, M)
    else:
        label = ''
    return sigma_R(R=R, k_arr=k_arr, PS_arr=PS_arr, plotting=plotting, label=label, scipy_opt=scipy_opt)


def sheth_tormen_fit_func(v, A=0.3222, p=0.3, q=0.707):
    '''
    Computes Sheth-Tormen fit function from year 1999. A = 0.3222, p = 0.3 and q = 0.707 by default.
    v := delta_c/sigma(M), thus dimensionless by definition.
    '''
    return A * (2.0 * q / np.pi) ** 0.5 * v * (1.0 + (q ** 0.5 * v) ** (-2.0 * p)) * np.exp((-q * v ** 2.0) / 2.)


def sheth_tormen_halo_bias(sigma, p=0.3, q=0.707, delta_c=1.686):
    v = delta_c / sigma
    term1 = (q * v ** 2 - 1.0) / delta_c
    term2 = (2.0 * p / delta_c) / (1.0 + (q * v ** 2) ** p)
    bias = 1.0 + term1 + term2
    return bias


def integral_for_derivative_term_for_HMF(PS_arr, k_arr, omega_0, M, plotting=False, scipy_opt=True):
    if not scipy_opt:
        assert (abs(np.log(k_arr[1] / k_arr[0]) - np.log(k_arr[-1] / k_arr[-2])) < EPSILON), \
            '''
            Caution: You applied an integration scheme, which works on an equidistant array of ln(k). 
            However, the k_arr you passed has first ln difference k_arr: %.5f and last ln difference in k_arr: %.5f
            and difference %.4e
            ''' % (np.log(k_arr[1] / k_arr[0]), np.log(k_arr[-1] / k_arr[-2]),
               (abs(np.log(k_arr[1] / k_arr[0]) - np.log(k_arr[-1] / k_arr[-2]))))
    assert len(PS_arr) == len(k_arr), \
        '''
        k array and power spectrum array are not of same dimension, 
        but length(PS_arr) = %i and length(k_arr) = %i
        ''' % (len(PS_arr), len(k_arr))
    R = radius_from_mass(M=M, omega_0=omega_0)
    integrand_arr = some_term_wrt_spherical_tophat(k_arr, R) * PS_arr / k_arr ** 2
    if plotting:
        label = 'R = %.4f Mpc/h from M = %.2e solar_mass/h' % (R, M)
        plt.plot(k_arr, integrand_arr, label='integrand')
        plt.xlabel('$k$ [$h$/Mpc]')
        plt.title('''Integrand part except power spectrum of integral of derivative term
        for %s''' % label)
        plt.legend()
        plt.show()
    if scipy_opt:
        integral_result = integrate.simps(y=integrand_arr, x=k_arr, axis=-1)
    else:
        dlnk = np.log(k_arr[1] / k_arr[0])
        integral_result = integrate.simps(y=integrand_arr * k_arr, dx=dlnk, axis=-1)
    return integral_result


def halo_mass_function_defining_equation(omega_0, M_arr, k_arr, PS_arr, plot_integrands=False,
                                         delta_c=1.686, scipy_opt=True):
    '''
    Computes the halo mass function, dn/dln(M) = n(M) dm/dln(M) = n(M)*M in units of (h/Mpc)^3.
    Omega_0 is the matter/CDM+baryon/CDM density and needs to be a float object,
    M_arr is mass array given in solar_mass/h and k_arr needs to be given in units h/Mpc and PS_arr in (Mpc/h)^3.
    '''
    if not scipy_opt:
        assert (abs(np.log(k_arr[1] / k_arr[0]) - np.log(k_arr[-1] / k_arr[-2])) < EPSILON), \
            '''
            Caution: You applied an integration scheme, which works on an equidistant array of ln(k). 
            However, the k_arr you passed has first ln difference k_arr: %.6f and last ln difference in k_arr: %.6f
            ''' % (np.log(k_arr[1] / k_arr[0]), np.log(k_arr[-1] / k_arr[-2]))
    sigma_squared_arr = sigma_M(M=M_arr, k_arr=k_arr, PS_arr=PS_arr, omega_0=omega_0, plotting=False,
                                scipy_opt=scipy_opt)
    integral_for_derivative_term = integral_for_derivative_term_for_HMF(PS_arr=PS_arr, k_arr=k_arr,
                                                                        omega_0=omega_0, M=M_arr, plotting=False,
                                                                        scipy_opt=scipy_opt)
    R_arr = radius_from_mass(M=M_arr, omega_0=omega_0)
    v_arr = delta_c / np.sqrt(sigma_squared_arr)
    ST_function = sheth_tormen_fit_func(v_arr)
    dlnsigma2_dlnM = np.array(integral_for_derivative_term) * (3.0 / (sigma_squared_arr * np.pi ** 2 * R_arr ** 4))
    dn_dlnM = -0.5 * omega_0 * rho_critical_at_redshift_zero() / M_arr * ST_function * dlnsigma2_dlnM
    hmf_dictionary = {'f(v)': ST_function,
                      'v': v_arr,
                      'M': M_arr,
                      'sigma_squared': sigma_squared_arr,
                      'dn_dlnM': dn_dlnM,
                      'I_sigma2': np.array(integral_for_derivative_term),
                      'dlnsigma2_dlnM': dlnsigma2_dlnM}
    return hmf_dictionary
