from hydrogen_halo_mass_relation import *


def hydrogen_density_profile_realspace(r, r_s, rho_0):
    return rho_0 * np.exp(-r / r_s)


def hydrogen_concentration_parameter(M, z, c_HI0, gamma, h):
    '''
    Returns HI concentration parameter c_HI(M,z) (dimensionless).
    M:                  Mass of dark matter halo in M_sun/h.
    z:                  Redshift
    c_HI0:              Initial HI concentration parameter, free parameter
    gamma:              "Free" parameter
    '''
    return c_HI0 * (M / (h * 1.0e11)) ** -0.109 * 4.0 / (1 + z) ** gamma


def hydrogen_density_profile_kspace(k_arr, M_arr, z, om_0, om_m, om_l, c_HI0, gamma, h):
    '''
    Return u_HI(k|M).
    k_arr:      Array of wavenumber k in h/Mpc
    M_arr:      Array of halo mass M in M_sun/h
    z:          Redshift z
    '''
    c_HI = hydrogen_concentration_parameter(M=M_arr, z=z, c_HI0=c_HI0, gamma=gamma, h=h)
    # get the virial radius in comoving coordinates in Mpc/h
    R_v = virial_radius_from_virial_mass(M_vir=M_arr, omega_0=om_0, omega_matter=om_m, omega_lambda=om_l, z=z,
                                         physical=False)
    r_s = R_v / c_HI
    a = 1.0 / (np.outer(k_arr, r_s))
    y = np.outer(k_arr, R_v)
    return (2.0 - np.exp(-c_HI) * (c_HI ** 2 + 2 * c_HI + 2.0)) ** (-1) * a ** 3 / ((a ** 2 + 1.0) ** 2) * \
           (np.exp(-c_HI) * (-a * (a ** 2 * y + a + y) * np.sin(y) - (a * (a * y + 2.0) + y) * np.cos(y) + np.sin(
               y)) + 2.0 * a)
    # Valid approximation for large c_HI
    # return 1.0/(1.0 + np.outer(k_arr, r_s) ** 2) ** 2
