from calculate_HMF import *
from cosmological_functions import *


def virial_radius_from_virial_mass(M_vir, omega_0, omega_matter, omega_lambda, z, physical):
    '''
    Returns the virial radius R_v in Mpc/h.

    M_vir:          Virial mass in M_sun/h
    physical:       Option to get virial radius in physical coordinates, if set True
    ...
    '''
    d = omega_matter * (1.0 + z) ** 3 / (efunc(z=z, om_m=omega_matter, om_v=omega_lambda, om_k=0, om_r=0) ** 2.) - 1.0
    Delta_v = 18 * np.pi ** 2 + 82 * d - 39 * d ** 2
    rho_0 = omega_0 * rho_critical_at_redshift_zero()
    Rvir_comov = (3 * M_vir / (4 * np.pi * Delta_v * rho_0)) ** (1.0 / 3.0)
    Rvir_phys = (3 * M_vir / (4 * np.pi * Delta_v * rho_0 * (1 + z) ** 3)) ** (1.0 / 3.0)
    if physical:
        return Rvir_phys
    else:
        return Rvir_comov


def circular_velocity_from_virial_mass(M, omega_0, omega_matter, omega_lambda, z):
    '''
    Returns the circular velocity from the virial halo mass in km/s.

    M:             Halo mass M in M_sun/h
    ...
    '''
    rvir_phys = virial_radius_from_virial_mass(M_vir=M, omega_0=omega_0, omega_matter=omega_matter,
                                               omega_lambda=omega_lambda, z=z, physical=True)
    term1 = const.G.to('km**2*Mpc/(Msun*s**2)').value * M / rvir_phys  # in km^2/s^2
    velc = term1 ** 0.5
    return velc
