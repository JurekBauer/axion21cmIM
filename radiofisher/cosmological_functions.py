import numpy as np
from scipy import integrate

C = 3e5  # Speed of light, km/s


def hubble_value(z, H0, om_m, om_v=None, om_k=None, om_r=None):
    '''
    Calculates the redshift dependent hubble parameter via the equation H(z) = H0 E(z).
    The E function is defined below in a FLRW metric.
    '''
    if not om_v:
        # print('setting default value for omega_v, omega_v = 1 - omega_m')
        om_v = 1.0 - om_m
    if not om_k:
        # print('setting default value for omega_k, omega_k = 0')
        om_k = 0
    if not om_r:
        # print('setting default value for omega_r, omega_r = 0')
        om_r = 0
    return H0 * efunc(z, om_m, om_v, om_k, om_r)


def efunc(z, om_m, om_v, om_k, om_r):
    '''
    This returns E(z), defined as H(z) = H0 E(z).
    '''
    return (om_m * (1.0 + z) ** 3 + om_v + om_k * (1.0 + z) ** 2 + om_r * (1.0 + z) ** 4) ** 0.5


def comoving_distance(z, om_m, om_v=None, om_r=None, nsamples=1000):
    '''
    Returns the comoving distance in Mpc/h in a flat universe (Omega_k = 0).
    z: float
        redshift
    om_m: float
        matter density parameter
    om_v: float
        vacuum/lambda density parameter
    om_r: float
        radiation density parameter
    nsamples: integer
        nsamples to integrate over
    '''
    if not om_v:
        # print('setting default value for omega_v, omega_v = 1 - omega_m')
        om_v = 1.0 - om_m
    if not om_r:
        # print('setting default value for omega_r, omega_r = 0')
        om_r = 0
    _z = np.linspace(0., z, nsamples)
    r_c = integrate.simps(y=1.0 / efunc(z=_z, om_m=om_m, om_v=om_v, om_k=0., om_r=om_r), x=_z)
    r_comov = (C / 100.) * r_c
    return r_comov
