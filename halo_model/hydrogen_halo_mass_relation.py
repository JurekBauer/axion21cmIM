from virialized_mass_velocity_radius import *


def padmanabhan_refregier_relation(M, om_b, om_m, om_comps, om_lambda, redshift, alpha=0.17, beta=-0.55,
                                   vc0=10 ** 1.57, primordial_helium_mass_fraction=0.24,
                                   print_params=False):
    '''
    Return M_HI(M) using the Padmanabhan relation in M_sun/h

    M:              Halo mass M in M_sun/h
    om_b:
    om_m:
    om_lambda:
    redshift:
    alpha:
    beta:
    vc0:            Minimum circular velocity in km/s
    primordial_helium_mass_fraction:
    print_params:
    '''
    if print_params:
        print('''
        Parameters passed to hydrogen-halo mass relation from Padmanabhan and Refregier:
        alpha = %f
        beta = %f
        vc0 = %f
        om_b = %f
        om_m = %f
        om_lambda = %f
        Y_p = %f
        z = %f
        ''' % (alpha, beta, vc0, om_b, om_m, om_lambda, primordial_helium_mass_fraction, redshift))
    cosmic_hydrogen_mass_fraction = (1.0 - primordial_helium_mass_fraction) * om_b / om_m
    vc = circular_velocity_from_virial_mass(M=M, omega_0=om_comps, omega_matter=om_m, omega_lambda=om_lambda,
                                            z=redshift)  # in km/s
    return alpha * cosmic_hydrogen_mass_fraction * M * (M / 10 ** 11) ** beta * np.exp(-(vc0 / vc) ** 3)


def bagla_et_al_relation(M, f3, redshift, h, vmin=30, vmax=200):
    def mass_from_vel(vel, z, h):
        return 1.0e10 * h * (vel / 60.0) ** 3.0 * ((1.0 + z) / 4.0) ** (-1.5)

    M_min = mass_from_vel(vel=vmin, z=redshift, h=h)
    M_max = mass_from_vel(vel=vmax, z=redshift, h=h)
    term = f3 * M / (1.0 + M / M_max)
    return np.where(M >= M_min, term, 0.0)


def hydrogen_halo_mass_relation(M, equation, **args):
    '''
    Returns the HI mass as a function of halo mass, M_HI(M), depending on the specific relation given, in M_sun/h

    M:              Halo mass, M, in M_sun/h
    equation:       Equation for M_HI(M) to use. Bull or PRA2017?
    args:           Additional arguments which need to be passed to the M_HI(M) function.
    '''
    print_params = False
    plot_Fig6 = False
    if equation == 'Bull':
        if print_params:
            print('''
        You specified the HI-halo mass relation from Bull\'s neutrino paper (based on Bagla et al. 2010).
        Your keyword arguments which will be passed to this function are: %r
        ''' % args)
        return bagla_et_al_relation(M, **args)
    elif equation == 'PRA2017':
        if print_params:
            print('''
        You specified the HI-halo mass relation from Padmanabhan, Refregier and Amara 2017 paper.
        Your keyword arguments which will be passed to this function are: %r
        ''' % args)
        if plot_Fig6:  # plots Fig6 of arXiv:1611:1611.06235
            M_arr = np.logspace(10, 15, 100)
            h = 0.69  # assuming h = 0.69 is fiducial value
            om_comps = args['om_comps']
            om_b = args['om_b']
            om_m = args['om_m']
            om_lambda = args['om_lambda']
            alpha = 0.09
            beta = -0.58
            vc0 = 10 ** 1.56
            colors = ['royalblue', 'g', 'tomato', 'c', 'mediumvioletred']
            f = plt.figure()
            f.set_size_inches(8, 8)
            ax = f.add_subplot(111)
            for i in range(0, 5):
                M_HI_arr = padmanabhan_refregier_relation(M=M_arr * h, om_b=om_b, om_comps=om_m, om_m=om_m,
                                                          om_lambda=om_lambda,
                                                          redshift=i, alpha=alpha, beta=beta, vc0=vc0)
                plt.plot(M_arr, M_HI_arr / h, label=r'$z = %i$' % i, c=colors[i], lw=2)
            plt.legend(loc='lower right')
            ax.set_xlim(1e10, 1e15)
            ax.set_ylim(1e6, 1e11)
            ax.set_xlabel(r'$M$ [$M_\odot$]')
            ax.set_ylabel(r'$M_{\mathrm{HI}}$ [$M_\odot$]')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.yaxis.set_ticks_position('both')
            ax.xaxis.set_ticks_position('both')
            ax.tick_params(direction='in')
            plt.show()
        return padmanabhan_refregier_relation(M, **args)
    else:
        print('''
        You specified an unknown hydrogen halo mass relation. Maybe you should implement it. 
        Currently implemented are \'Bull\' and \'PRA2017\'. 
        ''')
