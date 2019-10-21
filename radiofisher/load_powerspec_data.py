import matplotlib.pyplot as plt
import numpy as np


def primordial_PS(k, initial_amplitude, k_pivot, scalar_index, hubble):
    '''
    Returns the dimensionless primordial power spectrum.
    Input: Initial amplitude is dimensionless, k in h/Mpc and k_pivot in 1/Mpc.
    '''
    return initial_amplitude * (k * (hubble / 100.0) / k_pivot) ** (scalar_index - 1.0)


def transfer_to_PS(k, total_transfer, initial_amplitude, scalar_index, hubble, k_pivot, proportionality_factor=1):
    '''
    Transfer function is given in actual transfer function divided by k^2 with k in 1/Mpc,
    this returns the power spectrum in units Mpc^3/h^3.
    '''
    return total_transfer ** 2 * (k * hubble / 100) ** 4 * primordial_PS(k, initial_amplitude, k_pivot, scalar_index,
                                                                         hubble) * 2 * np.pi ** 2 / k ** 3 * proportionality_factor


def load_transfer_from_file(name, path_to_file, normalization):
    '''
    Loads transfer function output file: Total transfer function and k_transfer will be loaded automatically,
    specify other components (e.g. \'baryon\' to load the baryon transfer function).
    k will be in units h/Mpc and I assume the transfer functions to be dimensionless.
    '''
    components_dic = {'k': (0,), 'CDM': (1,), 'baryon': (2,), 'photon': (3,), 'massless neutrino': (4,),
                      'massive neutrino': (5,),
                      'axion': (6,), 'growth rate': (7,), 'total': (8,)}
    transfer_dic = {
        'k': np.loadtxt(path_to_file + name, unpack=True, usecols=components_dic['k']),
        'transfer_total': np.loadtxt(path_to_file + name, unpack=True, usecols=components_dic['total']) * normalization,
        'growth_rate': np.loadtxt(path_to_file + name, unpack=True, usecols=components_dic['growth rate']),
        'transfer_CDM': np.loadtxt(path_to_file + name, unpack=True, usecols=components_dic['CDM']),
        'transfer_baryon': np.loadtxt(path_to_file + name, unpack=True, usecols=components_dic['baryon']),
        'transfer_axion': np.loadtxt(path_to_file + name, unpack=True, usecols=components_dic['axion']),
        'transfer_massless-neutrino': np.loadtxt(path_to_file + name, unpack=True,
                                                 usecols=components_dic['massless neutrino']),
        'transfer_massive-neutrino': np.loadtxt(path_to_file + name, unpack=True,
                                                usecols=components_dic['massive neutrino']),
        # 'transfer_photon': np.loadtxt(path_to_file + name, unpack=True, usecols=components_dic['photon'])
    }
    return transfer_dic


def get_parameter_library(ini_filename, path_to_files):
    file = open(path_to_files + ini_filename).readlines()
    parameter_dic = {}
    for line in file:
        if 'transfer_interp_matterpower' in line:
            parameter_dic['transfer_interp_matterpower'] = line.split('=')[1].strip()
        elif 'scalar_spectral_index(1)' in line:
            parameter_dic['scalar_index'] = float(line.split('=')[1].strip())
        elif 'scalar_amp(1)' in line:
            parameter_dic['initial_amplitude'] = float(line.split('=')[1].strip())
        elif 'hubble' in line:
            parameter_dic['hubble'] = float(line.split('=')[1].strip())
        elif 'pivot_scalar' in line:
            parameter_dic['k_pivot'] = float(line.split('=')[1].strip())
        elif 'ombh2' in line:
            parameter_dic['omega_bh2'] = float(line.split('=')[1].strip())
        elif 'omch2' in line:
            parameter_dic['omega_ch2'] = float(line.split('=')[1].strip())
        elif 'omnuh2' in line:
            parameter_dic['omega_nuh2'] = float(line.split('=')[1].strip())
        elif 'omaxh2' in line:
            parameter_dic['omega_axh2'] = float(line.split('=')[1].strip())
        elif 'omk' in line:
            parameter_dic['omega_k'] = float(line.split('=')[1].strip())
        elif 'transfer_kmax' in line:
            parameter_dic['transfer_kmax'] = float(line.split('=')[1].strip())
        elif 'transfer_redshift(1)' in line:
            parameter_dic['z'] = float(line.split('=')[1].strip())
        elif 'output_root' in line:
            output_root = line.split('=')[1].strip()
            base_name = output_root + '_'
        elif 'transfer_filename(1)' in line:
            transfer_filename = line.split('=')[1].strip()
        elif 'transfer_matterpower(1)' in line:
            matterpower_filename = line.split('=')[1].strip()
    parameter_dic['transfer_filename'] = base_name + transfer_filename
    print('\thalo_model.get_parameter_library(): base_name = %s' % base_name)
    print('\thalo_model.get_parameter_library(): transfer_filename = %s ' % transfer_filename)
    parameter_dic['matterpower_filename'] = base_name + matterpower_filename
    return parameter_dic


def get_dict_key(comps):
    comps = list(comps)
    comps.sort(key=lambda v: v.upper())
    dict_key = ''
    for i in range(len(comps) - 1):
        dict_key += comps[i] + '+'
    dict_key += comps[-1]
    return dict_key


def get_omega_dict_key(cosmo):
    omega_dict_key = 0
    for i in cosmo['components_for_P']:
        if i == 'CDM':
            omega = cosmo['omega_cdm_0']
        elif i == 'baryon':
            omega = cosmo['omega_b_0']
        elif i == 'total':
            omega = cosmo['omega_M_0']
        elif i == 'axion':
            omega = cosmo['omega_ax_0']
        else:
            print('\tUnknown component %s specified for power spectrum.' % i)
            raise KeyError
        omega_dict_key += omega
    print('\tget_omega_dict_key(): omega_%s_0 = %1.4f' % (get_dict_key(cosmo['components_for_P']), omega_dict_key))
    return omega_dict_key


def recompute_PS_from_transfer_dictionary(transfer_dic, comps, test_PS_plotting=False):
    # define variables used later on
    comps = list(comps)
    comps.sort(key=lambda v: v.upper())
    initial_amplitude, scalar_index = transfer_dic['initial_amplitude'], transfer_dic['scalar_index']
    hubble, k_pivot = transfer_dic['hubble'], transfer_dic['k_pivot']
    om_bh2, om_ch2 = transfer_dic['omega_bh2'], transfer_dic['omega_ch2']
    om_axh2 = transfer_dic['omega_axh2']
    dict_key = get_dict_key(comps)
    try:
        # This needs to be adapted, if you want to specify other components than
        # ['total'], ['baryon', 'CDM'] or ['axion', 'baryon', 'CDM']
        if comps == ['baryon', 'CDM']:
            transfer_dic['transfer_' + dict_key] = (om_bh2 * transfer_dic['transfer_baryon'] + om_ch2 * transfer_dic[
                'transfer_CDM']) / (om_bh2 + om_ch2)
        elif comps == ['axion', 'baryon', 'CDM']:
            transfer_dic['transfer_' + dict_key] = (om_bh2 * transfer_dic['transfer_baryon'] + om_ch2 * transfer_dic[
                'transfer_CDM'] + om_axh2 * transfer_dic['transfer_axion']) \
                                                   / (om_bh2 + om_ch2 + om_axh2)
        elif comps == ['total']:
            pass
        else:
            print('Error: You specified a combination of components which is not implemented yet.')
            raise IOError
        k_transfer = transfer_dic['k']
        transfer_dic['PS_' + dict_key] = transfer_to_PS(k=k_transfer,
                                                        total_transfer=transfer_dic['transfer_' + dict_key],
                                                        initial_amplitude=initial_amplitude, scalar_index=scalar_index,
                                                        hubble=hubble, k_pivot=k_pivot)
        transfer_dic['PS_total'] = transfer_to_PS(k=k_transfer, total_transfer=transfer_dic['transfer_total'],
                                                  initial_amplitude=initial_amplitude, scalar_index=scalar_index,
                                                  hubble=hubble, k_pivot=k_pivot)
        if test_PS_plotting:
            plt.plot(k_transfer, transfer_dic['PS_' + dict_key], ls='-', label=dict_key)
            plt.plot(k_transfer, transfer_dic['PS_tot'], ls='--', label='total')
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel('$k$ [$h$ Mpc$^{-1}$]')
            plt.ylabel('$P(k)$')
            plt.title('Power spectrum')
            plt.legend()
            plt.show()
    except KeyError:
        print('\trecompute_PS_from_transfer_dictionary(): Necessary transfer components not in the given dictionary.\n'
              '\tI think you have to rerun axionCAMB and delete/move the cachefile.')
        raise KeyError
    return transfer_dic


def calc_PS_from_camb_output(path_to_files, ini_filename, comps, normalization, test_PS_plotting=False):
    parameter_dic = get_parameter_library(ini_filename=ini_filename, path_to_files=path_to_files)
    comps = list(comps)
    comps.sort(key=lambda v: v.upper())
    # assert that transfer function has same number of data points as power spectrum,
    # needed for check that PS generated by transfer fn fits PS given from axionCAMB.
    if test_PS_plotting:
        assert parameter_dic['transfer_interp_matterpower'] == 'F', \
            '''
            Matter power spectrum has been interpolated from transfer function and therefore both arrays do not have the
            same dimension.
            Set \'transfer_interp_matterpower = F\' and rerun axionCAMB to use this function.
            '''
    assert comps == ['baryon', 'CDM'] or comps == ['total'] or comps == ['axion', 'baryon', 'CDM'], \
        '''
        You have specified a list of components (%r) out of which the transfer function (and power spectrum) will be 
        constructed with, which is not yet implemented. You have the option to specify baryon+CDM or total. 
        Other components and combinations need to be implemented (specify the equation for the transfer function). 
        ''' % comps
    print('''
    -----------------------------------------------
    kmax for transfer function is set to %.1f . 
    A high enough value is important for convergence of the integration calcs.
    -----------------------------------------------''' % parameter_dic['transfer_kmax'])
    # define variables used later on
    transfer_filename, matterpower_filename = parameter_dic['transfer_filename'], parameter_dic['matterpower_filename']
    initial_amplitude, scalar_index = parameter_dic['initial_amplitude'], parameter_dic['scalar_index']
    hubble, k_pivot = parameter_dic['hubble'], parameter_dic['k_pivot']
    om_bh2, om_ch2 = parameter_dic['omega_bh2'], parameter_dic['omega_ch2']
    om_nuh2, om_axh2 = parameter_dic['omega_nuh2'], parameter_dic['omega_axh2']
    if transfer_filename[0] == '/':
        print('transfer_filename = %s' % transfer_filename)
        transfer_dic = load_transfer_from_file(transfer_filename, '', normalization)
    else:
        # transfer_dic = load_transfer_from_file(transfer_filename, '', normalization,
        #                                        *comps)
        path_to_transfer_file = path_to_files.replace(path_to_files.split('/')[-2] + '/', transfer_filename)
        print('path_to_transfer_file: %s' % path_to_transfer_file)
        print('transfer_filename: %s' % transfer_filename)
        transfer_dic = load_transfer_from_file(path_to_transfer_file, '', normalization)
    dict_key = get_dict_key(comps)
    # This needs to be adapted, if you want to specify other components than
    # ['total'], ['baryon', 'CDM'] or ['axion', 'baryon', 'CDM']
    if comps == ['baryon', 'CDM']:
        transfer_dic['transfer_' + dict_key] = (om_bh2 * transfer_dic['transfer_baryon'] + om_ch2 * transfer_dic[
            'transfer_CDM']) / (om_bh2 + om_ch2)
    elif comps == ['axion', 'baryon', 'CDM']:
        transfer_dic['transfer_' + dict_key] = (om_bh2 * transfer_dic['transfer_baryon'] + om_ch2 * transfer_dic[
            'transfer_CDM'] + om_axh2 * transfer_dic['transfer_axion']) \
                                               / (om_bh2 + om_ch2 + om_axh2)
    elif comps == ['total']:
        pass
    else:
        print('Error: You specified a combination of components which is not implemented yet.')
        raise IOError
    k_transfer = transfer_dic['k']
    parameter_dic['number_of_k_values'] = len(k_transfer)
    print('\tLength of k arr is %i.' % len(k_transfer))
    transfer_dic['PS_' + dict_key] = transfer_to_PS(k=k_transfer, total_transfer=transfer_dic['transfer_' + dict_key],
                                                    initial_amplitude=initial_amplitude, scalar_index=scalar_index,
                                                    hubble=hubble, k_pivot=k_pivot)
    transfer_dic['PS_total'] = transfer_to_PS(k=k_transfer, total_transfer=transfer_dic['transfer_total'],
                                              initial_amplitude=initial_amplitude, scalar_index=scalar_index,
                                              hubble=hubble, k_pivot=k_pivot)
    if test_PS_plotting:
        plt.plot(k_transfer, transfer_dic['PS_' + dict_key], ls='-', label=dict_key)
        plt.plot(k_transfer, transfer_dic['PS_tot'], ls='--', label='total')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('$k$ [$h$ Mpc$^{-1}$]')
        plt.ylabel('$P(k)$')
        plt.title('Power spectrum')
        plt.legend()
        plt.show()
    transfer_dic.update(parameter_dic)
    return transfer_dic
