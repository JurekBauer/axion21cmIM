'''
Load the input file into dictionaries used for the computation of the Fisher analysis
'''

def load_cosmology_input(file):
    cosmo_dic = {}
    analysis_dic = {}
    expt_dic = {}
    for line in file:
        line = line.split('=')
        if line[0] == '\n':
            continue
        elif line[0].strip() == 'noise':
            analysis_dic['noise'] = line[1].strip()
        elif line[0].strip() == 'experiment_number':
            expt_dic['k'] = eval(line[1].strip())
        elif line[0].strip() == 'ttot':
            expt_dic['ttot'] = eval(line[1].strip())
        elif line[0].strip() == 'M_min':
            analysis_dic['M_min'] = eval(line[1].strip())
        elif line[0].strip() == 'M_max':
            analysis_dic['M_max'] = eval(line[1].strip())
        elif line[0].strip() == 'N_mass_mesh':
            analysis_dic['N_mass_mesh'] = eval(line[1].strip())
        elif line[0].strip() == 'lmax':
            analysis_dic['lmax'] = eval(line[1].strip())
        elif line[0].strip() == 'lmin':
            analysis_dic['lmin'] = eval(line[1].strip())
        elif line[0].strip() == 'CAMB_kmax':
            analysis_dic['CAMB_kmax'] = eval(line[1].strip())
        elif line[0].strip() == 'transfer_k_perlogint':
            analysis_dic['transfer_k_perlogint'] = eval(line[1].strip())
        elif line[0].strip() == 'nonlinear_cutoff':
            analysis_dic['nonlinear_cutoff'] = eval(line[1].strip())        
        elif line[0].strip() == 'axfrac':
            cosmo_dic['axion_fraction'] = eval(line[1].strip())
            # axfrac = eval(line[1].strip())
        elif line[0].strip() == 'omega_d_0':
            cosmo_dic['omega_d_0'] = eval(line[1].strip())
            # omega_d_0 = eval(line[1].strip())
        elif line[0].split()[0] == '#':
            continue
        else:
            cosmo_dic[line[0].strip()] = eval(line[1].strip())
    return expt_dic, analysis_dic, cosmo_dic
