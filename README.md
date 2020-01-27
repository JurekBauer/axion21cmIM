****21cm Intenstiy Mapping with axionCAMB****

This code was written for my Master thesis. It calculates the Fisher matrices for a given 21cm intensity mapping experiment, including the axion fraction as forecasted parameter.

To get this code running, you need to have a running `axionCAMB` version installed. Change the `CAMB_EXEC` to the directory of your `axionCAMB` executable. 
It is recommended to run this code with the `mpi4py`-package if one is interested to scan many surveys (and axion masses).  

The calculation depends on different parameters. They are grouped in the following way:  

  * *cosmological parameters*: $$A_s$$, $$n_s$$, $$k_{\mathrm{piv}}$$, $$h$$, $$\Omega_c$$, $$\Omega_b$$, $$m_a$$, $$\Omega_a$$, $$m_\nu$$, $$N_{\mathrm{eff}}$$ (and $$\Omega_\Lambda = 1 - (\Omega_c + \Omega_b + \Omega_a + \Omega_\nu)$$, since a flat universe is assumed), matter components to use for $$P_{\mathrm{lin}}$$ and $$\rho_0$$ corresponding to the axion mass choice  
  * *HI halo model parameters*: $$v_{c,0}$$, $$\beta$$, $$\gamma$$, $$c_{\hi,0}$$ (depending on the chosen halo mass function) and $$\Omega_\mathrm{HI}$$  
  * *experimental specifications*: survey, mode, total observation time  
  * *analysis specifications*: $$\Delta z$$, $$\ell$$ range, noise expression and for convergence: $$k$$-range and number of points in linear power spectrum, minimum and maximum halo mass and number of halo mass points.  

Generally, the code works as follows:  

  1. An input file needs to be passed to the `full_experiment_inifile.py`, containing all the relevant information: the fiducial cosmological and astrophysical parameters (except $$N_{\mathrm{eff}} = 3.046$$, which is hard-coded), as well as the experimental and analysis specifications (except $$\Delta z$$, which is hard-coded). 
This input file is read in as dictionaries, called `cosmo` for the cosmological parameters, while the experimental and analysis specifications are contained in the `expt` and `analysis_specifications` dictionary, respectively.  
  2. Secondly, the redshift range of the survey is divided into redshift bins with width $$\Delta z$$ and the Fisher matrix for each redshift bin is calculated. Since each redshift bin is independent of the other (in the Limber approximation), this procedure was allowed to be run in parallel.
  3. Thirdly, the resulting Fisher matrix for the redshift bin is saved, as well as the fiducial $$C_\ell$$ and $$\Delta C_\ell$$.  

The sketch visualizes the pipeline to calculate the Fisher matrix of one redshift bin. `halo_model` and `radiofisher` are modules computing the quantities inside the boxes.

<img src="code_sketch.png" alt="code_sketch" width="750"/>






 
