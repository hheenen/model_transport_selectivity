#!/usr/bin/env python

"""

This script performs an estimation of catalyst roughness on basis of
(total) current density

> It should be noted, that this analysis should only be applied to 
different catalysts that have been measured in the same reaction setup
to exclude sensitivity of current to reaction design, pH, etc.

"""

import os, sys
import numpy as np
from scipy.optimize import curve_fit
import yaml

# plotting imports
sys.path.append("../.")
from model_Sel_CORAc.plotting import set_plotting_env, clrs, writefig
from scipy.constants import golden_ratio
import matplotlib.pyplot as plt


def load_data():
    """
      helper function to load data from `literature_data`
      hardcoded to load CuPd alloy, CuAg alloy and Huang data

    """
    key_dict = dict(
        CORR_CuPd=["Cu-NP", "Cu3.4Pd", "Cu0.3Pd", "CuPd", "d-CuPd"],
        CORR_CuAg=["CuAg_0%", "CuAg_0.25%", "CuAg_0.5%", "CuAg_1%"],
        CO2RR_SCs=["Huang-111","Huang-OD-Cu", "Huang-110", "Huang-100"], 
    )
    
    out = {}
    for dset in key_dict:
        # load full raw data (last two entries are always post-processed)
        dset_dat = {k:np.loadtxt("dat_CO2RR_CO_%s.txt"%k)[:,:-2] \
            for k in key_dict[dset]}
        out.update({dset:dset_dat})
    return out
            

def limit_curr(dj, clim):
    """
      helper function to remove data points w. current densities above 
      a certain clim

    """
    for k in dj:
        ind_ll = np.where(dj[k][:,1] <= clim)[0]
        dj[k] = dj[k][ind_ll,:]
    return dj


def remove_HER_from_curr(dj):
    """
      helper function to HER current from total current
      this would leave only CO2RR current
      NOTE: hard-coded to assume that the 3rd entry in dj[key] is Hydrogen FE

    """
    for k in dj:
        # compute partial current by factoring-out the hydrogen FE
        fe = dj[k][:,2:]
        fe_tot = fe.sum(axis=1)
        fe_cor = fe[:,1:].sum(axis=1)
        j_cor = dj[k][:,1] * (fe_cor/fe_tot)
        dj[k] = np.array([dj[k][:,0], j_cor]).T
    return dj


def fit_and_refit_exp_rgh(dj):
    """
        helper function to fit experimental U vs. J data against 
        y = a * np.exp(b*x)
        --
        note that alternatively one could also fit against 
        ln y = ln a + a*b; but equally robust and leads to over-emphasis of 
        low-overpotential region that behaves funny for e.g. SC data
        --
        dj = dictionary with data-sets for different catalyst
    
    """
    
    # (1) determine bs (exp coefficient) from full fit
    y = lambda x, a, b: a * np.exp(x*b) # full fit function
    popt_full = {}
    for k in dj:
        popt, pcov = curve_fit(y, dj[k][:,0], dj[k][:,1], p0=[0.0, -3.0])
        popt_full.update({k:popt})
    

    # (2) refit with different variants of b
    popt_refit = {}
    for kf in popt_full: # iterate full fits
        # function to fit with fixed exponential b
        cb = popt_full[kf][1]
        y = lambda x, a: a * np.exp(x*cb) # partial fit function
        ca = {}
        for k in dj:
            popt, pcov = curve_fit(y, dj[k][:,0], dj[k][:,1])
            ca.update({k:[popt[0], cb]})

        popt_refit.update({kf:ca})

    return popt_full, popt_refit


def average_roughness_from_fits(popt_refit):
    """
        helper function to determine relative roughness from prefactor a
        popt_refit (output from fit_and_refit) contains fitting parameters 
        from individual refits
    
    """
    ks = list(popt_refit.keys()) # fixed order keys
    rr = []
    for k in ks:
        # choose only prefact from fits with consistent expoential parameters
        ca = popt_refit[k]
        ca = np.array([ca[k][0] for k in ca])
        # normalize to relative roughness for each fit individually 
        #   little sensitivity to normalization after averaging
        rrgh = ca/ca.min() 
        rr.append(rrgh)
    
    # take mean and std of relative roughness & make dictionary
    rr = np.array(rr)
    rrgh = np.mean(rr, axis=0)
    rrgh = {ks[i]:rrgh[i] for i in range(len(ks))}
    rrgh_std = np.std(rr, axis=0)
    rrgh_std = {ks[i]:rrgh_std[i] for i in range(len(ks))}
    
    return rrgh, rrgh_std


def plot_exp_fit(fname, dj, dpopt, verbose=True, save_file=True):
    """
      function to plot data of exponential fits  
      
      Parameters
      ----------
      fname : string
        name for plot
      dj : dict
        raw data dict{key-catalyst:[U_SHE, current]}
      dpopt : list of dicts
        list with dictionaries for fitting parameters
        can also be single dict for single plot
      verbose : bool
        option to print to std-out
      save_file : bool
        option to save plots
    """
    # check if single set or list of sets - adjust figure layout
    if type(dpopt) is dict:
        nrows = 1; ncols = 1
    else:
        nrows = int(np.ceil(len(dpopt)/2))
        ncols = 2
    
    # setup plotting env
    height = 3.37/ golden_ratio *1.05 * 0.75*nrows
    set_plotting_env(width=3.37,height=height,
               lrbt=[0.16,0.98,0.35/height,0.98],fsize=9)
    fig, ax = plt.subplots(nrows, ncols)
    if ncols > 1:
        ax = ax.flatten()
    else:
        dpopt = [dpopt]
        ax = [ax]

    # plot indiviual fits
    for i in range(len(dpopt)):
        _plot_fit(ax[i], dj, dpopt[i])

    # output
    if save_file:
        if fname[0] != '_':
            fname = '_'+fname
        writefig("plot_fit%s"%fname, folder='roughness_estimates', \
            write_pdf=False, write_png=True)
    
    elif verbose:
        plt.show()
    plt.close(fig)


def _plot_fit(ax, dj, dpopt):
    """
        helper function to plot set of data against y = a * np.exp(b*x) fits
        dj = dictionary with arrays [U_SHE, current density] per catalyst
        dpopt = dictionary with fitting paramters per catalyst
    
    """
    # prepare color list from clrs
    lclrs = list(clrs.values())
    fclrs = {k:lclrs[i] for i,k in enumerate(dj.keys())}
    
    # plot
    for k in dj:
        xy = dj[k]
        # exp data points from dj
        flabel = "{:} a={:2.1e}, b={:2.1f}".format(k, *dpopt[k])
        ax.plot(xy[:,0], xy[:,1], marker='d', ls='None', color=fclrs[k], label=flabel)
        # fit from parameters
        xval = np.linspace(xy[:,0].min(), xy[:,0].max(), 100)
        y = lambda x, a, b: a * np.exp(x*b)
        ax.plot(xval, y(xval, *dpopt[k]), color=fclrs[k])

    # legend and labels
    ax.legend(loc=1, prop={'size':5}, handlelength=1.5, frameon=False)
    ax.set_xlabel(r"U vs. SHE (V)")
    ax.set_ylabel(r"j (mA cm$^{-2}$)")
    

def print_pretty_table(fname, rrgh, rrgh_std, m, verbose=True, save_file=True):
    ''' 
      function to output data in a ordered way
      
      Parameters
      ----------
      fname : string
        name-tag for file
      rrgh : dict
        output from `average_roughness_from_fits`
      rrgh_std : dict
        output from `average_roughness_from_fits`
      verbose : bool
        option to print to std-out
      save_file : bool
        option to save plots
    
    '''
    # header
    so = "{:^58}\n".format("-"*58)
    so += "{:^58}\n".format(fname)
    so += "{:^58}\n".format("-"*58)
    so += "| {:^20} | {:^14} | {:^14} |\n".format("key", "rel rho", "%i x rel. rho"%m)
    for k in rrgh:
        val = rrgh[k]; std = rrgh_std[k]
        so += "| {:20} | {:5.1f} +/- {:4.1f} | {:5.1f} +/- {:4.1f} |\n".format(k, val, std, val*m, std*m)
    so += "{:^58}\n".format("-"*58)

    # print
    if verbose:
        print(so)

    # save in file
    if save_file:
        if fname[0] != '_':
            fname = '_'+fname
        with open("roughness_estimates/roughness_fit%s.txt"%fname,'w') as out:
            out.write(so)
        

def rgh_estimates(dj, fname="", m=1, verbose=True, save_files=True):
    ''' 
      main function to perform roughness estimates
      
      Parameters
      ----------
      dj : dict of arrays
        dict = {key-cat-data:[Us, J]}
      fname : string
        name for plot/output file
      m : integer
        scaling factor for relative roughness to fit to model
      verbose : bool
        whether to print and show plots
      save_files : bool
        option to save plots and txt files

    '''
    # perform full and partial exp fit to current 
    popt_full, popt_refit = fit_and_refit_exp_rgh(dj) 

    # plot to document fit 
    plot_exp_fit(fname, dj, [popt_full]+[popt_refit[k] for k in popt_refit], \
        verbose=verbose, save_file=save_files)

    # obtain relative roughness averages and stds
    rrgh, rrgh_std = average_roughness_from_fits(popt_refit)

    print_pretty_table(fname, rrgh, rrgh_std, m=m, verbose=verbose, \
        save_file=save_files)

    return rrgh, rrgh_std


if __name__ == "__main__":

    ### remove HER from total current option ###
    # we left this option in, as it makes in principle sense for 
    # alloy data (Pd/Ag perfrom HER) - but it cannot be consistently applied to 
    # all data-sets (no HER data or SC which suffer from transport effects)
    remove_HER = False # left this option in, however, data inconsist

    # hardcoded m-values
    ms = {"CORR_CuPd":30, "CORR_CuAg":80, "CO2RR_SCs":1}
    alloy_data = {}

    # load and iterate through data
    out = load_data()
    for data_identifier in ["CORR_CuPd", "CORR_CuAg", "CO2RR_SCs"]:
        dj = out[data_identifier]
        
        if remove_HER:
            # reduce data to CO2RR current
            dj = remove_HER_from_curr(dj)
        else:     
            # reduce data to total current
            dj = {k:dj[k][:,[0,1]] for k in dj}
        
        # if H-cell, limit current density to 10 mA/cm2
        # due to mass-transport limitations
        if data_identifier == "CO2RR_SCs":
            dj = limit_curr(dj, 10) 
        
        # perform roughness estimate
        rrgh, rrgh_std = rgh_estimates(dj, data_identifier, \
            m=ms[data_identifier], verbose=True, save_files=True)

        # save scaled alloy data for plotting
        if data_identifier.find("Pd") != -1 or \
            data_identifier.find("Ag") != -1:
            alloy_data.update({k:float(rrgh[k]*ms[data_identifier]) \
                for k in rrgh})

    # save alloy-data for plotting in yaml file
    with open('roughness_estimates/alloy_data.yml', 'w') as outfile:
        yaml.dump(alloy_data, outfile, default_flow_style=False)


