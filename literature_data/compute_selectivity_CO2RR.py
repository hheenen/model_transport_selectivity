#!/usr/bin/env python

"""

This script performs the selectivity calculations for the CO2RR data:
CO (equation 18 in SI)
Ac (equation 19 in SI)
MeCHO (equation 20 in SI)

"""

import numpy as np
import sys
# plotting imports
sys.path.append("../.")
from model_Sel_CORAc.plotting import set_plotting_env, clrs, writefig
from scipy.constants import golden_ratio
import matplotlib.pyplot as plt


def load_data(filename):
    """
      helper function to load data from `literature_data`:
      use headers as keys in data-dict

    """
    # read numerical data
    dat = np.loadtxt(filename)[:,:-2] # don't load pre-computed values
    
    # read header --> 2nd line hardcoded, this would be easier using pandas...
    with open(filename, 'r') as f:
        lines = f.readlines()
        header = lines[1].strip("#").split()[:dat.shape[1]]

    # make dict for data-processing
    out = {}
    for i in range(dat[0,:].size):
        out.update({header[i]:dat[:,i]})

    return out


def compute_CO_selectivity(dat):
    """
      helper function to compute CO_selectivity
      dat = dictionary from load-data 

    """
    # non-FE keys
    okeys = ['U_SHE', 'I(mA/cm2)']
    # electron count of product in CO2RR
    nel = dict(CO=2, C2H4=10, EtOH=10, PropAldh=10, Prop=12, CH4=8, H2=2, \
        HCOO=2, Ac=6, H3CHO=8, C2H6=12, MeOH=6, HOH2CCHO=6, C2H4O2H2=8, \
        C3H5OH=10, C2H5CHO=10, Me2CO=10, MeCOCH2OH=8)
    # carbon count of product
    nC = dict(CO=1, C2H4=2, EtOH=2, PropAldh=3, Prop=3, CH4=1, H2=0, HCOO=1, \
        Ac=2, H3CHO=2, C2H6=2, MeOH=1, HOH2CCHO=2, C2H4O2H2=2, C3H5OH=3, C2H5CHO=3, Me2CO=3, MeCOCH2OH=3)
    # translating names
    tlm = {'Methane':'CH4', 'Ethylene':'C2H4', 'CO':'CO', 'Hydrogen':'H2','Ethane':'C2H6', 'Methanol':'MeOH',\
        'Formate':'HCOO', 'Ethanol':'EtOH', 'Acetate':'Ac', 'n-propanol':'Prop', 'Ac':'Ac', 'Acetaldehyde':'H3CHO',\
        'Glycolaldehyde':'HOH2CCHO', 'Ethylene-Glycol':'C2H4O2H2', 'Allyl-Alcohol':'C3H5OH', 'Propionaldehyde':'C2H5CHO',\
        'Acetone':'Me2CO', 'Hydroxyacetone':'MeCOCH2OH'}
    
    # now get molar (proportional) selectivity for CO vs. C1 and C2
    mkeys = [kk for kk in dat if kk not in okeys and kk not in ['Hydrogen', 'Formate']]
    # normalize to smth proportional to molar ratio
    # partial current
    jp = {kk:dat[kk] / nel[tlm[kk]] * nC[tlm[kk]] for kk in mkeys}
    # total CO+C1+C2 current
    jp_all = (np.array([jp[kk] for kk in jp]).T).sum(axis=1)
    sel_CO = jp['CO'] / jp_all
    return sel_CO


def compute_C2_selectivity(dat, Sads):
    """
      helper function to compute CO_selectivity
      dat = dictionary from load-data 

    """
    # electron count of products in CORR
    nel = {'Acetate':4, 'Ethanol':8, 'Ethylene':8, 'n-propanol':12, 'Acetaldehyde':6, \
        'Propionaldehyde':10, 'Ethane':10, 'Acetone':10, 'Allyl-Alcohol':10, \
        'Hydroxyacetone':6, 'Glycolaldehyde':4, 'Ethylene-Glycol':6}

    # C2+ products
    c2a = ['Ethylene', 'Acetate', 'Ethanol', 'Acetaldehyde', 'n-propanol',\
        'Allyl-Alcohol','Acetone','Ethane','Propionaldehyde','Allyl-alcohol',\
        'Hydroxyacetone','Ethylene-Glycol','Glycolaldehyde','Ethylene-glycol']

    # compute f_c2 (is only prop-to f_c2)
    f_c2 = np.zeros(dat['Ethanol'].shape)
    for ads in c2a:
        if ads in dat:
            f_c2 += dat[ads]/nel[ads]

    # compute selectivity
    sel = (dat[Sads] / nel[Sads]) / f_c2
    
    return sel


def show_sel_data(dsel, ylabel=r'selectivity (%)'):
    """
      function to plot selectivity data
      
      Parameters
      ----------
      dsel : dict of dicts
        contains selectivity data for various measurements
      ylabel : string
        ylabel
    
    """
    
    # setup plotting env
    height = 3.37/ golden_ratio *1.05 * 0.75
    set_plotting_env(width=3.37,height=height,
               lrbt=[0.16,0.98,0.35/height,0.98],fsize=9)
    
    fig, ax = plt.subplots()
    
    # prepare color list from clrs
    lclrs = list(clrs.values())
    fclrs = {k:lclrs[i] for i,k in enumerate(dsel.keys())}
    
    # plot # CONTINUE here
    for k in dsel:
        xy = dsel[k]
        ax.plot(xy[:,0], xy[:,1], marker='d', ls='-', color=fclrs[k], label=k)

    # legend and labels
    ax.legend(loc=1, prop={'size':5}, handlelength=1.0, frameon=False)
    ax.set_xlabel(r"U vs. SHE (V)")
    ax.set_ylabel(ylabel)
    
    plt.show()


if __name__ == "__main__":
    
    ########################
    ### CO selectivities ###
    ########################
    
    keys = ["Kanan-pc-Cu", "Kuhl-pc-Cu", "Huang-OD-Cu", \
                "Huang-111", "Huang-110", "Huang-100"]
    # iterate data and compute selectivities
    out = {}
    for key in keys:
        dat = load_data("dat_CO2RR_CO_%s.txt"%key)
        selCO = compute_CO_selectivity(dat)
        out.update({key:np.array([dat['U_SHE'], selCO]).T})
        
    show_sel_data(out, ylabel=r'selectivity CO / C$_1$ + C$_{2+}$ (%)')


    ########################
    ### Ac selectivities ###
    ########################

    keys = ["Cu-NP", "CuPd", "d-CuPd", "Cu3.4Pd", "Cu0.3Pd", \
        "CuAg_0%", "CuAg_0.25%", "CuAg_0.5%", "CuAg_1%"]

    # iterate data and compute selectivities
    out = {}
    for key in keys:
        dat = load_data("dat_CO2RR_CO_%s.txt"%key)
        selAc = compute_C2_selectivity(dat, Sads='Acetate')
        out.update({key:np.array([dat['U_SHE'], selAc]).T})

    show_sel_data(out, ylabel=r'selectivity Ac / C$_{2+}$ (%)')
    
    ##########################
    ### Acdh selectivities ###
    ##########################

    keys = ["COR@Cu-Flower", "COR@OD-Cu", "COR@pc-Cu"]

    # iterate data and compute selectivities
    out = {}
    for key in keys:
        dat = load_data("dat_CO2RR_CO_%s.txt"%key)
        selAcdh = compute_C2_selectivity(dat, Sads='Acetaldehyde')
        out.update({key:np.array([dat['U_SHE'], selAcdh]).T})

    show_sel_data(out, ylabel=r'selectivity Acdh / C$_{2+}$ (%)')

