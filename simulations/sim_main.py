#!/usr/bin/env python

"""
This script executes all part-scripts for publication-plot generation
Note: the numerical solution of the solver is not very efficient, the 
script may run up to 30 minutes for all solutions on a single core 

"""

import numpy as np


####################################################
#### (1) sampling pot and rghs and pH - Fig. 4a ####

from sim_current_selectivity import sample_polarization_curves, \
    plot_polarization_sel

msg = "run functions in sim_current_selectivity.py:\n"  \
    + "full multiscale simulations (Fig 4a)"
print(msg)

dat = sample_polarization_curves()
plot_polarization_sel('Fig4a_selectivity_pol', dat)

####################################################

print()

###############################################################
#### (2) concentration from simulations                    ####

from sim_concentrations import plot_conc_profile_ketene,\
    sample_local_pH_convection, plot_local_pH_convection, \
    get_surf_conc, plot_surf_concentrations

msg = "run function in sim_concentrations.py:\n" \
    + "concentration profiles, local concentrations (Fig 4cde, 5b, S10)"
print(msg)

#### (a) H2CCO concentration profile - rgh=15, pH=14       ####
d = sample_polarization_curves()
c = d[15][14]['Cnum'] # concentration data
x = np.linspace(0.0,1e-4,c[0][0,:].size) # x in Lx
plot_conc_profile_ketene('Fig5b_ketene_profile', x, c, d[15][14]['v'])

#### (b) plot local pH for different setups                ####
dat = sample_local_pH_convection()
plot_local_pH_convection('FigS10_local_pH_convection', dat)


#### (c) plot surface concentrations of H2CCO, OH, CO with ####
####     varying pH and roughness                          ####
csurf_kt, csurf_oh, csurf_co = get_surf_conc()
plot_surf_concentrations('Fig4cde_surface_concentrations', csurf_co, csurf_oh, csurf_kt)

#############################################################

print()

#################################################################
####  (3) model sensitivity                                  ####
from sim_model_sensitivity import sample_readsorption_energy, \
    plot_sel_vs_Gads, sample_rgh_pH, plot_selectivity_vs_rgh_pH,\
    sample_polarization_curves_facets, plot_polarization_facet

msg = "run function in sim_model_sensitivity.py:\n" \
    + "sensitivity to parameters and energyies (Fig S6, S11, S13)"
print(msg)

####  (a) sampling readsorption energy                       ####
engs = {100:{'de':0.36, 'ad':0.77},
        111:{'de':0.20, 'ad':1.05},
        110:{'de':0.37, 'ad':1.00},
        211:{'de':0.21, 'ad':0.52}}
dat = sample_readsorption_energy()
plot_sel_vs_Gads("FigS11_readsorption_vs_Gads", dat, {k:engs[k]['ad'] for k in engs})


####  (b) sampling rghs and pH vs. selectivity               ####
dph, drgh = sample_rgh_pH()
plot_selectivity_vs_rgh_pH('FigS6_selectivity_vs_rgh_pH', dph, drgh)


####  (c) exchange facet thermodynamics - keep SDS energies  ####
dat = sample_polarization_curves_facets()
plot_polarization_facet('FigS13a_selectivity_facet', dat,\
    adstag=False, legendtag=False)
dat = sample_polarization_curves_facets(include_ads=True)
plot_polarization_facet('FigS13b_selectivity_facet', dat,
    adstag=True, legendtag=True)
#################################################################


