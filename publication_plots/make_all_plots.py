#!/usr/bin/env python

"""

This script runs all individual scripts.
Individual scripts generate data and prepare plots

"""

from fig2a_ORR_Pt import make_plots_ORR_Pt
from fig2b_MeOHox_Pt import make_plots_MeOHox_Pt
from fig3a_CO2RR_Cu_SelCO import make_plot_CO2RR_CO_pot
from fig3b_CORR_Cu_SelAcHO import make_plot_CO2RR_Acdh_pot
from fig4_CORR_CuPd_SelAc import make_plot_CO2RR_Ac_rgh, make_plot_CO2RR_Ac_pot
from figS2_model_parameter_estimate import make_plot_model_parameter_selectivity

if __name__ == "__main__":
    # figure 2a, S3a
    make_plots_ORR_Pt()
    # figure 2b, S3b
    make_plots_MeOHox_Pt()

    no_errorbars=False
    # figure 3a
    make_plot_CO2RR_CO_pot(no_errorbars)
    # figure 3b
    make_plot_CO2RR_Acdh_pot(no_errorbars)
    # figure 4, S4b
    make_plot_CO2RR_Ac_rgh(no_errorbars)
    # figure S4a
    make_plot_CO2RR_Ac_pot(no_errorbars)
    
    # figure S2
    make_plot_model_parameter_selectivity()


