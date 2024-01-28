#!/usr/bin/env python

"""

This script runs all individual scripts.
Individual scripts generate data and prepare plots

"""

from fig2b_MeOHox_Pt import make_plots_MeOHox_Pt
from figS2_model_parameter_estimate import make_plot_model_parameter_selectivity

if __name__ == "__main__":
    # figure 2b
    make_plots_MeOHox_Pt
    # figure S2
    make_plot_model_parameter_selectivity()


