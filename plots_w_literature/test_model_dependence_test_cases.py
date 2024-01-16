#!/usr/bin/env python

"""

rerun test-case simulations and check for differences in mdls 

"""

from header import *


if __name__ == "__main__":
    datafile = "test_model_data_example.pkl"
    
    # run CO simulations
        
    # output: i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1
    dG_red_CO = 1.45; ddG_des = 0.165; cD_CO = 20.3e-10; cH_CO = 1100
    dGdes , dGads = get_adsdes_engs({'dG':0.1, 'ddG':ddG_des, 'H':cH_CO})
    args = dict(rdes=[dGdes], rads=[dGads], rred=[dG_red_CO], Ds=[cD_CO], Lxs=[50e-4])
    Us = [-1.6, -1.3]; rghs = [1, 50]
    
    dat_1 = sample_data(datafile, Us=Us, rghs=rghs, mdls=[1], **args)
    sel_1 = dat_1[:,13]/ dat_1[:,[13,14]].sum(axis=1)
    dat_2 = sample_data(datafile, Us=Us, rghs=rghs, mdls=[2], **args)
    sel_2 = dat_2[:,13]/ dat_2[:,[13,14]].sum(axis=1)
    
    # results are identical when converged
    print(sel_1)
    print(sel_2)

    # run Acetaldehyde simulations
    dG_red_CO = 0.96; ddG_des = 0.165; cD_Acdh = 13.75e-10; cH_Acdh = 0.0759
    dGdes , dGads = get_adsdes_engs({'dG':-0.14, 'ddG':ddG_des, 'H':cH_Acdh})
    args = dict(rdes=[dGdes], rads=[dGads], rred=[dG_red_CO], Ds=[cD_Acdh], Lxs=[50e-4])
    
    dat_1 = sample_data(datafile, Us=Us, rghs=rghs, mdls=[1], **args)
    sel_1 = dat_1[:,13]/ dat_1[:,[13,14]].sum(axis=1)
    dat_2 = sample_data(datafile, Us=Us, rghs=rghs, mdls=[2], **args)
    sel_2 = dat_2[:,13]/ dat_2[:,[13,14]].sum(axis=1)
    
    # results are identical when converged
    print(sel_1)
    print(sel_2)
