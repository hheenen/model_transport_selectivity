#!/usr/bin/env python

"""

This script executes a diffusion model selectivity thingy

"""

# general imports
import os
import numpy as np
import matplotlib.pyplot as plt
import pickle
from copy import deepcopy

# package imports
from model_Sel_CORAc.simtools import potsweep_rgh_pH
from model_Sel_CORAc.mkm.model_surfrct import get_CuAc_facet_mkm
from model_Sel_CORAc.transport.flux_conversion import f2j, convert_TOF2flux

from model_Sel_DRAR.transport_iterative import solve_MS_X_analytical
from model_Sel_DRAR.mkm_variable import get_mkm
from model_Sel_DRAR.model_plots import _plot_map


def convert_TOF2current(pr, roughness, nel_1, nel_2):
    pr_cur = np.array([f2j(pr[0]*roughness, nel_1), f2j(pr[1]*roughness, nel_2)])
    return(pr_cur)


def make_mkm_simple(eng_des, eng_ads, eng_red):
    mkm = get_CuAc_facet_mkm(100) # load standard model to manipulate
    mkm.engs0[3] = [0.0, eng_des, eng_des-eng_ads]
    mkm.engs0[4][1] = eng_red
    return mkm


def make_mkm_simple_2(eng_des, eng_ads, eng_red, k_dRI):
    mkm = get_mkm(100, k_dRI) # load standard model to manipulate
    mkm.engs0[3] = [0.0, eng_des, eng_des-eng_ads]
    mkm.engs0[4][1] = eng_red
    return mkm
    

def run_model(i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness):
    pH = 14
    i_diss_m = {1:2, 2:1}
    mkin = make_mkm_simple_2(eng_des, eng_ads, eng_red, k_dRI=i_model)
    
    ts, ys, pr, p0_norm, p0 = solve_MS_X_analytical(mkin, U_SHE, pH, \
        Dx, Lx, i_diss=i_diss_m[i_model], i_rads=7, roughness=roughness)
    output = {'cov':ys[-1,:], 'prod':pr, 'conc':p0}
    return(output)


def load_pickle_file(filename):
    with open(filename, 'rb') as pickle_file:
        data = pickle.load(pickle_file)
    return(data)


def write_pickle_file(filename, data):
    with open(filename, 'wb') as output:
        pickle.dump(data, output)
        

def get_data(data, *args):
    key = '-'.join([str(arg) for arg in args])
    if key in data:
        return(data[key])
    else:
        return(None)
            

def update_data(data, output, *args):
    key = '-'.join([str(arg) for arg in args])
    data.update({key: output})


def load_datafile(datafile):
    # if path, load file
    if os.path.isfile(datafile):
        sim_data = load_pickle_file(datafile)
    else:
        sim_data = {}
        # create backup-folder if not there
        datapath = '/'.join(datafile.split('/')[:-1])
        if not os.path.isdir(datapath):
            os.mkdir(datapath)

    return sim_data 


def sample_data(datafile, rdes, rads, rred, Us, Ds, Lxs, rghs, mdls):

    datafile = "bk_pkl_files/%s"%datafile
    sim_data = load_datafile(datafile)

    out_data = []
    n_count = 0

    for i_model in mdls:
        for roughness in rghs:
            for Dx in Ds:
                for Lx in Lxs:
                    for U_SHE in Us:
                        for eng_des in rdes:
                            for eng_ads in rads:
                                for eng_red in rred:
                                    args = (i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness)
                                    output = get_data(sim_data, *args)
                                    if output == None:
                                        print('running:', args)
                                        output = run_model(*args)
                                        update_data(sim_data, output, *args)
                                        if n_count%50 == 0:
                                            write_pickle_file(datafile, sim_data)
                                        n_count += 1
                             
                                    out_data.append([*args, *output['cov'], *output['prod'], output['conc'][0]])
    # i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1
    write_pickle_file(datafile, sim_data)
    if n_count != 0:
        print("evaluated %i data points"%n_count)

    out_data = np.array(out_data)
    
    return(out_data)


def sort_and_plot_colormap_energy(out_data):
    # probably should go for panda's to handle data!
    # plot data 
    for i_model in [1,2]:
        ind = np.where(out_data[:,0] == i_model)[0]
        dat = out_data[ind,:]
        for U_SHE in np.unique(dat[:,4]):
            ind = np.where(dat[:,4] == U_SHE)[0]
            di = dat[ind,:]
            for rgh in np.unique(di[:,7]):
                ind = np.where(di[:,7] == rgh)[0]
                dii = di[ind,:]
                for eng_red in np.unique(dii[:,3]):
                    ind = np.where(dii[:,3] == eng_red)[0]
                    d = dii[ind,:]
                    if d.size > 0:
             
                        # need to take np.nan out, set == 0.0 in grid
                        ind2 = np.isnan(d[:,13]).nonzero()[0]
                        d[ind2,13] = 0; d[ind2,14] = 1.0
                     
                        # i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1
                        xy = np.around(d[:,[1,2]], 4)
                        tag = r"$\Delta G_{\mathrm{CPET}}$=%.2f eV, U$_{SHE}$=%.1f, rgh=%.1f"%(d[0,3], U_SHE, rgh)
                        sel = d[:,13]/ d[:,[13,14]].sum(axis=1)
                     
                        xl = r"$\Delta G_{des}$ (eV)"
                        yl = r"$\Delta G_{ads}$ (eV)"
                        zl = r"selectivity"

                        lxy = np.array([[0.0, 0.4], [0.275, 0.0]])
                     
                        _plot_map("sel_clrmap_engs_model%i_%.2f_%.1f_%.0f"%(i_model, d[0,3], U_SHE, rgh), \
                            xy=xy, z=sel, ylabel=yl, xlabel=xl, zlabel=zl, showmarkers=False, tag=tag, line=lxy)


def _reduce_fixed(dat, fi, d_fixed):
    ind = np.where(np.around(dat[:,fi],5) == d_fixed[fi])[0]
    d = dat[ind,:]
    return d


def plot_effective_model(filename, dat, ixy, fixed, ylabel, xlabel):
    d = deepcopy(dat)
    for fi in fixed:
        d = _reduce_fixed(d, fi, d_fixed)
    ind2 = np.isnan(d[:,13]).nonzero()[0]; d[ind2,13] = 0; d[ind2,14] = 1.0
    d[:,6] = (d[:,5]*1e5)/(d[:,6]*1e-2) # Dx --> Dx/Lx
    d[:,2] = d[:,1]-d[:,2]
    #labels = {2:r"\Delta"
    xy = d[:,ixy]
    sel = d[:,13]/ d[:,[13,14]].sum(axis=1)
    zl = r"selectivity"
  # tag = r"$\Delta G_{\mathrm{des}}$=0.2 eV, $\Delta G_{\mathrm{ads}}$=0.3 eV, D=2e-10 cm2/s"
    _plot_map(filename, \
        xy=xy, z=sel, ylabel=yl, xlabel=xl, zlabel=zl, showmarkers=False, tag="")#, line=np.array([[1.0, -1.6], [1.3, -1.2]]))
    

# TODO:
# (x) adjust model to change energy adsorption desorption barrier
# (x) add diffusion model (no solution reaction); analytical
# (x) RDS dependence --> current --> molar flux! 
# (x) replace mkm with simpler mkm
# (x) add mkm with non-desorbing step
# (x) add pkl file
# (x) screen selectivity
# (x) figure out the des, ads barrier weirdness
# (x) check model 2 really exactly the same?
# (x) make wrapper functions for screening
# (x) get relation for ads/des energy --> does NOT dependent on c_surf
# (x) screen E_red vs. U_SHE; find out relation --> fully linear! Clear!
# (x) screening of rgh, D, Lx
# (x) show Dx/Lx plot selectivity depends on ratio
# (x) make plots showing all unrelated parameters: RDS vs. U_SHE; rgh vs. U_SHE; K_ads-des vs. U_SHE; Dx/Lx vs. U_SHE
# (x) get experimental data (selectivities to compare)
# (x) correct guessed data for better fit
# (x) think of fancy 2D plots: rgh vs. U_SHE --> different curves for different products
# ( ) need to replace U_SHE with a neutral expression, i.e. dG(U) --> for cases put by hand potential
# NOTE: check if assymetry in ORR example different than in general example - first generalize the e-barrier


# ( ) Selectivity maps: tags on the bottom left or top listing fixed effective paramters; red tags to display relations
# ( ) Diffusion coefficient with 1e-9 in label
# ( ) selectivity bars limited to 1.0 or transfer to percentage
# ( ) understand the 0.4/0.275 facdtor seems to not change with any transport behavior (maybe cause it's linear)
#     --> derive selectivity expression
# ( ) final plot selectivity map should have a cut along red-line as inset! show barrier competition

# ( ) collect plots showing relations: U_SHE vs. dG_red; dG_des vs. dG_ads; Dx vs Lx
    

    ### --> plots for manuscript: dG_ads vs. dG_des color map with inset!
    ### --> Dx/Lx vs. rgh --> need higher res.
    ### --> rgh vs USHE + Dx/Lx vs USHE --> double plot
    ### 2D plots sel vs. USHE for different cases

    ### SI: USHE vs. dG_CPET; Lx vs Dx
 

# parameters k_ads/k_des, Diff, Ld, rho

if __name__ == "__main__":

    pH = 13 # dummy, not important

    # how and what to screen --> reduce parameters:
    # (1) energy cirlce --> eng_des vs. eng_ads at different eng_red, different U_SHE, different roughness
    #       check if always same eng_des/eng_ads relation
    # (2) check if relation eng_red / U_SHE at different roughnesses same? 
    # (3) check at relevant parameters influence of Dx, Lx, roughness

    eng_des = 0.2
    eng_ads = 0.3
    eng_red = 1.1

    U_SHE = -1.4
    Lx = 1e-4 # 1 micrometer in cm
    Dx = 20e-10 # in m2/s
    roughness = 1
    

    ################################
    # (1) energy cirlce --> eng_des vs. eng_ads at different eng_red, different U_SHE, different roughness
  # datafile = "model_data.pkl"
  # out_data = sample_data(datafile, rdes=np.arange(0.0, 0.5, 0.1), rads=np.arange(0.0, 0.5, 0.1), \
  #             rred=np.arange(1.0, 1.3, 0.1), Us=[-1.2, -1.4, -1.6], Ds=[Dx], Lxs=[Lx], rghs=[1,50], mdls=[1,2])

  # sort_and_plot_colormap_energy(out_data)

  # # test change of linearity with transport properties
  # datafile = "model_data_test.pkl"
  # d = sample_data(datafile, rdes=np.arange(0.0, 0.5, 0.1), rads=np.arange(0.0, 0.5, 0.1), \
  #             rred=[1], Us=[-1.4], Ds=[Dx], Lxs=[Lx*100], rghs=[1], mdls=[1])
  #                     
  # # need to take np.nan out, set == 0.0 in grid
  # ind2 = np.isnan(d[:,13]).nonzero()[0]
  # d[ind2,13] = 0; d[ind2,14] = 1.0
  # xy = np.around(d[:,[1,2]], 4)
  # sel = d[:,13]/ d[:,[13,14]].sum(axis=1)
  # xl = r"$\Delta G_{des}$ (eV)"; yl = r"$\Delta G_{ads}$ (eV)"; zl = r"selectivity"
  # _plot_map("sel_clrmap_engs_model%i_%.2f_%.1f_%.0f_Lxfct10"%(1, d[0,3], -1.4, 1), \
  #     xy=xy, z=sel, ylabel=yl, xlabel=xl, zlabel=zl, showmarkers=False, line=np.array([[0.0, 0.4], [0.275, 0.0]]))

    ################################
    #### (2) check if relation eng_red / U_SHE at different roughnesses same? 
  # # i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1
                                    
  # datafile = "model_data_potred.pkl"
  # d = sample_data(datafile, rdes=[0.2], rads=[0.3], \
  #         rred=np.arange(1.0, 1.3, 0.1), Us=[-1.2, -1.4, -1.6], Ds=[Dx], Lxs=[Lx], rghs=[1], mdls=[1])
  # xy = np.around(d[:,[3,4]], 4)
  # sel = d[:,13]/ d[:,[13,14]].sum(axis=1)
  # xl = r"$\Delta G_{CPET}$ (eV)"; yl = r"U$_{\mathrm{SHE}}$ (V)"; zl = r"selectivity"
  # tag = r"$\Delta G_{\mathrm{des}}$=0.2 eV, $\Delta G_{\mathrm{ads}}$=0.3 eV, rgh=1"
  # _plot_map("sel_clrmap_redpot_model1", \
  #     xy=xy, z=sel, ylabel=yl, xlabel=xl, zlabel=zl, showmarkers=False, tag=tag, line=np.array([[1.0, -1.6], [1.3, -1.2]]))

    
    ################################
    #### (3) check at relevant parameters influence of Dx, Lx, roughness
  # # i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1
    # typical range Diff-COeff: 8e-10 - 60e-10 # in m2/s 
    #       also tortuosity --> up to 0.2 --> 2e-10 - 60e-10
    # typical diffusion lengths 1-100 micrometer - 1e-4-1e-2

  # datafile = "model_data_transport.pkl"
  # dat = sample_data(datafile, rdes=[0.2], rads=[0.3], \
  #         rred=[1.1], Us=[-1.3], Ds=np.linspace(2,60,10)*1e-10, Lxs=np.logspace(-4,-2,10), rghs=[1,5,10,25,50], mdls=[1])
  # for rgh in [1,50]:
  #     ind = np.where(dat[:,7] == rgh)[0]; d = dat[ind,:]
  #     ind2 = np.isnan(d[:,13]).nonzero()[0]; d[ind2,13] = 0; d[ind2,14] = 1.0
  #     xy = d[:,[5,6]]
  #     sel = d[:,13]/ d[:,[13,14]].sum(axis=1)
  #     xl = r"$D_{eff}$ (m$^2$ s$^{-1}$)"; yl = r"$L_x$ (cm)"; zl = r"selectivity"
  #     tag = r"$\Delta G_{\mathrm{des}}$=0.2 eV, $\Delta G_{\mathrm{ads}}$=0.3 eV, rgh=%i"%rgh
  #     _plot_map("sel_clrmap_tr_D_Lx_rgh%i_model1"%rgh, \
  #         xy=xy, z=sel, ylabel=yl, xlabel=xl, zlabel=zl, showmarkers=False, tag=tag, \
  #         line=[np.array([[0.0, 0.0], [6e-9, 0.006]]), np.array([[0.0,0.0], [6e-9, 0.003]])])

  # ind = np.where(np.around(dat[:,5],11) == 2.13e-09)[0]; d = dat[ind,:]
  # ind2 = np.isnan(d[:,13]).nonzero()[0]; d[ind2,13] = 0; d[ind2,14] = 1.0
  # xy = d[:,[7,6]]
  # sel = d[:,13]/ d[:,[13,14]].sum(axis=1)
  # yl = r"$L_x$ (cm)"; xl = r"roughness"; zl = r"selectivity"
  # tag = r"$\Delta G_{\mathrm{des}}$=0.2 eV, $\Delta G_{\mathrm{ads}}$=0.3 eV, D=2e-10 cm2/s"
  # _plot_map("sel_clrmap_tr_rgh_Lx_model1", \
  #     xy=xy, z=sel, ylabel=yl, xlabel=xl, zlabel=zl, showmarkers=False, tag=tag, line=np.array([[1.0, -1.6], [1.3, -1.2]]))

    # rgh also shows asymptotic behavior

    ################################
    #### (4) plot effective relations: U_SHE/dG_red ; dG_des/dG_ads ; n * D/Lx ; rgh
    # sampling
 #  datafile = "model_data_effective.pkl"
 #  dat = sample_data(datafile, rdes=[0.2], rads=np.arange(0.0, 0.4, 0.05), \
 #          rred=[1.1], Us=np.arange(-1.5,-1.15,0.05), Ds=[20e-10], Lxs=np.logspace(-3,-2,6), rghs=[1,5,10,25], mdls=[1]) #,25,50], mdls=[1])

 #  # i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1
 #  d_fixed = {2:0.3, 4:-1.3, 6:0.001, 7:1}
 #  # (a) U_SHE vs. roughness
 #  filename = "sel_clrmap_effmodel_USHE_rgh"; ixy = [4,7]; fixed = [2, 6]
 #  yl = r"$\rho$"; xl = r"U$_{\mathrm{SHE}}$ (V)"
 #  plot_effective_model(filename, dat, ixy, fixed, ylabel=yl, xlabel=xl)
 #  # (b) U_SHE vs. Dx/Lx
 #  filename = "sel_clrmap_effmodel_USHE_DxLx"; ixy = [4,6]; fixed = [2, 7]
 #  yl = r"$D_x / L_x$ (10$^5 \times$m s$^{-1}$)"; xl = r"U$_{\mathrm{SHE}}$ (V)"
 #  plot_effective_model(filename, dat, ixy, fixed, ylabel=yl, xlabel=xl)
 #  # (c) U_SHE vs. dG_equ
 #  filename = "sel_clrmap_effmodel_USHE_dG"; ixy = [4,2]; fixed = [6, 7]
 #  yl = r"$\Delta G_{\mathrm{ads-des}}$ (eV)"; xl = r"U$_{\mathrm{SHE}}$ (V)"
 #  plot_effective_model(filename, dat, ixy, fixed, ylabel=yl, xlabel=xl)
 #  # (d) rgh vs. Dx/Lx
 #  filename = "sel_clrmap_effmodel_rgh_DxLx"; ixy = [7,6]; fixed = [2, 4]
 #  yl = r"$D_x / L_x$ (10$^5 \times$m s$^{-1}$)"; xl = r"$\rho$"
 #  plot_effective_model(filename, dat, ixy, fixed, ylabel=yl, xlabel=xl)
 #  # (e) rgh vs. dG_equ
 #  filename = "sel_clrmap_effmodel_rgh_dG"; ixy = [7,2]; fixed = [6, 4]
 #  yl = r"$\Delta G_{\mathrm{ads-des}}$ (eV)"; xl = r"$\rho$"
 #  plot_effective_model(filename, dat, ixy, fixed, ylabel=yl, xlabel=xl)

 


 #  args = (1, 0.2, 0.3, 1.1, -1.4, Dx, Lx, roughness)
 #  output = run_model(*args)

 #  print(output)
 #  pr = output['prod']
 #  print(pr[0]/pr.sum())

#       sel = pr[0] / pr.sum()
#       #print(U_SHE, eng_red, sel)
#       #print(U_SHE, sel)
#       print( pr, pr.sum(), p0_norm) # pr is in TOF

    # mol/cm2
    # TOF --> mA/cm2
 #  print( pr, pr.sum(), p0_norm) # pr is in TOF
 #  pr_f = convert_TOF2flux(deepcopy(pr), roughness)
 #  print(pr_f, pr_f.sum())
 #  pr_j = convert_TOF2current(deepcopy(pr), roughness, nel_1=2, nel_2=4)
 #  print(pr_j, pr_j.sum())
    


 #  c_D_kt = 14.0224e-10    # m2/s +/- 3.16e-10 (10.1021/ie201944h0)
 #  c_D_oh = 52.7e-10       # m2/s              (10.1073/pnas.1713164114)
 #  c_k_sol = 5.26e4 * 1e-3 # m3/(mol*s)        (10.1139/v00-032)

 #  c_H_co = 1100           # Henry's constant for CO   (10.5194/acp-15-4399-2015)
 #  c_D_co = 20.3e-10       # m2/s              (Cussler, E. L. (1997). Diffusion: Mass Transfer in Fluid Systems)


    
  # ################################################
  # #### sampling pot and rghs and pH - Fig. 4a ####
  # dat = sample_polarization_curves()
  # plot_polarization_sel('Fig4a_selectivity_pol', dat)
  # ################################################
   

 #  # estimate d_levic --> for comparison of tested diffusion layer thicknesses
 #  Ds = [15e-10, 20e-10]
 #  #Ds = [1e-10, 60e-10]
 #  #D_ketene = 14.0224e-6 * 1e-4 # cm2/s --> m2/s
 #  visc =  (0.9980 +  1.1200)/2.0 * 1e-3 # mPa s --> Nm2/s for KOH 10.1021/je000019h
 #  omega = [261.8, 157.1, 2.6] # ~ 1500, 25 rpm
 #  dlevic = [1.6126 * D**(1./3.) * visc**(1./6.) * o**(-1./2.) *1e2 \
 #      for o in omega for D in Ds] # in cm
 #  print(dlevic)
 #  print('estimate of Levich diffusion layer: 8&60e-10 cm2/s = 2500 rpm: %.2e %.2e 1500 rpm: %.2e %.2e 25 rpm %.2e %.2e (cm)'%tuple(dlevic))

