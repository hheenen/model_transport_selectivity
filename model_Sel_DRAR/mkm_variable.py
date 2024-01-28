""" 
    Module containing details of generic model to simulate the 
    desorption--re-adsorption--reaction (DRAR) rate-determining step
    (see manuscript)

"""

import numpy as np
from model_Sel_CORAc.mkm.ODE_mkm_solver import ODE_mkm_solver
from model_Sel_CORAc.tools import load_pickle_file, write_pickle_file
from model_Sel_CORAc.data.energies_mkm import prep_surf_energies


##########################################################
### functions describing two models of DRAR for ODEint ###
##########################################################

##############################################################################
### description of dRI and dRI_2                                           ###
###                                                                        ###
### function dRI and dRI_2 contain the ODE's to solve (dy/dx = ...) of mkm ###
### takes c (concentrations), t (time), *args (rate constants)             ###
### [compare e.g.                                                          ###
###  https://sam-dolan.staff.shef.ac.uk/mas212/notebooks/ODE_Example.html] ###
###                                                                        ###
###   The reaction system is hard-coded for the desorption--re-adsorption- ###
###   -reaction model (see manuscript)                                     ###
###   args[:16] = reaction constants for surface reactions                 ###
###     Ag --> A*  [adsorption A] args[0:2]                                ###
###     A* --> B*  [RDS] args[2:4]                                         ###
###     B* --> Cg  [reaction to desorbed C] args[4:6] - dRI(_1)            ###
###     B* --> C*  [reduction to C] args[4:6] - dRI_2                      ###
###     C* --> Cg  [de-sorption C] args[6:8]                               ###
###     C* --> D*  [reduction of C] args[8:10]                             ###
###     D* --> Dg  [desorption of product] args[10:12]                     ###
###   args[12] = surface activity Cg                                       ###
###   args[13] = surface activity Ag                                       ### 
##############################################################################


def dRI(c,t,*args):
    """
      Parameters
      ----------
      c : array (1D)
        vector containing the centrations [*, A, B, C, D, E]
      t : array (1D)
        array of timesteps to evaluate (necessary for ODEint)
      args : array (1D)
        array containing reaction rates and gas-pressures (see above)
          
      Returns
      -------
      docc : array 
        change in surface concentration of [*, A, B, C, D, E]

      > function for ODE describes differential equations for coverages
      > A-D and *=1-SUM(A-D) according to model 1
      [*,A,B,C,D] = coverage vector
      dzero = d[*]/dt
      docc:
        dA/dt = +/- (Ag + * <-> A*) -/+ (A* <-> B*)
        dB/dt = +/- (A* <-> B*) -/+ (B* <-> Cg + *)
        dC/dt = +/- (C_g + * <-> C*) -/+ (C* <-> E*)
        dD/dt = +/- (C* <-> D*) -/+ (D* <-> Dg)

    """
    docc = [_f2(np.array([c[0]*args[13]]+c[1:3].tolist()), *args[0:4]), # Ag, A, B
            _f2(np.array(c[[1,2]].tolist()+[args[12]]), *args[2:6]), # A, B, Cg
            _f2(np.array([c[0]*args[12]] + c[[3,4]].tolist()), *args[6:8][::-1], *args[8:10]), # Cg, C, D
            _f2(c[[3,4]].tolist()+[0.0], *args[8:10], *args[10:12]), # C, D, Dg
            ]
    dzero = -1 * sum(docc)
    return [dzero, *docc]


def dRI_2(c,t,*args):
    """
      Parameters
      ----------
      c : array (1D)
        vector containing the centrations [*, A, B, C, D, E]
      t : array (1D)
        array of timesteps to evaluate (necessary for ODEint)
      args : array (1D)
        array containing reaction rates and gas-pressures (see above)
          
      Returns
      -------
      docc : array 
        change in surface concentration of [*, A, B, C, D, E]

      > function for ODE describes differential equations for coverages
      > A-D and *=1-SUM(A-D) according to model 1
      [*,A,B,C,D] = coverage vector
      dzero = d[*]/dt
      docc:
        dA/dt = +/- (Ag + * <-> A*) -/+ (A* <-> B*)
        dB/dt = +/- (A* <-> B*) -/+ (B* <-> C*)
        dC/dt = +/- (B* <-> C*) +/- (C_g + * <-> C*) -/+ (C* <-> E*)
        dD/dt = +/- (C* <-> D*) -/+ (D* <-> Dg)

    """
    docc = [_f2(np.array([c[0]*args[13]]+c[1:3].tolist()), *args[0:4]),  # Ag, A, B
            _f2(c[[1,2,3]], *args[2:6]), # A, B, C
            _f3(np.array(c[[2,3]].tolist() + [c[0]*args[12] ,c[4]]), *args[4:6], *args[6:8], *args[8:10]), # B, C, Cg, D
            _f2(c[[3,4]].tolist()+[0.0], *args[8:10], *args[10:12]), # C, D, Dg
            ]
    dzero = -1 * sum(docc)
    return [dzero, *docc]
    

def _f2(c, kf_a, kb_a, kf_b, kb_b):
    """
      helper-function to solve compute dc[1]/dt for reaction 
      with two connecting pionts: c[0] <-k_a-> c[1] <-k_b-> c[2]

    """
    return kf_a * c[0] - kb_a * c[1] - kf_b * c[1] + kb_b * c[2]


def _f3(c, kf_a, kb_a, kf_b, kb_b, kf_c, kb_c):
    """
      helper-function to solve compute dc[1]/dt for reaction 
      with three connecting pionts: c[0] <-k_a-> c[1] <-k_b-> c[2] 
                                                 c[1] <-k_c-> c[3]

    """
    return kf_a * c[0] - kb_a * c[1] - kf_b * c[1] + kb_b * c[2] - kf_c * c[1] + kb_c * c[3]



#########################################################################
###         wrapper functions to create ready-to-use mkm              ###
#########################################################################

def get_mkm(surf, k_dRI=1, datkey=None):
    """
      function to prepare mkm model (see mkm.ODE_mkm_solver.py)
      based on input of free energies
      
      Parameters
      ----------
      surf : int
        surface miller index (100, 111, 110, 211)
      k_dRI : int
        1 or 2 for model 1 or model 2 according to dRI()/dRI_2()
      datkey : str
        key for which data to load
          
      Returns
      -------
      mkm : mkm.ODE_mkm_solver object
        mkm object with all pre-loaded data

    """
    ##  import DFT energies from data module (Cu data from prev. model)
    ##  --> to be overwritten by input data as listed in paper
    dat = prep_surf_energies()

    engs0 = []
    tags = dat[surf]['tags']
    engs = dat[surf]['engs']

    ## convert energy dict to energy list
    ## truncate mechanism in data.rct_data further to:
    ## A -> B       C -> E -> Eg  (also include C -> D -> E; but blocked below)
    ##      |-> Cg -^
    
    ## fixed order for short mechanism conversion of data.rct_data 
    eng_lst = [[0,0,1], [2,4,1], [5,5,1], [6,6,-1], [7,7,1], [8,8,1]]

    ## use eng_lst to make energy list
    for i in range(len(eng_lst)):
        # function to list IS-TS-FS ... 
        a = np.array([_return_order_engs(engs[l]) \
            for l in range(eng_lst[i][0], eng_lst[i][1]+1)]).flatten()
        if eng_lst[i][2] < 0:
            a = a[::-1]
        engs0.append(_reduce_ee(a))
    
    # set RDS - as indicated in manuscript
    engs0[1][1] = engs0[1][0] + 1.45 # general RDS

    # re-referencing for first step to 0 (does not affect result)
    # better readability for debugging
    for i in range(len(engs0)):
        engs0[i][1] -= engs0[i][0]
        engs0[i][2] -= engs0[i][0]
        engs0[i][0] = 0.0
    return(_simple_model_std(engs0, k_dRI))


def _simple_model_std(engs0, k_dRI):
    """
      helper-function to prepare DRAR model based on engs0
      hard coded --> echem, nint, p_ind, r_ind (i_diss into transport solver)
    """
    echem = [False] + [True] * 2 + [False] + [True] * 2 # change for 'chemical' steps
    nint = 4 # number of intermediates
    # Desorption needs to be first rate to evaluate flux
    if k_dRI == 1:
        p_ind = [2,3,4] # intermediate(s) leading to a product [*, A, B, C, D]
        r_ind = [4,6,10] # rate leading to product
        i_diss = 2 # needs to be plugged into transport solver
        dRI_m = dRI
    elif k_dRI == 2:
        p_ind = [3,4] # intermediate(s) leading to a product [*, A, B, C, D]
        r_ind = [6,10] # rate leading to product
        i_diss = 1 # needs to be plugged into transport solver
        dRI_m = dRI_2

    mkin = ODE_mkm_solver(dRI_m, echem, engs0, nint, p_ind, r_ind)
    return(mkin)


def _return_order_engs(ee):
    """
      helper-function to prepare energy dict
    """
    if len(ee) == 3:
        return([ee['IS'], ee['TS'], ee['FS']])
    else:
        return([ee['IS'], ee['FS']])


def _reduce_ee(el):
    """
      helper-function to prepare energy list
    """
    if len(el) > 2:
        return([el[0], max(el[1:-1]), el[-1]])
    else:
        return([el[0], max(el[0],el[-1]), el[-1]])


