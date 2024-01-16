""" 
    Module containing details of generic model

"""

import numpy as np
from model_Sel_CORAc.mkm.ODE_mkm_solver import ODE_mkm_solver
from model_Sel_CORAc.tools import load_pickle_file, write_pickle_file
from model_Sel_CORAc.data.energies_mkm import prep_surf_energies


#########################################################################
### functions describing ketene pathway reaction necessary for ODEint ###
#########################################################################

def dRI(c,t,*args):
    """
      function containing the ODE's to solve (dy/dx = ...) of mkm
      takes c (concentrations), t (time), *args (rate constants)
      [compare https://sam-dolan.staff.shef.ac.uk/mas212/notebooks/ODE_Example.html]
      
      Parameters
      ----------
      c : array (1D)
        vector containing the centrations [*, A, B, C, D, E]
      t : array (1D)
        array of timesteps to evaluate (necessary for ODEint)
      args : array (1D)
        array containing reaction rates and gas-pressures following
        reaction system (see explanation below)
          
      Returns
      -------
      docc : array 
        change in surface concentration of [*, A, B, C, D, E]

      The reaction system is hard-coded for the ketene pathway: 
      args[:16] = reaction constants for surface reactions
        g --> A  [CO adsorption] args[0:2]
        A --> B  [RDS effective] args[2:4]
        B --> Cg [reaction to ketene] args[4:6]
        C --> Cg [de-sorption ketene] args[6:8]
        C --> E  [protonation of adsorbed ketene] args[8:10]
        B --> D  [alternative route non-ketene] args[10:12] ***
        D --> E  [alternative route non-ketene] args[12:14] ***
        E --> Eg [desorption C2 product] args[14:16]
      *** (blocked by large barrier)
      args[16] = surface concentration Cg (ketene)
      args[17] = surface concentration Ag (CO)

      [*,A,B,C,D,E] = coverage vector
      dzero = d[*]/dt
      docc:
        dA/dt = +/- (Ag + * <-> A*) -/+ (A* <-> B*)
        dB/dt = +/- (A* <-> B*) -/+ (B* <-> Cg + *) -/+ (B* <-> D*)
        dC/dt = +/- (C_g + * <-> C*) -/+ (C* <-> E*)
        dD/dt = +/- (B* <-> D*) -/+ (D* <-> E*) (blocked by large barrier) --> taken out!!!
        dE/dt = +/- (C* <-> E*) +/- (D* <-> E*) -/+ (E* <-> Eg)


      In our multiscale approach the mkm is solved iteratively with transport
      NOTE: as an alternative to renormalization of total flux in the iterative 
        procedure, only the flux between C* <-> Cg can be normalized like:
        >> corr = 0
        >> if args[18] != 0: # normalization when a flux is set 
        >>    corr = args[18] - (c[2]*args[4] + c[3]*args[6] - args[16]*args[5] - args[16]*args[7])
        (...)
        >> _f2(np.array([args[16]] + c[[3,5]].tolist()), *args[6:8][::-1], *args[8:10]) - corr, # Cg, C, E --> towards C (subtract corr from C)
        which leads to a constraint solution of the mkm with very similar results as normalization of total flux. This, however,
        appears physically less correct


    """
    docc = [_f2(np.array([c[0]*args[13]]+c[1:3].tolist()), *args[0:4]), # A -> B
            _f2(np.array(c[[1,2]].tolist()+[args[12]]), *args[2:6]), # A B Cg
            _f2(np.array([c[0]*args[12]] + c[[3,4]].tolist()), *args[6:8][::-1], *args[8:10]), # Cg, C, E
            _f2(c[[3,4]].tolist()+[0.0], *args[8:10], *args[10:12]),
            ]
    dzero = -1 * sum(docc)
    return [dzero, *docc]


def dRI_2(c,t,*args):
    docc = [_f2(np.array([c[0]*args[13]]+c[1:3].tolist()), *args[0:4]), # A -> B
            _f2(c[[1,2,3]], *args[2:6]), # A B Cg
            _f3(np.array(c[[2,3]].tolist() + [c[0]*args[12] ,c[4]]), *args[4:6], *args[6:8], *args[8:10]), # B, C, Cg, E
            _f2(c[[3,4]].tolist()+[0.0], *args[8:10], *args[10:12]),
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
      based on DFT data its basic free energy corrections as obtained 
      from the data-module
      
      Parameters
      ----------
      surf : int
        surface miller index (100, 111, 110, 211)
      datkey : str
        key for which data to load
          
      Returns
      -------
      mkm : mkm.ODE_mkm_solver object
        mkm object with all pre-loaded data according to 
        data type and Cu-surface

    """
    ##  import DFT energies from data module:
    ##  --> free energies (with corrections as obtained from DFT
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
    
    #########################################################
    #### add energy corrections according to manuscript #####
    ## block barrier for step 5 (C -> D -> E)
    # engs0[5][1] += 2.0
    
    ## replace RDS in here as well (replace desorption barriers etc outside)
    engs0[1][1] = engs0[1][0] + 1.45
    #########################################################

    # re-referencing for first step to 0 (does not affect result)
    # better readability for debugging
    for i in range(len(engs0)):
        engs0[i][1] -= engs0[i][0]
        engs0[i][2] -= engs0[i][0]
        engs0[i][0] = 0.0
    return(_simple_model_std(engs0, k_dRI))


def _simple_model_std(engs0, k_dRI):
    """
      helper-function to prepare ketene pathway model based on engs0
      hard coded --> echem, nint, p_ind, r_ind (i_diss into transport solver)
    """
    echem = [False] + [True] * 2 + [False] + [True] * 2 # change for 'chemical' steps
    nint = 4 # number of intermediates
    # Desorption needs to be first rate to evaluate flux
    if k_dRI == 1:
        p_ind = [2,3,4] # intermediate(s) leading to a product [*, A, B, C, E]
        r_ind = [4,6,10] # rate leading to product
        i_diss = 2 # needs to be plugged into transport solver
        dRI_m = dRI
    elif k_dRI == 2:
        p_ind = [3,4] # intermediate(s) leading to a product [*, A, B, C, E]
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


