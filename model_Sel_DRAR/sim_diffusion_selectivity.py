#!/usr/bin/env python

"""

    This script mainly includes high-level function to run 
    the mkm-diffusion coupled selectivity model of the 
    desorption--re-adsorption--reaction (DRAR)

"""

# general imports
import os
import numpy as np
import pickle

# package imports
from model_Sel_DRAR.transport_iterative import solve_MS_X_analytical
from model_Sel_DRAR.mkm_variable import get_mkm



def sample_data(datafile, rdes, rads, rred, Us, Ds, Lxs, rghs, mdls):
    """
      function run model of DRAR reaction model based on 
      ranges of input parameters
      
      Parameters
      ----------
      datafile : string
        name of datafile to store data in
      rdes : iterable of floats
        values of desorption barries of C* to use
      rads : iterable of floats
        values of adsorption barries of C* to use
      rred : iterable of floats
        values of reduction barriers of C* to use
      Us : iterable of floats
        values U_SHE to use
      Ds : iterable of floats
        values diffusion coefficients to use [cm2/s]
      Lxs : iterable of floats
        values diffusion length to use [cm]
      rghs : iterable of floats
        values roughnesses to use
      rghs : iterable of ints
        [1, 2] for model 1 or 2 or both
          
      Returns
      -------
      out_data : dict 
        output data corresponding to input variables

    """

    # load pkl datafile (in hard-coded bk_pkl_files folder
    datafile = "bk_pkl_files/%s"%datafile
    sim_data = load_datafile(datafile)

    out_data = []
    n_count = 0

    # iterate through input variables
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


def run_model(i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness):
    """
      helper-function to DRAR reaction model for one set of parameters

    """
    pH = 14 # doesn't matter
    i_diss_m = {1:2, 2:1}
    mkin = make_mkm_simple_2(eng_des, eng_ads, eng_red, k_dRI=i_model)
    
    ts, ys, pr, p0_norm, p0 = solve_MS_X_analytical(mkin, U_SHE, pH, \
        Dx, Lx, i_diss=i_diss_m[i_model], i_rads=7, roughness=roughness)
    output = {'cov':ys[-1,:], 'prod':pr, 'conc':p0}
    return(output)


def make_mkm_simple_2(eng_des, eng_ads, eng_red, k_dRI):
    """
      helper-function to setup mkm with input energies

    """
    mkm = get_mkm(100, k_dRI) # load standard model to manipulate
    mkm.engs0[3] = [0.0, eng_des, eng_des-eng_ads]
    mkm.engs0[4][1] = eng_red
    return mkm


def get_data(data, *args):
    """
      helper-function to retrieve from data-dict sub-dict
      with model settings according to args:
      args = (i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness)

    """
    key = '-'.join([str(arg) for arg in args])
    if key in data:
        return(data[key])
    else:
        return(None)
            

def update_data(data, output, *args):
    """
      helper-function to save sub-dict to data-dict
      ky corresponds to model settings according to args

    """
    key = '-'.join([str(arg) for arg in args])
    data.update({key: output})


def load_datafile(datafile):
    """
      helper-function to load central data file 
      holding data-dict of model

    """
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


def load_pickle_file(filename):
    """
      helper-function to load pickle file

    """
    with open(filename, 'rb') as pickle_file:
        data = pickle.load(pickle_file)
    return(data)


def write_pickle_file(filename, data):
    """
      helper-function to write pickle file

    """
    with open(filename, 'wb') as output:
        pickle.dump(data, output)


if __name__ == "__main__":
    pass

