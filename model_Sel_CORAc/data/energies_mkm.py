""" 
    Module to prepare energetics for microkinetic model on basis 
    of standard catmap energy input files

"""

from model_COR_Ac.data.eval_thermodynamics import _apply_rune_corrections, \
    _evaluate_all_vib_contributions, _relate_kW_CHE, get_rxn_state


def _prep_FE_CHE_from_efile(engdat, pressures, mechanism, dat_rxn, surf, pH, voltage, rune=False):
    """
      function to apply functions from `data.eval_thermodynamics` to evaluate 
      free energies of reactions
      
      Parameters
      ----------
      engdat : dict
        nested dictionary containing energies and vibrational frequencies
      pressures : dict
        dictionary with gas pressures used for free energy evaluation
      mechanism : dict
        dictionary with mechanism (list) for each facet
      dat_rxn : list
        list with strings describing reaction network (catmap format)
      surf : int
        integer for surface i.e. 100, 111, 110, 211
      pH : float
        pH to evaluate energies at
      voltage : float
        voltage (vs. SHE) to evaluate energies at
      rune : bool
        whether rune corrections should be applied
        
          
      Returns
      -------
      engs : array 
        array with relative reaction energetics
      tags : array
        array with according surface state

    """
    # evaluate free energy contributions according to entropy changes and pH
    e_dat = _prep_edat(engdat, rune)
    
    # hardcoded dictionaries (not necessary in this work)
    sigma_dat = {'sigma_params':{}, 'sigma_input':[1.0,0.0]}
    sites_dict = {'g': ['None', 'gas'], 't': ['Cu', surf], 'sl': ['Cu', 'dl']}
    
    # evaluate final free energies of reactions including pressures CHE, etc. 
    # pH and U_SHE = 0 so effectively only pressures and thermalization considered
    tags = []; engs = []
    for n_rxn in mechanism:
        rtags, rengs = get_rxn_state(n_rxn, dat_rxn, sites_dict, \
            e_dat, pH, voltage, sigma_dat, pressures, contributions=['fE','vc','pH_CHE','U_CHE','pc'])
        tags.append(rtags); engs.append(rengs)
    engs = _combine_rcts(engs)
    return(engs, tags)


def _combine_rcts(engs):
    """
      helper-function to relate all reaction energies to each other

    """
    estart = 0.0
    for i in range(0,len(engs)):
        engs[i] = {k:engs[i][k]-engs[i]['IS']+estart for k in engs[i]}
        estart = engs[i]['FS']
    return(engs)


def _prep_edat(e_dat, rune=False):
    """
      helper-function wrapping some pre-computations
      before full free energy evaluation in `get_rxn_state`

    """
    if rune:
        e_dat = _apply_rune_corrections(e_dat)

    # add vib contributions
    e_dat = _evaluate_all_vib_contributions(e_dat)
    
    # precalc ele, OH_g and H_g and add to edat
    e_dat = _relate_kW_CHE(e_dat) # is at pH=0
    return(e_dat)


def prep_surf_energies():
    """
      wrapper function to prepare dictionary compatible with
      mkm-class required energetics - energies returned correspond
      to unchanged DFT energies
      apart from typical free energy corrections, and rune corrections
          
      Returns
      -------
      fedat : dict
        array with energetics for reactions on 111, 100, 110, 211

    """
    from model_COR_Ac.data.rct_data import rxn, pressures, e_dat, mechs
    
    pH = 0; voltage = 0.0 # 0.5 RHE

    fedat = {s:{} for s in mechs}
    for s in mechs:
        engs, tags = _prep_FE_CHE_from_efile(e_dat, pressures, mechanism=mechs[s], dat_rxn=rxn, surf=str(s), pH=pH, voltage=voltage, rune=True)
        fedat[s].update({'engs': engs})
        fedat[s].update({'tags':tags})
    return(fedat)


