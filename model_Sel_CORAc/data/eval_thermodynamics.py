""" 
    Module to with functions to evaluate free energies from
    DFT data following standard surface science approaches
    -- note that some (finally) unused compabilities are included
    here: CHE, pH and field-effect evaluation. Maybe used to 
    extend the simple mkm

"""


import numpy as np
from copy import deepcopy
from ase import units
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from ase.build import molecule


######################################################################
####  Functions to correct and evaluate free energy contributions ####
######################################################################

def _apply_rune_corrections(mkm_engs):
    ''' 
      helper function to apply `Rune`-corrections for systematic DFT errors
      following 10.1039/c5cy01332a:
      correction of 0.15 eV per C=O and 0.1 eV per H2
      to be applied to catmap-like energy dict

    '''
    e_co = 0.15 #e_h2 = 0.666666 * e_co
    co_dict = {
     'CO2':2, 'H2':0.6666667, 'CH3COOH':1, 'CH3CH2OH':0, 'CH2CO':1, 'H2CCO':1, 'H2O':0, 'CH4':0, 'O2':0, 'ele':0, 'C2H4':0,
     'clean':0, 'OH':0, 'CO':1, 'H':0, 'CO-':1, 'CO2-':2, 'O':0,
     'COH':0, 'CHO':1, 'CH':0, 'CH2':0,'CH3':0,
     'H-CO':1,
     'CH2O':1, 'CHOH':0, 'CH2OH':0, 'C':0,
     'H-H':0, 'H2O-ele':0, 'H-ele':0, 'H-H2O-ele':0, 'H2-ele':0,
     'COOH':1, 'CO-OH':1, 'COO-H-ele':1, 'COOH-OH-ele':1, 'CO-OH-ele':1, 
     'COOH-OH':1, 'COOH-OH2-ele':1, 'H2O-CO-ele':1, 'H2O-CO-ele':1,
     'OCCO':2, 'OCCHO':2, 'OCCOH':1, 'OC-CO':2, 'OCCO-H2O-ele':2, 'OCC-OH-ele':1,
     'OCC':1, 'OCCH':1, 'OCCH2':1, 'HCCHO':1,
     'CC':0, 'CCH2':0, 'CCH2O':0, 'CCHO':1, 'CCOH':0, 'HCCH':0, 'HCCH2O':0, 'HCCHOH':0, 
     'HCCOH':0, 'HOCCH2O':0, 'HOCCHO':1, 'HOCCHOH':0, 'HOCO':1, 'OCCH2O':1, 'OCCH2OH':1,
     'OCCHOH':1, 'OHCCH2':1, 'H3CCHO':1, 
     'OCCH3':1, 'HOCCH2':0, 'HOCCOH':0, 'H2CCH2O':0, 'H2CCHO':1, 'H3CCH2O':0, 'H3CCOH':0,
     'CH3COO':1, 'OCCH2-':1,
     'OCC-H2O-ele':1, 'OCCH-H2O-ele':1, 'OCCH2-H2O-ele':1, 'HCOC-H2O-ele':1, 'OCCH2-OH2-ele':1, 'OCCH2-ele-H2O':1,
    }
    mkm_engs = deepcopy(mkm_engs)
    for surf in mkm_engs:
        for site in mkm_engs[surf]:
            for ads in mkm_engs[surf][site]:
                eng = float(mkm_engs[surf][site][ads]['fE']) + e_co*co_dict[ads]
                mkm_engs[surf][site][ads]['fE'] = eng 
    return(mkm_engs)


def _evaluate_all_vib_contributions(e_dat):
    ''' 
      helper function to evaluate vibrational free energy contributions

    '''
    out_dat = {surf:{site:{} for site in e_dat[surf]} for surf in e_dat}
    # compute surface vibs
    for surf in e_dat:
        if surf != 'None':
            for site in e_dat[surf]:
                for ads in e_dat[surf][site]:
                    dads = e_dat[surf][site][ads]
                    vc = 0.0
                    if len(dads['vib']) > 2:
                        vc = get_dG_from_vib(dads['vib'])
                    out_dat[surf][site].update({ads:{'fE':dads['fE'],'vc':vc}})
    for ads in e_dat['None']['gas']:
        dads = e_dat['None']['gas'][ads]
        vc = 0.0
        if len(dads['vib']) > 2 and ads not in ['OH']: # are replaced later
            vc = get_gas_dG_from_vib(dads['vib'], ads)
        out_dat['None']['gas'].update({ads:{'fE':dads['fE'],'vc':vc}})
    return(out_dat)


def get_dG_from_vib(vibfreq):
    ''' 
      helper function to evaluate free energy for harmonic adsorbate
      using ASE's thermopackage

    '''
    vibfreq = _vibstring_to_list(vibfreq)
    gibbs = HarmonicThermo(vib_energies = np.array(vibfreq)*units.invcm)
    eng_corr = gibbs.get_helmholtz_energy(300, verbose=False)
    return(eng_corr)


def get_gas_dG_from_vib(vibfreq, ads):
    ''' 
      helper function to evaluate free energy for ideal gas 
      using ASE's thermopackage

    '''
    vibfreq = _vibstring_to_list(vibfreq)
    vibenergies = np.array(np.zeros(6).tolist() + vibfreq)*units.invcm
    
    # hardcoded parameters
    # linear / symmetry number
    gas_prop = {'H2':['linear',2], 'CO':['linear',1], 'CH3COOH':['nonlinear',2], 'CH4':['nonlinear',12],\
        'H2O':['nonlinear',2], 'CO2':['linear',2], 'CH2CO':['nonlinear',2], 'OCCH2':['nonlinear',2], \
        'CH3CH2OH':['nonlinear',1], 'C2H4':['nonlinear',6]}
    # get atoms
    if ads not in ['OCCH2', 'CH2CO']:
        atoms = molecule(ads)
    else:
        atoms = molecule('H2CCO')

    gibbs = IdealGasThermo(vib_energies=vibenergies, geometry=gas_prop[ads][0], \
                        spin=0, symmetrynumber=gas_prop[ads][1], atoms=atoms)
    eng_corr = gibbs.get_gibbs_energy(300, pressure=101325, verbose=False) # 1.0 atm
    return(eng_corr)


def _vibstring_to_list(vibstring):
    ''' 
      helper function to convert saved vib-string to list

    '''
    if vibstring[-2:] == ',]':
        vibstring = vibstring[:-2]+']'
    slist = [float(v) for v in vibstring[1:-1].split(',')]
    return(slist)
    

def _relate_kW_CHE(e_dat):
    ''' 
      helper function to compute water equilibrium to relate pH etc.

    '''
    # precalc ele, OH_g and H_g and add to edat is at pH=0
    gas_dat = e_dat['None']['gas']
    gas_dat.update({'ele':{'fE':0.0, 'vc':0.0}})
    gas_dat.update({'H':{k:gas_dat['H2'][k]/2. for k in gas_dat['H2']}})
    gas_dat.update({'OH':{k:gas_dat['H2O'][k]-gas_dat['H'][k] for k in gas_dat['H2O']}}) #add comprehension like above
    e_dat['None']['gas'].update(gas_dat)
    return(e_dat)
    

######################################################################
####  Functions to handle and interprete catmap-like rcts strings ####
######################################################################
    
def get_rxn_state(n_rxn, dat_rxn, sites_dict, e_dat, pH, U, sigma_dat, \
                  pressures, contributions):
    # decompose rxn string
    drct = _decompose_reaction(dat_rxn[n_rxn])
    
    eng = {}
    # get IS energies
    d_IS = get_state(e_dat, sites_dict, \
        _filter_rct(drct['IS'], ['2*_t','*_t','*_dl','*_sl','*_sdl'], exclude=True),
        pH, U, sigma_dat, pressures)
    eng.update({'IS':sum([d_IS[c] for c in contributions])})
    
    # get FS energies
    d_FS = get_state(e_dat, sites_dict, \
        _filter_rct(drct['FS'], ['2*_t','*_t','*_dl','*_sl','*_sdl'], exclude=True),
        pH, U, sigma_dat, pressures)
    eng.update({'FS':sum([d_FS[c] for c in contributions])})
    
    # get TS energies
    if 'TS' in drct:
        d_TS = get_state(e_dat, sites_dict, \
            _filter_rct(drct['TS'], ['2*_t','*_t','*_dl','*_sl','*_sdl'], exclude=True),
        pH, U, sigma_dat, pressures)
        eng.update({'TS':sum([d_TS[c] for c in contributions])})

    return({st:'+'.join(_filter_rct(drct[st], \
            ['2*_t','*_t','*_dl','ele_g','OH_g','H_g','H2O_g'])) for st in drct},
            {st:eng[st] for st in drct})

        
def get_state(e_dat, sites_dict, rxn_str, pH, U, sigma_dat, pressures):
    """
      separately calculate all various contribution to a state
      fE = formation (electronic energy)
      vc = vibration contribution 
      pc = pressure contribution 
      pH_CHE = pH effect in electrochemistry
      U_CHE = potential contribution
      qc = (surface) charge contribution (field-effect)
    """
    kB = 8.617e-5
    out = {k:0.0 for k in ['fE','vc','pc','pH_CHE','U_CHE','qc']}
    for im in rxn_str:
        # decompose intermedates
        spec, tag = im.split('_')
        spec = spec.split('*')[0]
        # for multiplication with stoichiometry
        nstoich, spec = _str_int(spec)
        
        sident = sites_dict[tag]
        
        # get electronic energy
        out['fE'] += nstoich * e_dat[sident[0]][sident[1]][spec]['fE']

        # get vibration contriubution
        out['vc'] += nstoich * e_dat[sident[0]][sident[1]][spec]['vc']
        
        # get pressure contribution
        if spec+'_g' in pressures and sident[1] == 'gas':
            out['pc'] += nstoich * 300*kB*np.log(pressures[spec+'_g'])
        
        # get CHE contribution
        if spec == 'OH' and sident[1] == 'gas':
            out['pH_CHE'] += nstoich * .0592*pH
        elif spec == 'H' and sident[1] == 'gas':
            out['pH_CHE'] -= nstoich * .0592*pH
        if spec == 'ele' and sident[1] == 'gas':
            out['U_CHE'] -= nstoich * U
        if im.find('ele') != -1 and spec != 'ele':
            out['U_CHE'] -= nstoich * 0.5 * U # is this right?
        
        # get surface charge / field contribution
        if im in sigma_dat['sigma_params']:
            sigma=np.poly1d(sigma_dat['sigma_input'])(U)
            sparams = sigma_dat['sigma_params'][im]
            p = np.poly1d(sparams)
            out['qc'] += p(sigma) - sparams[-1]
    return(out)


def _decompose_reaction(rxn):
    '''
      decompose the reaction string
    
    '''
    dat_rct = {}
    rxn_seq = [rx.strip() for rx in rxn.split(';')[0].split('<->')]

    # append IS, FS, TS
    dat_rct.update({"IS":[im.strip() for im in rxn_seq[0].split('+')]})
    dat_rct.update({"FS":[im.strip() for im in rxn_seq[-1].split('+')]})
    if len(rxn_seq) == 3:
        dat_rct.update({"TS":[im.strip() for im in rxn_seq[1].split('+')]})

    # throw out accidental '' following a '+' too much
    dat_rct = {st:[im for im in dat_rct[st] if im != ''] for st in dat_rct}
    
    # read beta
    beta = 0.5 #standard assumed
    if len(rxn.split(';')) > 1:
        beta = float(rxn.split(';')[1].split('=')[1])
    return(dat_rct)


def _filter_rct(rct_list, okeys=[], exclude=True):
    if exclude:
        return([im for im in rct_list if _str_int(im)[1] not in okeys])
    if not exlude:
        return([im for im in rct_list if _str_int(im) in okeys])

        
def _str_int(string):
    """ 
      helper function to remove leading int in str
    
    """
    try:
        nstoich = int(string[0])
        nstring = string[1:]
        return(nstoich, nstring)
    except ValueError:
        return(1, string)

