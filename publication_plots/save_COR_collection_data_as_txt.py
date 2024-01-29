#!/usr/bin/env python

# imports experimental data
sys.path.append("../experimental_reference")
from read_data import *


def save_data(selCO, dall, seltag):
    for k in selCO:
        ks = list(dall[k].keys()); 
        ks.remove('I(mA/cm2)'); ks.remove('V_RHE')
        ks = ['I(mA/cm2)'] + [kk for kk in ks if type(dall[k][kk]) != str]
        kkdat = [dall[k][kk][:,0] for kk in ks]
        kkdat = [selCO[k][:,0]] + kkdat + [selCO[k][:,1]]
        kkdat = np.array(kkdat).T
        ks = ['U_SHE'] + ks + [seltag]
        hdl1 = "data from %s with roughness %.1f, faradaic efficiencies given if not otherwise indicated\n"
        np.savetxt('../literature_data/dat_CO2RR_CO_%s.txt'%k, kkdat, \
            header=hdl1%(dall[k]['doi'], rghCO[k]) + "  ".join(ks), fmt='%.6e')


if __name__ == "__main__":

    # fig 3a data
    
    selCO, rghCO, phCO, dall = load_CO_Cu_data_collection()
    kout = ['Huang-110', 'Huang-100', 'Wang-NSph-Cu+pH14.0']
    selCO = {k:selCO[k] for k in selCO if k not in kout}
    rghCO['Huang-OD-Cu'] = 32 # Huang-OD-Cu is roughly roughness = 30 estimated from current; Huang_100 is also much higher in roughness...!
    rghCO['Kuhl-pc-Cu'] = 1.1
    rghCO['Huang-111'] = 1.2
    rghCO['Huang-110'] = 1.3 # roughness is actually 4 according to estimate
    rghCO['Huang-100'] = 1.4 # roughness is actually 4 according to estimate
    rghCO['Kanan-pc-Cu'] = 30.1 # from total current totally visible in paper that roughness = 30
    
    rghCO = {k:float(rghCO[k]) for k in rghCO}
    selCO = {k:selCO[k] for k in selCO if phCO[k] < 10.}
    #selCO = {k:selCO[k] for k in selCO if phCO[k] > 10.}
    rghCO = {k:rghCO[k] for k in selCO}
    phCO = {k:phCO[k] for k in selCO}

    # print-out data - make function for other examples
    
    seltag = 'sel CO / C1+C2'
    save_data(selCO, dall, seltag)


    # fig 3b data
