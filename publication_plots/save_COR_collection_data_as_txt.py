#!/usr/bin/env python

# imports experimental data
import sys
sys.path.append("../experimental_reference")
from read_data import *


def save_data(selCO, dall, seltag, rghdat):
    for k in selCO:
        ks = list(dall[k].keys()); 
        ks.remove('I(mA/cm2)'); ks.remove('V_RHE')
        ks = ['I(mA/cm2)'] + [kk for kk in ks if type(dall[k][kk]) != str]
        kkdat = [dall[k][kk][:,0] for kk in ks]
        kkdat = [selCO[k][:,0]] + kkdat + [selCO[k][:,1]] + [selCO[k][:,2]]
        kkdat = np.array(kkdat).T
        ks = ['U_SHE'] + ks + [seltag] + ['sig-sel']
        hdl1 = "data from %s with roughness %.1f, faradaic efficiencies given if not otherwise indicated\n"
        np.savetxt('../literature_data/dat_CO2RR_CO_%s.txt'%k.split(" ")[0], kkdat, \
            header=hdl1%(dall[k]['doi'], rghdat[k]) + "  ".join(ks), fmt='%.6e')


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
    save_data(selCO, dall, seltag, rghCO)


    # fig 3b data
    dat, d_all = load_Acdh_Cu_data()
    del d_all['CO2R@Cu-Sputtered $\\rho$=None']
    rghAc = {k:float(d_all[k]['roughness']) for k in d_all}

    # print-out data - make function for other examples
    seltag = 'sel Acdh / C2'
    save_data(dat, d_all, seltag, rghAc)


    # fig 4 data
    r_cupd, d_cupd, f_cupd = load_Ac_alloy_data("Ji_COR_CuPd", -1.4)
    f_cupd = {'-'.join(k.split('-')[1:]).split('+')[0]:f_cupd[k] for k in f_cupd}
    for k in f_cupd:
        f_cupd[k]['roughness'] = str(r_cupd[k][3] * 0.4)
    r_cuag, d_cuag, f_cuag= load_Ac_alloy_data("Sargent_COR_CuAg", -1.6, ulim=300)
    f_cuag= {'-'.join(k.split('-')[1:]).split('+')[0]:f_cuag[k] for k in f_cuag}
    for k in f_cuag:
        f_cuag[k]['roughness'] = str(r_cuag[k][3] * 5.0)
    # merge data
    dsel = d_cupd
    dsel.update(d_cuag)
    fdat = f_cupd
    fdat.update(f_cuag)
    rghdat = {}
    for k in r_cupd:
        rghdat.update({k:r_cupd[k][3]*0.4})
    for k in r_cuag:
        rghdat.update({k:r_cuag[k][3]*5.0})

    print({k:r_cupd[k][0] for k in r_cupd})
    print(rghdat)

    # print-out data - make function for other examples
    seltag = 'sel Ac / C2'
    save_data(dsel, fdat, seltag, rghdat)

