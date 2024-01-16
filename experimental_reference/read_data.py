#!/usr/bin/env python

"""

This script reads preprepared and digitized data

"""

import os, sys
import numpy as np
from copy import deepcopy
sys.path.append("/Users/heenen/Nextcloud/DTU_related/Publications/Acetate_Selectivity/co_acetate/experimental_data/COR_data_paper")
from plot_paper_data import _load_pickle_file
sys.path.append("/Users/heenen/Nextcloud/DTU_related/Publications/Acetate_Selectivity/co_acetate/experimental_data")
from get_c2_selectivities import compute_C2_selectivity
from scipy.optimize import curve_fit
import pandas as pd 

def load_CO_Cu_data():
    abs_path = "/Users/heenen/Nextcloud/Projects_FHI/kinetics_transport/model_perspective/experimental_reference" 
   #dat = np.loadtxt(abs_path+"/dat_CO_Cu.txt")
   #header = ["rgh", "CO", "C2H4", "EtOH", "PropAldh", "Prop", "CH4"]
    dat = np.loadtxt(abs_path+"/dat_CO_Cu_full.txt")
    header = ["rgh", "CO", "C2H4", "EtOH", "PropAldh", "Prop", "CH4", "H2", "HCOO"]
    nel = dict(CO=2, C2H4=10, EtOH=10, PropAldh=10, Prop=12, CH4=8)#, H2=2, HCOO=2) # not interesting
    nC = dict(CO=1, C2H4=2, EtOH=2, PropAldh=3, Prop=3, CH4=1)#, H2=0, HCOO=1)

    curr = dat[:,1:].sum(axis=1) * dat[:,0]
    irgh = curr/ curr.min()
  # print(curr)
  # print(curr / curr.min())

    # sort data
    ddat = {header[i]:dat[:,i] for i in range(len(header))}
    del ddat['H2']; del ddat['HCOO']
    # normalize data to CO molar flow
    for k in nel:
        ddat[k] /= nel[k]
        ddat[k] *= nC[k]

    # compute selectivity
    dtot = np.array([ddat[k] for k in nel]).T.sum(axis=1)
    sel = ddat['CO'] / dtot
    return np.array([ddat['rgh'], sel]).T, np.array([irgh, sel]).T
   

def load_CO_Cu_data_collection():
    ffull = ["Bertheussen_COR-pcCu", "Bertheussen_COR",  "Huang_CO2R",
        "Kanan_CO2R-ODCu", "Kanan_COR-ODCu", "Kuhl_CO2", "Wang_COR", "Wang_COR-Cuflower", "Raciti_COR",
        "Jouny_COR", "Luc_COR", "Ma_CO2R",
        "Gregorio_CO2R", "Sargent_CO2R-CuN-C", "Zuettel_CO2R", "Kanan_COR-GDE", "Ji_COR_CuPd", "Li_COR-GDE", "Wang_CO2R_GDE"]

    dfull = _load_data_paper(ffull)
    # get folders that are CO2R
    dat = {k:dfull[k] for k in dfull if dfull[k]['mode'] == 'CO2R'}
    # hand out-select bad ones
    # del dat['Kuhl-pc-Cu']; del dat['Kanan-pc-Cu']


    d_sel = {}; d_rgh = {}; d_ph = {}
    for k in dat:
        if 'CO' in dat[k]:
            # compute CO selectivity:
            sel_CO = compute_CO_selectivity(dat[k])
            U_SHE = dat[k]['V_RHE'][:,0] - (0.059*float(dat[k]['pH']))
            
            d_sel.update({k:np.array([U_SHE, sel_CO]).T})
            #d_sel.update({k:np.array([dat[k]['V_RHE'][:,0], sel_CO]).T})
            d_rgh.update({k:dat[k]['roughness']})
            d_ph.update({k:float(dat[k]['pH'])})

    return(d_sel, d_rgh, d_ph)
        
        
def compute_CO_selectivity(dat):
    # std. dicts required
    stdkeys = ['V_RHE', 'I(mA/cm2)', 'mode', 'pH', 'catalyst', 'roughness', 'cell', 'doi']
    nel = dict(CO=2, C2H4=10, EtOH=10, PropAldh=10, Prop=12, CH4=8, H2=2, \
        HCOO=2, Ac=6, H3CHO=8, C2H6=12, MeOH=6, HOH2CCHO=6, C2H4O2H2=8, \
        C3H5OH=10, C2H5CHO=10, Me2CO=10, MeCOCH2OH=8)
    nC = dict(CO=1, C2H4=2, EtOH=2, PropAldh=3, Prop=3, CH4=1, H2=0, HCOO=1, \
        Ac=2, H3CHO=2, C2H6=2, MeOH=1, HOH2CCHO=2, C2H4O2H2=2, C3H5OH=3, C2H5CHO=3, Me2CO=3, MeCOCH2OH=3)
    tlm = {'Methane':'CH4', 'Ethylene':'C2H4', 'CO':'CO', 'Hydrogen':'H2','Ethane':'C2H6', 'Methanol':'MeOH',\
        'Formate':'HCOO', 'Ethanol':'EtOH', 'Acetate':'Ac', 'n-propanol':'Prop', 'Ac':'Ac', 'Acetaldehyde':'H3CHO',\
        'Glycolaldehyde':'HOH2CCHO', 'Ethylene-Glycol':'C2H4O2H2', 'Allyl-Alcohol':'C3H5OH', 'Propionaldehyde':'C2H5CHO',\
        'Acetone':'Me2CO', 'Hydroxyacetone':'MeCOCH2OH'}
    
    # now get molar selectivity (compare other CO data) for CO vs. C1 and C2
    mkeys = [kk for kk in dat if kk not in stdkeys and kk not in ['Hydrogen', 'Formate']]
    # normalize to smth proportional to molar ratio
    # (enough to be proportional for selectivity)
    jp = {kk:dat[kk] / nel[tlm[kk]] * nC[tlm[kk]] for kk in mkeys}
    jp_all = (np.array([jp[kk][:,0] for kk in jp]).T).sum(axis=1)
    sel_CO = jp['CO'][:,0] / jp_all
    return sel_CO


def load_Acdh_Cu_data_silly():
    abs_path = "/Users/heenen/Nextcloud/Projects_FHI/kinetics_transport/model_perspective/experimental_reference" 
    dat = np.loadtxt(abs_path+"/dat_AcdH_Cu.txt")
    header = ["U_RHE", "AcdH", "Ac", "EtOH"]
    # sort data
    ddat = {header[i]:dat[:,i] for i in range(len(header))}
    # normalize FE to proportional sel
    nel = dict(EtOH=8, Ac=4, AcdH=6)
    for sp in nel:
        ddat[sp] /= nel[sp]
    # compute selectivity
    dtot = np.array([ddat[k] for k in nel]).T.sum(axis=1)
    sel = ddat['AcdH'] / dtot
    return np.array([ddat['U_RHE'], sel]).T      
    

def load_Ac_alloy_data(key, curr_lim, ulim=None):
    dat = _load_data_paper([key])
    #dat = {_rename_key(dat[k]):dat[k] for k in dat}
    
    ### need to feed dummy roughness values
    for k in dat:
        if dat[k]['roughness'] == 'None':
            dat[k]['roughness'] = 1.0 # dummy value

    Ack = 'Acetate'
    sdat = compute_C2_selectivity(dat, Sads=Ack)
    dout = {k:np.array([sdat[k]['U_SHE'][:,0], sdat[k]['SC2_%s'%Ack][:,0]]).T for k in sdat}
    
    rrgh = {k:float(dat[k]['roughness']) if dat[k]['roughness'] != 'None' else np.nan for k in dat}
    rgh_min = min(rrgh.values()); rrgh = {k:rrgh[k]/rgh_min for k in rrgh}

    curr_i = {k:np.where(dout[k][:,0] > curr_lim)[0][-1] for k in dout}
    #print({k:dout[k][curr_i[k],0] for k in curr_i})
    jc = {k:dat[k]['I(mA/cm2)'][curr_i[k],0] for k in dat}
    #jc = {k:dat[k]['I(mA/cm2)'][curr_i[k],0]*(1-dat[k]['Hydrogen'][curr_i[k],0]) for k in dat}
    jc_min = min(jc.values()); jc = {k:jc[k]/jc_min for k in jc}

    # NOTE: not used, can probably be deleted
    #curr_i2 = {k:np.where(dout[k][:,0] > -1.5)[0][-1] for k in dout}
    curr_i2 = {k:np.absolute(dout[k][:,0] + 1.5).argmin() for k in dout}
    jc2 = {k:dat[k]['I(mA/cm2)'][curr_i2[k],0] for k in dat}
    jc_min2 = min(jc2.values()); jc2 = {k:jc2[k]/jc_min2 for k in jc2}
    
    # make exp fit to current for better estimate - seems also robust for Pd alloys
    dj = {k:np.array([sdat[k]['U_SHE'][:,0], dat[k]['I(mA/cm2)'][:,0]]).T for k in dat}
    #dj = {k:np.array([get_USHE(dat[k])[:,0], compute_CO2RR_partial(dat[k])]).T for k in dat} ## --> doesn't help
    # cut current to lim
    if ulim is not None:
        dj_full = deepcopy(dj)
        ilim = {k:np.where(dj[k][:,1] <= ulim)[0] for k in dj}
        dj = {k:dj[k][ilim[k],:] for k in dj}

    icrgh, icrgh_std = _fit_average_exp_rgh(dj)
    # _show_exp_fit(dj_full, ulim)

    # note that disordered Pd alloys have a lot of inactive surface (capacitance)
    # which is CO2RR inactive: total currents don't fit roughness values at all
    # roughness even appears to be on the order of ordered CuPd; compare jc and jc2
    # transformation of d-CuPd to CuPd in CO atmosphere tells us that Cu are pulled to the surface
    # it appears that average total current estimate (-1.4 V / most negative) gives a good estimate corresponding to selectivity
    # at -1.4 V very comparable potential values: 'Ji-Cu-NP+pH14.0': -1.396, 'Ji-CuPd+pH14.0': -1.3860000000000001, 'Ji-d-CuPd+pH14.0': -1.366, 'Ji-Cu3.4Pd+pH14.0': -1.376, 'Ji-Cu0.3Pd+pH14.0': -1.396} >>> use jc
   #for k in dout:
   #    print(k)
   #    print("%.2f   %.2f   %.2f   %.2f"%(float(dat[k]['roughness']), rrgh[k], jc[k], jc2[k]), dout[k][:,1])
   #    #print(dout[k])
   #assert False

    # selectivity vs. potential
    #dout = {k:dout[k] for k in ['Ji-Cu-NP+pH14.0', 'Ji-CuPd+pH14.0']}
    dout = {'-'.join(k.split('-')[1:]).split('+')[0]:dout[k] for k in dout}
    # roughness values
    drgh = {k:[float(dat[k]['roughness']), rrgh[k], jc[k], icrgh[k]] for k in rrgh}
    drgh = {'-'.join(k.split('-')[1:]).split('+')[0]:drgh[k] for k in drgh}
    return(drgh, dout)


def _fit_average_exp_rgh(dj):
    # determine b - range
    cc = _fit_exp_full(dj)
    print('tafel', cc)

    # fit with different variants of b
    rr = []; ks = list(dj.keys())
    #for cb in [min(cc), np.mean(cc), max(cc)]:
    for cb in cc:
        ca = []
        y = lambda x, a: a + np.exp(x*cb)
        for k in ks:
            popt, pcov = curve_fit(y, dj[k][:,0], dj[k][:,1])
            ca.append(popt[0])
        ca = np.array(ca)
        # relative roughness
        rrgh = ca-ca.min()+1
        rr.append(rrgh)
    # take mean of relative roughness
    rr = np.array(rr)
    print('rgh', np.mean(rr, axis=0), np.std(rr, axis=0), np.std(rr, axis=0)/np.mean(rr, axis=0))
    #print('rgh', np.mean(rr, axis=0)*0.4, np.std(rr, axis=0)*0.4, np.std(rr, axis=0)/np.mean(rr, axis=0))
    #print('rgh', np.mean(rr, axis=0)*5, np.std(rr, axis=0)*5, np.std(rr, axis=0)/np.mean(rr, axis=0))
    rrgh = np.mean(rr, axis=0)
    rrgh = {ks[i]:rrgh[i] for i in range(len(ks))}
    rrgh_std = np.std(rr, axis=0)
    rrgh_std = {ks[i]:rrgh_std[i] for i in range(len(ks))}
    return rrgh, rrgh_std


def _fit_exp_full(dj):
    # determine b - max/min
    y = lambda x, a, b: a + np.exp(x*b)
    cc = []
    for k in dj:
        popt, pcov = curve_fit(y, dj[k][:,0], dj[k][:,1], p0=[0.0, -3.0])
        cc.append(popt[1])
    return cc


def _show_exp_fit(dj, ulim, llim):

    ilim = {k:np.arange(dj[k][:,0].size).tolist() for k in dj}
    if ulim is not None:
        ilim = {k:np.where(dj[k][:,1] <= ulim)[0] for k in dj}
    if llim is not None:
        ilim = {k:ilim[k][np.where(dj[k][ilim[k],1] >= llim)[0]] for k in dj}
    dj_cut = {k:dj[k][ilim[k],:] for k in dj}
    cb = np.mean(_fit_exp_full(dj_cut))
    y = lambda x, a: a + np.exp(x*cb)
    
    import matplotlib.pyplot as plt
    # TODO: make use of xy_plot and make proper plot

    clr = {'Sarge-CuAg_0%':'r', 'Sarge-CuAg_0.25%':'y', 'Sarge-CuAg_0.5%':'g', 'Sarge-CuAg_1%':'b'}
    clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']; ks = list(dj_cut.keys())
    clr = {ks[i]:clrs[i%len(clrs)] for i in range(len(ks))}
    
    fig, ax = plt.subplots(2,1)
    for k in dj:
        popt, pcov = curve_fit(y, dj_cut[k][:,0], dj_cut[k][:,1])
        yy = [y(x, *popt) for x in dj[k][:,0]]
        ax[0].plot(dj[k][:,0], dj[k][:,1], color=clr[k])
        ax[0].plot(dj[k][:,0], yy, color=clr[k], ls=':')
        ax[1].plot(dj[k][:,0], dj[k][:,1], color=clr[k])
        ax[1].plot(dj[k][:,0], yy, color=clr[k], ls=':')

    ax[1].set_yscale('log')
    plt.show()


def load_Acdh_Cu_data():
    ffull = ["Bertheussen_COR-pcCu", "Bertheussen_COR",  "Huang_CO2R",
        "Kanan_CO2R-ODCu", "Kanan_COR-ODCu", "Kuhl_CO2", "Wang_COR", "Wang_COR-Cuflower", "Raciti_COR",
        "Jouny_COR", "Luc_COR", "Ma_CO2R",
        "Gregorio_CO2R", "Sargent_CO2R-CuN-C", "Zuettel_CO2R", "Kanan_COR-GDE",
        ]
    dfull = _load_data_paper(ffull)
    # get folders that have Acetaldehyde
    Ack = 'Acetaldehyde'
    dat = {k:dfull[k] for k in dfull if Ack in dfull[k]}
    # hand out-select bad ones
    del dat['Kuhl-pc-Cu']; del dat['Kanan-pc-Cu']

    print(dat.keys())
    
    dat = {_rename_key(dat[k]):dat[k] for k in dat}
    sdat = compute_C2_selectivity(dat, Sads=Ack)
    dout = {k:np.array([sdat[k]['U_SHE'][:,0], sdat[k]['SC2_%s'%Ack][:,0]]).T for k in sdat}
    dout = {k:_remove_zerox(dout[k]) for k in dout}
    return(dout)

def _remove_zerox(dat):
    indz = np.where(dat[:,0] != 0.0)[0]
    return(dat[indz,:])

def _rename_key(dat):
    key = r"%s@%s $\rho$=%s"%(dat['mode'], dat['catalyst'], dat['roughness'])
    #key = r"%s@%s $\rho$=%s $^{[]}$"%(dat['mode'], dat['catalyst'], dat['roughness'])
    return(key)


def _load_data_paper(folders):
    bpath = "/Users/heenen/Nextcloud/DTU_related/Publications/Acetate_Selectivity/co_acetate/experimental_data"
    data = {}
    for f in folders:
        dat = _load_pickle_file('%s/COR_data_paper/%s/%s.pkl'%(bpath,f,f))
        
        tag = f.split('_')[0][:5]
        for k in dat:
            if k != 'Cu_oh': # faulty potential
                data.update({tag+'-'+k:dat[k]})
    return(data)
    

def load_H2CO_Pt_data():
    abs_path = "/Users/heenen/Nextcloud/Projects_FHI/kinetics_transport/model_transport_selectivity/experimental_reference" 
    dat = np.loadtxt(abs_path+"/dat_H2CO_Pt.txt")
    dat2 = np.loadtxt(abs_path+"/dat_H2CO_Pt_2.txt")
    header = ["Pt loading", "H2CO", "HCOOH", "CO2"]
    # sort data
    ddat = {header[i]:dat[:,i] for i in range(len(header))}
    # no normalization necessary
    return dat[:,:2], dat2[:,:2]
    

def load_ORR_Pt_data():
    abs_path = "/Users/heenen/Nextcloud/Projects_FHI/kinetics_transport/model_perspective/experimental_reference" 
    dat1 = np.loadtxt(abs_path+"/dat_ORR_Pt_Pc.txt")
    dat2 = np.loadtxt(abs_path+"/dat_ORR_Pt_22.txt")
    dat3 = np.loadtxt(abs_path+"/dat_ORR_Pt_10.txt")
    rgh = [1., 0.22, 0.1]

    ind2 = []; ind3 = []
    for v in dat1[:,0]:
        ind2.append(np.absolute(dat2[:,0] - v).argmin())
        ind3.append(np.absolute(dat3[:,0] - v).argmin())
    dat2 = dat2[ind2,:]; dat3 = dat3[ind3,:]
    v = np.array([dat1[:,0], dat2[:,0], dat3[:,0]]).T
    v = np.around(np.mean(np.array([dat1[:,0], dat2[:,0], dat3[:,0]]).T, axis=1),4)

    out = {}
    # out = np.array([v, dat1[:,1], dat2[:,1], dat3[:,1]]).T
    # print(out[ind,:])
    #ind = [2, 9, 13]
    ind = [8, 9, 10, 11]
    #ind = range(dat1[:,0].size)
    for ii in ind:
        out.update({v[ii]:np.array([rgh, [dat1[ii,1], dat2[ii,1], dat3[ii,1]]]).T})
    return out
    

def give_rgh_estimates():
    ''' 
        hard-coded to yield roughnesses based on current estimates,
        includes alloy data, Huang-SC data, decaying data

    '''
  ##read_degr_data()
  ##assert False
    inputs = {"Ji_COR_CuPd":[None,None], "Huang_CO2R":[2.0,1000], "Sargent_COR_CuAg":[None,300]}
    out = {}
    for dk in inputs:
        print(dk)
        dat = _load_data_paper([dk])
        llim, ulim = inputs[dk]
    
        # make exp fit to current for better estimate
        dj = {k:np.array([get_USHE(dat[k])[:,0], dat[k]['I(mA/cm2)'][:,0]]).T for k in dat}
        #dj = {k:np.array([get_USHE(dat[k])[:,0], compute_CO2RR_partial(dat[k])]).T for k in dat} ## --> doesn't help
    
        # cut current to lim
        ilim = {k:np.arange(dj[k][:,0].size).tolist() for k in dj}
        if ulim is not None:
            ilim = {k:np.where(dj[k][:,1] <= ulim)[0] for k in dj}
        if llim is not None:
            ilim = {k:ilim[k][np.where(dj[k][ilim[k],1] >= llim)[0]] for k in dj}
        dj = {k:dj[k][ilim[k],:] for k in dj}
 
        icrgh, icrgh_std = _fit_average_exp_rgh(dj)
        out.update({dk:[icrgh, icrgh_std]})
        if dk == "Huang_CO2R":
            _show_exp_fit(dj, ulim=ulim, llim=llim)
    print(40*"#")
    for dk in out:
        print(dk, ' || '.join(["'%s' = %.2f +/- %.2f"%(k, out[dk][0][k], out[dk][1][k]) for k in out[dk][0]]))
    print(40*"#")
    print("Ji_COR_CuPd", [393.0, 109.5, 68.3, 188.4, 3.0])
    print("Sargent_COR_CuAg", [37.6, 3.8, 12.2, 1.0])
    # Need to leave out Cu(110) from Huang --> OD-Cu rgh=20, others 3.5


def get_USHE(data):
    urhe = data['V_RHE']
    ushe = deepcopy(urhe)
    ushe[:,0] -= 0.059 * float(data['pH'])
    return ushe


def compute_CO2RR_partial(dat):
    nel = dict(CO=2, C2H4=10, EtOH=10, PropAldh=10, Prop=12, CH4=8, H2=2, \
        HCOO=2, Ac=6, H3CHO=8, C2H6=12, MeOH=6, HOH2CCHO=6, C2H4O2H2=8, \
        C3H5OH=10, C2H5CHO=10, Me2CO=10, MeCOCH2OH=8)
    tlm = {'Methane':'CH4', 'Ethylene':'C2H4', 'CO':'CO', 'Hydrogen':'H2','Ethane':'C2H6', 'Methanol':'MeOH',\
        'Formate':'HCOO', 'Ethanol':'EtOH', 'Acetate':'Ac', 'n-propanol':'Prop', 'Ac':'Ac', 'Acetaldehyde':'H3CHO',\
        'Glycolaldehyde':'HOH2CCHO', 'Ethylene-Glycol':'C2H4O2H2', 'Allyl-Alcohol':'C3H5OH', 'Propionaldehyde':'C2H5CHO',\
        'Acetone':'Me2CO', 'Hydroxyacetone':'MeCOCH2OH'}
    jco2rr = np.zeros(dat['Ethanol'][:,0].size)
    for ads in tlm:
        if ads in dat:
            pj = dat['I(mA/cm2)'][:,0] * dat[ads][:,0]
            jco2rr += pj

    return jco2rr


def read_degr_data():
    # TODO: continue here!
    df = pd.read_csv('../experimental_reference/buonsanti.csv') # roughle at -1.5 V vs. SHE
    # --> can't fit this, but should consistently determine it from total J
    print(df)
    print(df.columns.values)
    assert False
    d_sCO_t = {}; d_rgh_t = {}; d_sCO_rgh = {}
    for NPs in [16, 41, 65]:
        d_sCO_t.update({"NP %inm"%NPs:np.array([df['Time (Hours)'], df['S_CO %i'%NPs]*100.]).T})
        rgh_r = df['J C2H4 %i'%NPs]/df['J C2H4 16'][0]
        d_rgh_t.update({"NP %inm"%NPs:np.array([df['Time (Hours)'], rgh_r]).T})
        d_sCO_rgh.update({"NP %inm"%NPs:np.array([rgh_r, df['S_CO %i'%NPs]*100.]).T})


def get_CuZn_alloy_CO():
    dat = _load_data_paper(["Koper_CO2R_CuZn"])
    #dat = {_rename_key(dat[k]):dat[k] for k in dat}
    
    d_sel = {}; d_rgh = {}; d_rghI = {}
    for k in dat:
        if 'CO' in dat[k]:
            # compute CO selectivity:
            sel_CO = compute_CO_selectivity(dat[k])
            U_SHE = dat[k]['V_RHE'][:,0] - (0.059*float(dat[k]['pH']))
            
            d_sel.update({k:np.array([U_SHE, sel_CO]).T})
            #d_sel.update({k:np.array([dat[k]['V_RHE'][:,0], sel_CO]).T})
            d_rgh.update({k:float(dat[k]['roughness'])})
        
            # no current data
           #dj = {k:np.array([get_USHE(dat[k])[:,0], dat[k]['I(mA/cm2)'][:,0]]).T for k in dat}
           #icrgh, icrgh_std = _fit_average_exp_rgh(dj)
           #d_rghI.update({k:[icrgh, icrgh_std]})


    ### checkin
 #  for k in dat:
 #      print(k, dat[k]['catalyst'], dat[k]['roughness'], d_sel[k][2,:])

    return d_sel, d_rgh



if __name__ == "__main__":
    
  # selCOr, selCOir = load_CO_Cu_data()
  # print(selCO)
  # print("rgh from current is be:",[1., 1.248, 2.70, 3.97])

  # selAcdH = load_Acdh_Cu_data()
  # print(selAcdH)

  # selORR = load_ORR_Pt_data() 
  # print(selORR)
    
  # selH2CO = load_H2CO_Pt_data()
  # print(selH2CO)

  # selAcCuPd = load_Ac_alloy_data("Ji_COR_CuPd", -1.4)
  # print(selAcCuPd)
    
  # selAcCuAg = load_Ac_alloy_data("Sargent_COR_CuAg", -1.6, ulim=300)
  # print(selAcCuAg)

  # selCO, rghCO, pHCO = load_CO_Cu_data_collection()
  # print(selCO)

    give_rgh_estimates()
    
    d_sel, d_rgh = get_CuZn_alloy_CO()



