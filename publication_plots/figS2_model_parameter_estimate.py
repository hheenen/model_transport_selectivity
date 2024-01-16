#!/usr/bin/env python

"""

This script only produces 1D cuts for each parameter

"""

from header import *


def screen_model_parameters(param_std, ranges):
    """
      function to run model with different model parameters
      # parameters are taken from global dict `param_std`
      and `ranges`
 
      Parameters
      ----------
        param_std : dict
          dictionary with standard parameters hold constant during screening
        ranges : dict
          dictionary with ranges for screening

      Returns
      -------
        out : dict 
          complex output {model:{param:[val_param, val_selectivity]}}
        pstd_plot : dict 
          complex output of standard values {param:[val_param, val_selectivity]}}
    
    """

    out = {}
    for mdl in [1,2]: # include model dependence in parameter screening
        out_plot = {}
        for p in ['cD', 'Lx', 'dG', 'ddG_des', 'ddG_red', 'rgh']:
            if p not in ['dG', 'ddG_des']: # complex workaround
                dGdes , dGads = get_adsdes_engs(param_std)
                dGdes = [dGdes]; dGads = [dGads]
            else:
                dGdes = []; dGads = []
                for e in ranges[p]:
                    pstd = deepcopy(param_std); pstd.update({p:[e]})
                    ddGdes , ddGads = get_adsdes_engs(pstd)
                    dGdes.append(ddGdes); dGads.append(ddGads)
                    
 
            params = deepcopy(param_std)
            params.update({p:ranges[p]})
 
            datafile = "model_data_screen_%s.pkl"%p
            # output: i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1
            if p not in ['dG', 'ddG_des']:
                dat = sample_data(datafile, rdes=dGdes, rads=dGads, 
                    rred=params['ddG_red'], Ds=params['cD'], Lxs=params['Lx'], 
                    Us=params['U'], rghs=params['rgh'], mdls=[mdl])
            else: # complex workaround for matched values
                dat = np.zeros((len(params[p]), 16))
                for i in range(len(params[p])):
                    dat[i,:] = sample_data(datafile, rdes=[dGdes[i]], rads=[dGads[i]], 
                        rred=params['ddG_red'], Ds=params['cD'], Lxs=params['Lx'], 
                        Us=params['U'], rghs=params['rgh'], mdls=[mdl])[0].tolist()
 
            # save roughness vs. sel (in %)
            sel = (dat[:,13]/ dat[:,[13,14]].sum(axis=1)) * 100.
            # print(p, params[p], sel)
 
 
            params[p] = adjust_param_for_plot(p, params[p])
            
            out_plot.update({p:np.array([params[p], sel]).T})
        
        out.update({mdl:out_plot})


    pstd_plot = {p:adjust_param_for_plot(p, deepcopy(param_std[p])) for p in param_std}

    return pstd_plot, out


def adjust_param_for_plot(p, param):
    """
      helper-function to adjust units for plotting

    """
    if p == 'dG':
        param = np.array(param) * -1 # let's call it adsorption free energy (previously was desorption
    elif p == 'ddG_red':
        param = np.array(param) * param_std['U']*0.5
    elif p == 'Lx':
        param = np.array(param) * 1e4
    elif p == 'cD':
        param = np.array(param) * 1e10
    return param


def plot_model_parameter_selectivity_1D(filename, dat, pstd, datm2={}):
    """
      function to make complex plot
      
      Parameters
      ----------
        dat : dict 
          complex output {param:[val_param, val_selectivity}}
        pstd : dict 
          complex output of standard values {param:[val_param, val_selectivity]}}
        datm2 : dict 
          complex output {param:[val_param, val_selectivity}}

    """
    set_plotting_env(width=3.37*2,height=3.37/ golden_ratio *1.05,\
               lrbt=[0.08,0.98,0.19,0.98],fsize=9.0)
        
    # Reminder - output dat: 
    # i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1

    # rgh, Dx, Lx, Gads, dGdes, dGred*U*0.5
    klabel = {'dG':r'$\Delta G _{\mathrm{ads}}$ (eV)', 
              'ddG_des':r'$\Delta G ^{\ddag} _{\mathrm{des}}$ (eV)',
              'ddG_red':r'$\Delta G ^{\ddag} _{\mathrm{red}} (U, \alpha)$ (eV)', 
              'cD':r'$D_{\mathrm{diff}}$ ($10^{10}$ $\frac{\mathrm{m}^2}{s}$)', #m$^2$ s$^{-1}$)', 
              'Lx':r'$L_{\mathrm{diff}}$ ($\mathrm{\mu}$m)', 
              'rgh':r'roughness $\rho$' }

    tag = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

    fig, ax = plt.subplots(1,6)

    for i, k in enumerate(['dG', 'ddG_des', 'ddG_red', 'cD', 'Lx', 'rgh']):
        d = dat[k]
        ax[i].plot(d[:,0], d[:,1], color='k')
        ax[i].set_xlabel(klabel[k])
        ax[i].set_ylim(-0.02*100., 1.02*100.)
        ax[i].plot(pstd[k], 0.514*100., marker='o', color='k', markersize=2)
        ax[i].annotate(tag[i], xy=(0.1,0.85), xycoords='axes fraction')
        if k in datm2: # other model
            dd = datm2[k]
            ax[i].plot(dd[:,0], dd[:,1], color='k', ls='--')
    
    [ax[i].set_yticklabels([]) for i in range(1,6)]
    plt.subplots_adjust(wspace=0.05)

    if len(datm2) > 0:
        ax[1].annotate(r'model b', xy=(0.7,0.8), xytext=(0.1,0.7), size=7, \
            xycoords='axes fraction', arrowprops={"arrowstyle":'-', "color":'k', 'lw':1.0})

    ax[0].set_ylabel(r'selectivity $C_g$ (%)')
    
    writefig(filename ,folder='output', write_pdf=False, write_png=True) 


def make_plot_model_parameter_selectivity():
    """
      helper-function to make full plot

    """
    pstd_plot, out = screen_model_parameters(param_std, ranges)
    # def out - model 1/2
    # out[2] == direct reduction to C*
    # out[1] == reduction and desorptoin to C_aq
    plot_model_parameter_selectivity_1D('FigS2_model_parameter_sel_1D', \
        out[2], pstd_plot, out[1])
  

def eval_sel_diffusion_coeff():
    """
      helper-function to show selectivities for diffusion coefficient range 
      as mentioned in paper

    """
    # specific parameter output for selectivity change from 13.75e-10 to 20.3e-10
    datafile = "model_data_screen_%s.pkl"%'cD'
    Ds_eval = [13.75e-10, 20.3e-10]
    
    params = deepcopy(param_std)
    dGdes , dGads = get_adsdes_engs(param_std)
    dGdes = [dGdes]; dGads = [dGads]
    dat = sample_data(datafile, rdes=dGdes, rads=dGads, 
        rred=params['ddG_red'], Ds=Ds_eval, Lxs=params['Lx'], 
        Us=params['U'], rghs=params['rgh'], mdls=[2])
    sel = (dat[:,13]/ dat[:,[13,14]].sum(axis=1)) * 100.
    print('Selectivity for diffusion coeffs', Ds_eval)
    print(sel)


def make_levich_estimate():
    '''
      function to make diffusion length estimate from levich equation
      mentioned in supporting info
    
    '''
    # estimate d_levic --> for comparison of tested diffusion layer thicknesses
    D_CO = 2.03e-5 * 1e-4 # cm2/s --> m2/s
    D_range = np.array([1.0, 2.2]) * 1e-5 * 1e-4 # 
    # dynamic --> kinematic viscosities (mPa s) / (g/cm3) = m2/s *1e-6
    visc =  (0.9980/1.04652 +  1.1200/1.091221)/2.0 * 1e-6 # mPa s --> Ns/m2 for KOH 10.1021/je000019h
    v_range = np.array([0.5, 1.5]) * 1e-6 # --> temperature range of water

    rpms = [100, 400, 800, 1600, 2500]
    ds = []
    for rpm in rpms:
        for D in D_range:
            for v in v_range:
                o = rpm * np.pi*2/60
                dlevic = 1.6126 * D_CO**(1./3.) * visc**(1./6.) * o**(-1./2.) *1e2 # in cm
                ds.append(dlevic * 1e-2 *1e6)

    print("Estimate for diffusion length from Levich equation")
    print(min(ds), max(ds), 'in mu')


########################################
### presets to screen model behavior ###
########################################

param_std = dict(
    cH = [1],  # Henry's constant = 1
    cD = [20.0e-10],    # m2/s
    Lx = [50e-4],       # cm (50e-6m)
    dG = [-0.25],
    ddG_des = [0.15],
    ddG_red = [0.815],
    rgh = [1],
    U = [-1.0],
    )


# ranges - rough later refine in case there are sensitive regions with maximum effect
ranges = dict(
    cD = np.arange(10,61,5)*1e-10, # typical D(aq), \cite{Cussler}
    Lx = np.arange(10,100,5)*1e-4,    # ?! need range from levich equation estimate
    dG = np.arange(-0.45, -0.05, 0.05),
    ddG_des = np.arange(0.1,0.5,0.05), # around 0.15 \cite{Heenen2022}
    ddG_red = np.arange(0.65,1.05,0.05),
    rgh = np.array([1,2,3,5,8,10,15,20,25,30,35,40,50,60]),
    )



if __name__ == "__main__":
    make_plot_model_parameter_selectivity()
    eval_sel_diffusion_coeff()
    make_levich_estimate()



