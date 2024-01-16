""" 
    Module containing functions to solve kinetic model via ODE solver

"""

import sys, os, logging
import numpy as np
from scipy.integrate import odeint

from model_COR_Ac.mkm.mkm_energetics import adjust_SHE_engs, make_rates


class ODE_mkm_solver(object):
    """ Object for solving microkinetic model (mkm) via scipy ODEint
        solves for steady-state solution (can be changed)
        uses (free) energies of rct-system and takes U_SHE and pH as
        variables to solve model
      
      Parameters
      ----------
      rct_equ : function
        function containing the ODE's to solve (dy/dx = ...) of mkm
        takes c (concentrations), t (time), *args (rate constants)
        see model_surfrct.py
        [compare https://sam-dolan.staff.shef.ac.uk/mas212/notebooks/ODE_Example.html]
      echem : list of booleans
        Boolean list to indicate which steps are electrochemical
      engs : list of lists / array [n,3] (n = number of rcts)
        energies for mkm (to be transferred to rate constants)
      n_int : int
        number of intermediates (non-empty cov)
      p_ind : list
        list of indices of coverages leading to products
      r_ind : list
        list of indices of rates leading to products
    
    """
    def __init__(self, rct_equ, echem, engs, n_int, p_ind, r_ind):
        self.rct_equ = rct_equ # function representing the kinetic model to be solved via ODE solver
        self.echem = echem # Boolean list to indicate which steps are electrochemical
        self.engs0 = engs # energies for kinetic model
        self.n_int = n_int # number of intermediates (non-empty cov)
        self.p_ind = [p for p in p_ind] # coverages leading to products
        self.r_ind = r_ind # rates leading to products
        assert len(p_ind) == len(r_ind)

    
    def _solve_ODE(self, tend, y0, rts):
        """
          helper-method to solve system with rates `rts` and 
          starting coverated `y0` until tend (seconds)

        """
        ts = np.linspace(0,tend,100)
        ys = odeint(self.rct_equ, y0, ts, args=tuple(rts))
        ys = np.array(ys)
        return(ts, ys)

    
    def _converge_time(self, y0, rts):
        """
          helper-method to use _solve_ODE until steady-state

        """
        for te in np.arange(-4,14)[::-1]:
            ts, ys = self._solve_ODE(10**(-1.*te), y0, rts)
            dy = np.array([np.gradient(ys[:,i])/ys[:,i] \
                for i in range(ys[0,:].size)]).T
            l = dy.shape[0]//2
            if np.all(np.absolute(dy[l:,:]) < 1e-4) and not np.all(dy == 0.0):
                return(10**(-1.*te))
        # return last in list
        logging.debug("MKMinner: warning, no steady state convergence")
        return(10**(-1.*te))

    
    def run_kinetic_ODE(self, y0, U_SHE, pH, p0=[]):
        """
          method to obtain steady-state solution using ODE solver
 
          Parameters
          ----------
          y0 : array
            array of starting concentrations (usually np.zeros(self.n_int+1))
          U_SHE : float
            value of applied potential vs. SHE
          pH : float
            value of applied pH
          p0 : list
            list of surface concentrations to be used in mkm (defined in rct_equ)
 
          Returns
          -------
          ts : array 
            time evolution until steady-state
          ys : array
            coverages at timesteps (time, cov); y[-1,:] for steady-state
          pR : array
            production rate as defined through r_ind and p_ind
        
        """
        engsU = adjust_SHE_engs(self.engs0, U_SHE, self.echem)
        rts = make_rates(engsU, pH, self.echem, PD_alk=True)

        rts = np.array(rts.tolist() + p0)
 
        # solve coupled ODE
        tend = self._converge_time(y0, rts) # get time when steady-state
        ts, ys = self._solve_ODE(tend, y0, rts)

 
        # give steady-state production rate; final concentrations hard-coded
        pR = [ys[-1,self.p_ind[i]] * rts[self.r_ind[i]] \
            for i in range(len(self.r_ind))]
 
        return(ts, ys, pR)


if __name__ == "__main__":
    pass


