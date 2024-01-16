""" 
    Module containing the diffrct1d2cmp object for numerical solutions of a 
    1D-reaction-diffusion system with 2 components

"""

import sys, os
import numpy as np
from copy import deepcopy
from time import time
import warnings, logging


class diffrct1d2cmp(object):
    """ 
        Object for numerical solution of 1D-reaction-diffusion system
        with 2 components; all units are SI units in this module

      Parameters
      ----------
      D : list
        list of diffusion constants of species
      kexp : float
        reaction constant for solution reaction
      Lx : float
        diffusion length in m 
      nV : int
        number of discrete elements
      dx : float
        length of elements
      dt : float
        timestep for solver
      imax : int
        max number of solver iterations
      log : bool
        flag for logger
      tol : float
        tolerance for numerical solver to converge
    
    """
    def __init__(self, D, kexp, Lx, nV=50, dx=None, dt=None, imax=10000000, log=False, tol=5e-6):
        self.D = D # list with value per species 
        self.kexp = kexp # reaction constant for solution reaction(s)
        self.ns = 2 # number of species (hard coded)
        if dx == None:
            self.nV = nV # number of discretized volume units
        else:
            self.nV = int(Lx // dx)
        self.dx = Lx / self.nV # length of discretized volume units
        self.dt = dt # timestep in s
        if self.dt == None: # if not given assumed stable step width
            #self.dt = 4950.*(self.dx**2.0) #timestep in s
            self.dt = 4950*1e4*(self.dx**2.0)*(5e-9/D[1]) #timestep in s
        self.imax = imax # maximum number of iterations
        self.init_clear_arrays()
        self.log = log
        if tol != None:
            self.tol = tol # convergence tolerance
        else:
            self.tol = 5e-6 # default for long diffusion lengths
    
    
    def init_clear_arrays(self):
        """
          helper-method to clear internal arrays
 
        """
        self.tC = [] # general concentration container
        self.ts = [] # timestep tracker
        self.Cx = None # current concentration


    def _integrate_internal(self):
        """
          helper-method integrating internal elements following a 
          central difference scheme including diffusion and solution reaction
 
        """
        # first step in integration (adjust internal to outer)
        C0 = self.tC[-1][0,:]; C1 = self.tC[-1][1,:]
        Cx = self.Cx
        dt, dx, D, kexp = self.dt, self.dx, self.D, self.kexp
        for i in range(1,self.nV-1):
            # central difference 2nd order
            Cx[0,i] = C0[i] + \
                dt * D[0] * (C0[i+1] - 2* C0[i] + C0[i-1]) / (dx**2.0) \
                - dt * kexp[0] * C0[i] * C1[i]
            Cx[1,i] = C1[i] + \
                dt * D[1] * (C1[i+1] - 2* C1[i] + C1[i-1]) / (dx**2.0) \
                - dt * kexp[1] * C0[i] * C1[i]

    
    def _integrate_endpoints(self):
        """
          helper-method integrating fist and last element following a 
          central difference scheme including diffusion and solution reaction
 
        """
        # integrate endpoints - overwritten by boundary
        Cx = self.Cx
        C0 = self.tC[-1][0,:]; C1 = self.tC[-1][1,:]
        dt, dx, D, kexp = self.dt, self.dx, self.D, self.kexp
        # first point
        Cx[0,0] = C0[0] + dt * D[0] * (C0[1] - C0[0]) / (dx**2.0) - dt * kexp[0] * C0[0] * C1[0]
        Cx[1,0] = C1[0] + dt * D[1] * (C1[1] - C1[0]) / (dx**2.0) - dt * kexp[1] * C0[0] * C1[0]
        # last point
        Cx[0,-1] = C0[-1] + dt * D[0] * (C0[-1] - C0[-2]) / (dx**2.0) - dt * kexp[0] * C0[-1] * C1[-1]
        Cx[1,-1] = C1[-1] + dt * D[1] * (C1[-1] - C1[-2]) / (dx**2.0) - dt * kexp[1] * C0[-1] * C1[-1]


    def converge_num(self, flux, cB, initC=np.array([]), stout=False):
        """
          main routine converging the transport model to a flux and a
          target bulk concentration of B
      
          Parameters
          ----------
          flux : list
            list of target fluxes for both species
          cB : float
            bulk concentrationn of B (boundary condition)
          initC : np.array
            initial guess for solution
          stout : bool
            flag for output
 
          NOTE the out-commented try/except statements supress
          warnings

        """
        #warnings.filterwarnings('error') # catch bad timesteps
        # init arrays
        self.init_clear_arrays()
        if initC.size == 0:
            initC = np.ones((2,self.nV))
            initC[0,:] *= cB[0]; initC[1,:] *= cB[1]
        self.tC = [initC]; self.ts = [0]
        
        dt, dx, D, kexp = self.dt, self.dx, self.D, self.kexp
        t0 = time()
        for i in range(0, self.imax):
            self.Cx = np.zeros((2,self.nV))
         #### - catch errornous integration/boundaries
         ###try:
            # integrate
            self._integrate_internal()
            self._integrate_endpoints()
            
            # enforce boundary conditions
            # (1) bulk
            self.Cx[0,-1] = cB[0]
            self.Cx[1,-1] = cB[1]
            # (2) source
            # full flux as global/boundary condition
            fn = np.trapz(self.Cx[0,:]*self.Cx[1,:]*kexp[0], dx=dx)
            self.Cx[0,0] += (flux[0]-fn)/(kexp[0]*dx*cB[1]) # necessary correction
            # for now use derivative (ignore OH consumption in first step)
            self.Cx[1,0] = (flux[1]/D[1])*dx + self.Cx[1,1]
            
         ###except RuntimeWarning: # abort problematic integration
         ###    self.Cx[:] = np.nan
         ###    break
            
          ##self.tC.append(self.Cx) # log starting vector
          ##self.ts.append(self.ts[-1]+dt) # timestep
            self.tC = [self.tC[-1], self.Cx] # log starting vector
            self.ts = [self.ts[-1]+dt] # timestep
           
            # maximum deviation from previous run
            dev_prev = max([np.absolute(self.tC[-1][i,:]-self.tC[-2][i,:]).max() \
                for i in range(self.Cx.shape[0])])
            #if dev_prev < 1e-12: # 1e-12 mol/cm3 as minium (for 5-50 mu Lx)
            if dev_prev < self.tol: # 1e-6 mol/m3 as minium (for 5-50 mu Lx) | 5e-6 much faster - and normally enough but not at low fluxes
                break
        #warnings.resetwarnings()
    
        if i == self.imax-1:
            if stout:
                print('not converged - maximum number of steps reached')
            elif self.log:
                logging.debug("DIFFNUM: not converged - maximum number of steps reached")
            self.Cx = np.zeros((2,self.nV))
        elif stout:
            print('converged within %i steps (%.1fs)'%(i, time()-t0))
            #print('flux condition: %.3e vs. %.3e; diff condition: %.3e vs. %.3e'%(fn, fluxK, (self.tC[-1][1,0] - self.tC[-1][1,1])/dx*D_oh, fluxOH))
        elif self.log:
            logging.debug("DIFFNUM: converged within %i steps (%.1fs)"%(i, time()-t0))


    def get_concentration_profile(self):
        """
          helper-method returning the final result
 
        """
        # return last result
        if not np.all(np.isnan(self.Cx)):
            assert np.array_equal(self.Cx, self.tC[-1])
        return({'x':np.linspace(0,self.dx*self.nV,self.nV), \
            'conc':self.Cx})
    

if __name__ == "__main__":
    pass


