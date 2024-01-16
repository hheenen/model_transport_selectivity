""" 
    physical constants and unit conversion used throughout the code
    hard-coded

"""

####################################################
####    functions and imports for constants     ####
from ase.build import fcc110, fcc111, fcc100, fcc211

def _get_cell_xy(a):
    """
      helper function: get area of slab model with z=vacuum in ASE-atoms object

    """
    A = a.get_volume()/a.cell[2,2]
    return(A)
####################################################



####    pysical constants   ####
c_kB = 8.6173331e-5     # Boltzmann constant [eV/K]
c_Na = 6.02214e23       # Avogadro's constant 
c_qe = 1.602e-19        # elementary charge [C]


####    chemical constants  ####
c_lat_Cu = 3.6149       # lattice constant Cu
d_AperSite = {          # dictionary with area per site
    110:_get_cell_xy(fcc110('Cu', a=c_lat_Cu, size=(1,1,9), vacuum=10.0, orthogonal=True)),\
    111:_get_cell_xy(fcc111('Cu', a=c_lat_Cu, size=(1,1,4), vacuum=10.0, orthogonal=False)),\
    100:_get_cell_xy(fcc100('Cu', a=c_lat_Cu, size=(1,1,5), vacuum=10.0, orthogonal=True)),\
    211:_get_cell_xy(fcc211('Cu', a=c_lat_Cu, size=(3,1,4), vacuum=10.0, orthogonal=True)),\
    }
c_D_kt = 14.0224e-10    # m2/s +/- 3.16e-10 (10.1021/ie201944h0)
c_D_oh = 52.7e-10       # m2/s              (10.1073/pnas.1713164114)
c_D_co = 20.3e-10       # m2/s              (Cussler, E. L. (1997). Diffusion: Mass Transfer in Fluid Systems)
c_k_sol = 5.26e4 * 1e-3 # m3/(mol*s)        (10.1139/v00-032)

c_H_co = 1100           # Henry's constant for CO   (10.5194/acp-15-4399-2015)


####    unit conversion    ####
c_A2cm = 1e-8


