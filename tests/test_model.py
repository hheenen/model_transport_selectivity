import unittest

import numpy as np
import os


from model_Sel_CORAc.simtools import potsweep_rgh_pH, make_facet_mkin_model_gpaw
from model_Sel_CORAc.tools import load_pickle_file


class TestModel(unittest.TestCase):
    # regression tests if everything works and results are still the same
    
    def test_run_model(self):
        ### fixed model-run to compare
        mkin = make_facet_mkin_model_gpaw(100)
        Us = [-1.7, -1.4]; pHs=[14.3]; rghs = [15]

        dat = potsweep_rgh_pH(Us=Us, pHs=pHs, \
            rghs=rghs, Lx=[0.35e-4, 1e-4], diffmode='diff-rct-num', mkin=mkin, savekey='test_sim_gpaw', tol_diff=1e-8)
        
        ### clean-up data
        ofiles = [f for f in os.listdir('.') if f.split('.')[0] == "dat_test_sim_gpaw"]
        for f in ofiles:
            os.remove('%s'%f)

        ### comparison data
        datr = load_pickle_file("data/test_run_model.pkl")
        
        d = dat[15][14.3]; dr = datr[15][14.3]
        for k in d: #numerical differences based on rounding
            self.assertTrue(np.array_equal(np.around(d[k],4),np.around(dr[k],4)))


if __name__ == "__main__":
    unittest.main()
