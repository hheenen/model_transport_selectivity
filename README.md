Model desorption--re-adsorption--reaction
=========================================

This repository contains the code and data for the transport coupled 
microkinetic model to 
1) simulate a minimal model for the desorption--re-adsorption--reaction
and
2) The acetate selectivity mechanism also contained in https://
doi.org/10.5281/zenodo.5013854.

### Content of the repository

- ``model_Sel_DRAR`` is a python script-collection for the 
general model of the desorption--re-adsorption--reaction
- ``model_Sel_CORAc`` is a python script-collection for the transport model 
to describe acetate selectivity during CORR and according to 
https://doi.org/10.5281/zenodo.5013854
- ``tests`` contains a regression test for the ``model_Sel_CORAc``
- ``publication_plots`` contains scripts to run simulations and plot figures
for the manuscript 10.26434/chemrxiv-2023-cxbgr
- ``literature_data`` contains experimental data in txt files used to 
cmopare the model in the manuscript 10.26434/chemrxiv-2023-cxbgr

### To get started 

The scripts contain minimal dependencies. Besides numpy, scipy, ase-dtu is 
required. 
You can obtain the code by cloning from this page:
```bash
cd PATHTOREPO
git clone https://github.com/hheenen/model_transport_selectivity.git
```
No installation script is provided, you can run the model through the plotting
scripts:
```bash
cd publication_plots
python make_all_plots.py # optionally add error-bars for experimental data
```

### Dependencies

- numpy
- scipy
- ASE (atomic simulation environment)

