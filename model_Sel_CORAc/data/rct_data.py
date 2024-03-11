""" 
    hard-coded parameters to prepare energetics for mkm-model from
    DFT data stored in catmap-style data units

"""


################################################
####   catmap-notation of reaction network  ####
####   more complex than actually used      ####
####   effective mechanism see `mechs` list ####
################################################

rxn = [
    'CO_g + *_t <-> CO*_t',	                                    #0
   #'2CO_t <-> OCCO_t',                                             #1 - skipped 
   #'OCCO_t + H2O_g + ele_g <-> OCCOH_t + OH_g',                    #2 - skipped
    '2CO_t + H2O_g + ele_g <-> OCCOH_t + OH_g',                     #1
    'OCCOH_t + ele_g <-> OCC_t + OH_g',                             #2 
    'OCC_t + H2O_g + ele_g <-> OCCH_t + OH_g',                      #3
    
    'OCCH_t + H2O_g + ele_g <-> OCCH-H2O-ele_t <-> OCCH2_g + OH_g + 2*_t; beta=0.5',   #4
    'OCCH2_g + 2*_t <-> OCCH2-_t  <-> OCCH2_t',                     #5
    'OCCH2_t + H2O_g + ele_g <-> OCCH2-H2O-ele_t <-> OHCCH2_t + OH_g; beta=0.5',        #6
    'OHCCH2_t + 3H2O_g + 3ele_g <-> CH3CH2OH_g + 2*_t + 3OH_g',     #7

    'OCCH2_g + H2O_g <-> CH3COOH_g',                                #8
    
    'OCCH_t + H2O_g + ele_g <-> HCCHO_t + OH_g; beta=0.5',          #9
    'HCCHO_t + H2O_g + ele_g <-> OHCCH2_t + OH_g',                  #10
    
    'OCCH2_t <-> CH2CO_g + 2*_t',                               #11
    'CH2CO_g + H2O_g <-> CH3COOH_g',                            #12

    '2CO_t + 7H2O_g + 8ele_g <-> CH3CH2OH_g + 2*_t + 8OH_g',    #13
    
    'OCCH2_t + H2O_g + ele_g <-> OCCH2-OH2-ele_t <-> OCCH3_t + OH_g; beta=0.5',        #6 for 111
    'OCCH2_t + H2O_g + ele_g <-> OCCH2-ele-H2O_t <-> HOCCH2_t + OH_g; beta=0.5',        #6 for 211
    ]


####################################################
####  effective mechanism compare rct_data.rxn  ####
####  small difference between facets due to    ####
####  OHCCH2, OCCH3, HOCCH2 as most stable      ####
####################################################

mechs = {111: [0,0,1,2,3,4,5,14,7],
         100: [0,0,1,2,3,4,5,6,7],
         110: [0,0,1,2,3,4,5,6,7],
         211: [0,0,1,2,3,4,5,15,7],
    }
    

####################################################
####  pressures as std use --> note that over   ####
####  written by transport module for CO        ####
####################################################

pressures = {
    'CO_g' :1.0,
    'H2_g' :0.0,
    'H2O_g':0.035,
    'CH3COO_g':0.0,
    }


####################################################
####  function and internal path to energy file ####
####  in catmap format --> can be replaced with ####
####  other energy files                        #### 
####################################################

def read_catmap_input(infile):
    """
      helper-function to read standard catmap file 
      (see https://catmap.readthedocs.io/en/latest/tutorials/generating_an_input_file.html)

    """
    with open(infile,'r') as ifile:
        lines = ifile.readlines()
    dat = {}
    for i in range(1,len(lines)):
        lf = lines[i].split('\t')
        if lf[0] not in dat: #add surface
            dat.update({lf[0]:{}})
        if lf[1] not in dat[lf[0]]: #add facet
            dat[lf[0]].update({lf[1]:{}})
        dat[lf[0]][lf[1]].update({lf[2]:{'fE':float(lf[3]),\
                                'vib':lf[4],'tag':lf[5].strip()}})
    return(dat)


########################################################
####  read catmap-style file and store as `e_dat`   ####
####  file can be replaced by other energetics      ####
########################################################

import model_Sel_CORAc
p = model_Sel_CORAc.__path__[0]
engfile = '%s/data/mkm_files/catin_allfacets_4.400_GPAW_SJM_wbarriers.txt'%p
e_dat = read_catmap_input(engfile)

