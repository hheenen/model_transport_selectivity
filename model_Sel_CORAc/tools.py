""" 
    Module containing some basic functions
"""

import pickle
import hashlib, base64
import numpy as np


def load_pickle_file(filename,py23=False):
    pickle_file = open(filename,'rb')
    if py23:
        data = pickle.load(pickle_file, encoding='latin1')
    else:
        data = pickle.load(pickle_file)
    pickle_file.close()
    return(data)


def write_pickle_file(filename,data,py2=False):
    ptcl = 0
    if py2:
        filename = filename.split('.')[0]+'_py2.'+filename.split('.')[1]
        ptcl = 2
    output = open(filename, 'wb')
    pickle.dump(data, output, protocol=ptcl)
    output.close()


def _make_hash(array):
    """
      helper-function to make hash from array

    """
    array[np.where(array == 0)] = 0
    hashed = hashlib.md5(array).digest()
    hashed = base64.urlsafe_b64encode(hashed)
    return(hashed)


