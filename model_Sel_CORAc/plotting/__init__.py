""" 
    helper functions for plots

"""

#!/usr/bin/env python

import sys, os, logging
import numpy as np

import matplotlib.pyplot as plt
from scipy.constants import golden_ratio
from matplotlib import rcParams


def set_plotting_env(width=3.37,height=3.37/ golden_ratio *1.5/2,\
                    lrbt=[0.135,0.80,0.25,0.95],fsize=9.0):
    """
      function to set some default parameters for formatting of standard plots

      Parameters
      ----------
      width : float
          width of the plot in inches
      height : float
          height of the plot in inches
      lrbt : list/tuple with 4 entries
          left, right, bottom, top (fractional) margin of plot
      fsize : float
          font size in the plot
 
    """
    # set plot geometry
    rcParams['figure.figsize'] = (width, height) # x,y
    rcParams['font.size'] = fsize
    rcParams['figure.subplot.left'] = lrbt[0]  # the left side of the subplots of the figure
    rcParams['figure.subplot.right'] = lrbt[1] #0.965 # the right side of the subplots of the figure
    rcParams['figure.subplot.bottom'] = lrbt[2] # the bottom of the subplots of the figure
    rcParams['figure.subplot.top'] = lrbt[3] # the bottom of the subplots of the figure

    rcParams['xtick.top'] = True
    rcParams['xtick.direction'] = 'in'
    rcParams['ytick.right'] = True
    rcParams['ytick.direction'] = 'in'

    rcParams['legend.fancybox'] = False
    rcParams['legend.edgecolor'] = 'k'


clrs = {'orange':[0.8901960784313725, 0.4470588235294118, 0.13333333333333333],
        'darkgray':[0.34509803921568627, 0.34509803921568627, 0.35294117647058826],
        'lightblue':[0.596078431372549, 0.7764705882352941, 0.9176470588235294],
        'lightblue2':[0.39215686274509803, 0.6274509803921569, 0.7843137254901961],
        'lightblue3': [0.0, 0.6, 1.0],
        'deepblue':[0.0, 0.396078431372549, 0.7411764705882353],
        'azurblue':[0., 0.833, 1.],
        'darkblue':[0., 0., 0.5], 
        'byellow':[1., 0.90, 0.], 
        'darkyellow':[1.0, 0.706, 0.0],
        'orange2':[0.8980392156862745, 0.20392156862745098, 0.09411764705882353],
        'darkred':[0.5, 0., 0.],
        'green':[0.6352941176470588, 0.6784313725490196, 0.0],}


def writefig(filename, folder='output',  write_pdf=True, write_eps=False, write_png=False):
    """
      wrapper for creating figures
      
      Parameters
      ----------
      filename : string
        name of the produced figure (without extention)
      folder : string 
        subfolder in which to write the figure (default = output)
      write_eps : bool
        whether to create an eps figure (+ the usual pdf)

    """
    # folder for output
    if not os.path.isdir(folder):
        os.makedirs(folder)

    fileloc = os.path.join(folder, filename)
    if write_pdf:
        logging.info("writing {}".format(fileloc+'.pdf'))
        plt.savefig(fileloc+'.pdf')
    
    if write_eps:
        logging.info("writing {}".format(fileloc+'.eps'))
        plt.savefig(fileloc+'.eps')
    
    if write_png:
        logging.info("writing {}".format(fileloc+'.png'))
        plt.savefig(fileloc+'.png', dpi=300)
    
    plt.close(plt.gcf())


