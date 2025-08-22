#!/usr/bin/env python3
""" This is a simple lat/lon plotter
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse
import glob

import sys
sys.path.insert(0,'/home/ridley/Software/PyITM/')

from pyitm.fileio import gitmio
from pyitm.modeldata import utils
from pyitm.plotting import plotutils
from pyitm.plotting import plotters

# Needed to run main script as the default executable from the command line
if __name__ == '__main__':

    filelist = glob.glob('x*/2DGEL_t021221_000500.bin')

    iPed_ = 4
    iAvee_ = 6
    iEflux_ = 7
    varsToRead = [iPed_, iAvee_, iEflux_]
    
    allData = gitmio.read_gitm_all_files(filelist, varsToRead)

    ped = allData['data'][:, 0, :, :, :]
    avee = allData['data'][:, 1, :, :, :]
    eflux = allData['data'][:, 2, :, :, :]

    avee1d = avee[eflux > 0.5]
    ped1d = ped[eflux > 0.5]
    eflux1d = eflux[eflux > 0.5]

    pedRob = 40.0 * avee1d / (16.0 + avee1d**2) * np.sqrt(eflux1d) 
    
    dpi = 200
    fig = plt.figure(figsize=(10, 10), dpi = dpi)
    ax = fig.add_axes([0.07, 0.06, 0.9, 0.9])
    ax.scatter(avee1d, ped1d, color = 'r', label = 'GITM')
    ax.scatter(avee1d, pedRob, color = 'b', label = 'Robinson')
    ax.set_ylabel('Pedersen Conductance (mhos)')
    ax.set_xlabel('Average Energy (keV)')
    ax.legend()
    
    outFile = 'ped_vs_avee.png'
    print(" ==> Writing file : ", outFile)
    fig.savefig(outFile, dpi = dpi)
    plt.close(fig)


    fig = plt.figure(figsize=(10, 10), dpi = dpi)
    ax = fig.add_axes([0.07, 0.06, 0.9, 0.9])
    ax.scatter(pedRob, ped1d, color = 'k')
    ax.plot([0,30],[0,30])
    
    ax.set_ylabel('GITM Ped. Conductance (mhos)')
    ax.set_xlabel('Robinson Ped. Conductance (mhos)')
    
    outFile = 'gitm_vs_rob.png'
    print(" ==> Writing file : ", outFile)
    fig.savefig(outFile, dpi = dpi)
    plt.close(fig)
    
