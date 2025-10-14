#!/usr/bin/env python3
""" This is for satellite plots.
"""

import argparse
import os
import numpy as np
from pyitm.fileio import satelliteio
from pyitm.fileio.util import read_all_files
from pyitm.modeldata import satellite
import matplotlib.pyplot as plt

def makesatplot(satDataDict, savepath, verbose=False, saveName=None):

    fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    didSmooth = 'smoothed_sat_rho' in satDataDict.keys() 
    
    if didSmooth:
        axs[0].plot(satDataDict['times'], satDataDict['model_rho'], 'b', 
                    alpha = 0.1)
        axs[0].plot(satDataDict['times'], satDataDict['sat_rho'], 'r', 
                    alpha = 0.1)
        axs[0].plot(satDataDict['times'], satDataDict['smoothed_model_rho'], 'b', 
                    label="Modeled Density")
        axs[0].plot(satDataDict['times'], satDataDict['smoothed_sat_rho'], 'r', 
                    label=satDataDict['sat_name'].upper() + " Density")
    else:
        axs[0].plot(satDataDict['times'], satDataDict['model_rho'], 'b', 
                    label="Modeled Density")
        axs[0].plot(satDataDict['times'], satDataDict['sat_rho'], 'r', 
                    label=satDataDict['sat_name'].upper() + " Density")
    axs[0].legend()

    
    if didSmooth:
        axs[1].plot(satDataDict['times'], 
                    100*(satDataDict['sat_rho']- satDataDict['model_rho'])
                    /satDataDict['sat_rho'],
                    alpha = 0.1, color='k')
        axs[1].plot(satDataDict['times'], 
                    100*(satDataDict['smoothed_sat_rho']-satDataDict['smoothed_model_rho'])
                    /satDataDict['smoothed_sat_rho'],
                    color='k')
    else:
        axs[1].plot(satDataDict['times'], 
                    100*(satDataDict['sat_rho']- satDataDict['model_rho'])
                    /satDataDict['sat_rho'],
                    color='k')
    
    axs[1].hlines(0, min(satDataDict['times']), max(satDataDict['times']),
                  linestyle='--', color='k')

    axs[1].set_ylabel(' Density (kg/m3)')
    axs[1].set_ylabel(' Diff (%)')
    fig.suptitle(satDataDict['sat_name'].upper())

    i = 0
    outdate = satDataDict['times'][0].strftime('%Y%m%d')
    savename = os.path.join(savepath, f"{satDataDict['sat_name']}_density_{outdate}_000.png")
    
    old_end = str(i).rjust(3, '0') + '.png'
    while os.path.exists(savename):
        # start adding zeros!
        i += 1
        new_end = str(i).rjust(3, '0') + '.png'
        savename = savename.replace(old_end, new_end)
        old_end = new_end
    
    print('--> Saving plot as: ' + savename)
    plt.savefig(savename)

    return
