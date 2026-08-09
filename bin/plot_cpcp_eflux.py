#!/usr/bin/env python3
""" read in 2d electrodynamic files and plot CPCP and Hemispheric Power
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse

from pyitm.fileio import util, controlfile, logfile
from pyitm.general import geometry
from pyitm.plotting import line_plots, axes

# ----------------------------------------------------------------------------
# Get arguments as inputs into the code
#-----------------------------------------------------------------------------

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Plot CPCP and Hemispheric Power from electrodynamics')
    
    parser.add_argument('-list',  \
                        action='store_true', default = False, \
                        help = 'list variables in file')
    
    parser.add_argument('-plotfile',
                        help = 'output file for plot',
                        default = 'cpcp_hp.png')

    # Get the files to plot:
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------------
# Needed to run main script as the default executable from the command line
# ----------------------------------------------------------------------------

if __name__ == '__main__':

    # Get the input arguments
    args = get_args()
    plt.rcParams.update({'font.size': 14})

    filelist = args.filelist
    varToPlot = ['cpcp', 'eflux']

    allData3D = util.read_all_files(filelist, varToPlot)

    lons2d = allData3D['lons'][:, :, 0]
    lats2d = allData3D['lats'][:, :, 0]
    Altitude = 120.0
    area2d = geometry.calc_areas(lons2d, lats2d, Altitude)
    
    if (np.max(lats2d) > 80.0):
        containsNorth = True
        maskNorth = lats2d > 0.0
    if (np.min(lats2d) < -80.0):
        containsSouth = True
        maskSouth = lats2d < 0.0
        
    logData = { 'times': [], \
                'alt': Altitude, \
                'northCPCP': [], \
                'southCPCP': [], \
                'northPower': [], \
                'southPower': []}
                
    for iTime, time in enumerate(allData3D['times']):
        pot2d = allData3D['data'][iTime, 0, :, :, 0]
        eflux2d = allData3D['data'][iTime, 1, :, :, 0]
        logData['times'].append(time)
        if (containsNorth):
            logData['northCPCP'].append( \
                (np.max(pot2d[maskNorth]) - \
                np.min(pot2d[maskNorth]))/1000.0)
            logData['northPower'].append( \
                np.sum(area2d[maskNorth] * eflux2d[maskNorth]) / \
                (1000.0 * 1e9 * 1e4))
        else:
            logData['northCPCP'].append(0.0)
            logData['northPower'].append(0.0)
        if (containsSouth):
            logData['southCPCP'].append( \
                (np.max(pot2d[maskSouth]) - \
                np.min(pot2d[maskSouth]))/1000.0)
            logData['southPower'].append( \
                np.sum(area2d[maskSouth] * eflux2d[maskSouth]) / \
                (1000.0 * 1e9 * 1e4))
        else:
            logData['northCPCP'].append(0.0)
            logData['northPower'].append(0.0)

    logfile.write_log(logData, fileHeader = 'cpcp_hp', \
                      message = 'CPCP and HP extracted using plot_cpcp_eflux')

    dpi = 120
    fig = plt.figure(figsize=(10, 10), dpi = dpi)

    yBot = 0.06
    yTop = 0.05
    yBuf = 0.06
    ax = axes.get_axes_one_column(fig,
                                  2,
                                  yBot,
                                  yTop,
                                  yBuf)
    outFile = 'cpcp_hp_' + logData['times'][0].strftime('%Y%m%d') + '.png'
    line_plots.lineplot_data(logData, vars = ['northCPCP'],
                             fig = fig, \
                             ax = [ax[0]], \
                             linewidth = 1.0, \
                             color = 'k', \
                             label = 'North CPCP', \
                             linestyle = None,
                             ylabel = 'Cross Polar Cap Potential (kV)')
    line_plots.lineplot_data(logData, vars = ['southCPCP'],
                             fig = fig, \
                             ax = [ax[0]], \
                             linewidth = 1.0, \
                             color = 'r', \
                             linestyle = None,
                             label = ['South CPCP'], \
                             ylabel = 'Cross Polar Cap Potential (kV)', \
                             title = None)
    line_plots.lineplot_data(logData, vars = ['northPower'],
                             fig = fig, \
                             ax = [ax[1]], \
                             label = ['North HP'], \
                             ylabel = 'Hemispheric Power (GW)', \
                             linewidth = 1.0, \
                             color = 'k')
    line_plots.lineplot_data(logData, vars = ['southPower'],
                             fig = fig, \
                             ax = [ax[1]], \
                             label = ['South HP'], \
                             ylabel = 'Hemispheric Power (GW)', \
                             linewidth = 1.0, \
                             color = 'r', \
                             outFile = outFile)
    
