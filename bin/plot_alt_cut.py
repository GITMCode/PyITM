#!/usr/bin/env python3
""" This is a simple lat/lon plotter
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse

import sys
sys.path.insert(0,'/home/ridley/Software/PyITM/')

from pyitm.fileio import gitmio
from pyitm.modeldata import utils
from pyitm.plotting import plotutils

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Plot Aether / GITM model results')
    
    # select altitude to plot:
    parser.add_argument('-alt', metavar = 'alt', default = 400.0, type = float, \
                        help = 'altitude :  alt in km (closest)') 

    # variable to plot as a number
    parser.add_argument('-ivar',  \
                        default = 3, type = int, \
                        help = 'variable to plot (number)')

    # User can set max and min of the colorbar:
    parser.add_argument('-mini',  default = 1e32, type = float, \
                        help = 'manually set the minimum value for the plots')
    parser.add_argument('-maxi',  default = -1e32, type = float, \
                        help = 'manually set the maxiumum value for the plots')
    
    parser.add_argument('-latmin',  default = -90, type = float, \
                        help = 'manually set the minimum latitude for the plots')
    parser.add_argument('-latmax',  default = 90, type = float, \
                        help = 'manually set the maxiumum latitude for the plots')
    parser.add_argument('-lonmin',  default = 0, type = float, \
                        help = 'manually set the minimum longitude for the plots')
    parser.add_argument('-lonmax',  default = 360, type = float, \
                        help = 'manually set the maxiumum latitude for the plots')

    # Get the files to plot:
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')

    args = parser.parse_args()

    return args


# Needed to run main script as the default executable from the command line
if __name__ == '__main__':

    # Get the input arguments
    args = get_args()
    filelist = args.filelist
    varToPlot = [args.ivar]
    altGoal = args.alt

    allData = gitmio.read_gitm_all_files(filelist, varToPlot)

    alts1d = allData['alts'][0, 0, :]
    diff = np.abs(alts1d - altGoal)
    iAlt = np.argmin(diff)
    realAlt = alts1d[iAlt]

    lons1d = allData['lons'][:, 0, iAlt]
    lats1d = allData['lats'][0, :, iAlt]
    lonsEdge = utils.move_centers_to_edges(lons1d)
    latsEdge = utils.move_centers_to_edges(lats1d)

    allSlices = utils.data_slice(allData, iAlt = iAlt)
    allTimes = allData['times']

    # get min and max values, plus color table:
    dataMinMax = plotutils.get_min_max_data(allSlices, lats1d, \
                     yMin = args.latmin, yMax = args.latmax, \
                     color = 'red', \
                     minVal = 1e32, maxVal = -1e32)

    dpi = 120
    varName = allData['vars'][0]
    sVarNum = 'var%03d_' % args.ivar
    sAltNum = 'alt%04d_' % int(realAlt)

    for iTime, uTime in enumerate(allTimes):

        fig = plt.figure(figsize=(10, 5.5), dpi = dpi)
        ax = fig.add_axes([0.07, 0.06, 0.97, 0.9])

        title = uTime.strftime("%d %b %Y %H:%M:%S UT")
        title = title + '; Alt: %.0f km' % realAlt
        value2d = allSlices[iTime, :, :].transpose()
        con = ax.pcolormesh(lonsEdge, latsEdge, value2d, \
                            cmap = dataMinMax['cmap'], \
                            vmin = dataMinMax['mini'], \
                            vmax = dataMinMax['maxi'])
        ax.set_ylim(args.latmin, args.latmax)
        ax.set_xlim(args.lonmin, args.lonmax)
        ax.set_aspect(1.0)
        ax.set_title(title)

        cbar = fig.colorbar(con, ax = ax, shrink = 0.75, pad = 0.02)
        cbar.set_label(varName, rotation=90)

        sTimeOut = uTime.strftime('%y%m%d_%H%M%S')
        outFile = sVarNum + sAltNum + sTimeOut + '.png'
        print(" ==> Writing file : ", outFile)
        fig.savefig(outFile, dpi = dpi)
        plt.close(fig)
