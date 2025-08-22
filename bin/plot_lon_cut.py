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

from pyitm.fileio import gitmio, util
from pyitm.modeldata import utils
from pyitm.plotting import plotutils

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Plot Aether / GITM model results')
    
    # select altitude to plot:
    parser.add_argument('-lon', metavar = 'lon', default = 180.0, type = float, \
                        help = 'longitude in deg (closest, -1 for zonal average)') 

    # variable to plot as a number
    parser.add_argument('-var',  \
                        default = '3', \
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
    parser.add_argument('-altmin',  default = 0, type = float, \
                        help = 'manually set the minimum altitude for the plots')
    parser.add_argument('-altmax',  default = 360, type = float, \
                        help = 'manually set the maxiumum altitude for the plots')

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
    varToPlot = args.var

    allData = util.read_all_files(filelist, varToPlot)

    if (args.lon >= 0):
        lonGoal = args.lon
        lons1d = allData['lons'][:, 0, 0]
        diff = np.abs(lons1d - lonGoal)
        iLon = np.argmin(diff)
        realLon = lons1d[iLon]
    else:
        realLon = -1.0
        iLon = -1

    alts1d = allData['alts'][iLon, 0, :]
    lats1d = allData['lats'][iLon, :, 0]
    altsEdge = utils.move_centers_to_edges(alts1d)
    latsEdge = utils.move_centers_to_edges(lats1d)

    allSlices = utils.data_slice(allData, iLon = iLon)
    varName = allData['longname'][0]
    sVarNum = allData['shortname'][0] + '_'
    
    if (args.lon < 0):
        sLonNum = 'lonZave'
    else:
        sLonNum = 'lon%04d_' % int(realLon)

    allTimes = allData['times']

    # get min and max values, plus color table:
    dataMinMax = plotutils.get_min_max_data(allSlices, alts1d, \
                     yMin = args.altmin, yMax = args.altmax, \
                     color = 'red', \
                     minVal = 1e32, maxVal = -1e32)

    dpi = 120

    for iTime, uTime in enumerate(allTimes):

        fig = plt.figure(figsize=(10, 5.5), dpi = dpi)
        ax = fig.add_axes([0.07, 0.06, 0.97, 0.9])

        title = uTime.strftime("%d %b %Y %H:%M:%S UT")
        if (iLon > -1):
            title = title + '; Lon: %.0f deg' % realLon
        else:
            title = title + '; Zonal Average'
        value2d = allSlices[iTime, :, :].transpose()
        con = ax.pcolormesh(latsEdge, altsEdge, value2d, \
                            cmap = dataMinMax['cmap'], \
                            vmin = dataMinMax['mini'], \
                            vmax = dataMinMax['maxi'])
        ax.set_ylim(args.altmin, args.altmax)
        ax.set_xlim(args.latmin, args.latmax)
        #ax.set_aspect(1.0)
        ax.set_title(title)

        cbar = fig.colorbar(con, ax = ax, shrink = 0.75, pad = 0.02)
        cbar.set_label(varName, rotation=90)

        sTimeOut = uTime.strftime('%y%m%d_%H%M%S')
        outFile = sVarNum + sLonNum + sTimeOut + '.png'
        print(" ==> Writing file : ", outFile)
        fig.savefig(outFile, dpi = dpi)
        plt.close(fig)
