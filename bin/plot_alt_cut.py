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

from pyitm.fileio import util
from pyitm.modeldata import utils
from pyitm.plotting import plotutils
from pyitm.plotting import plotters

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Plot Aether / GITM model results')
    
    # select altitude to plot:
    parser.add_argument('-alt', metavar = 'alt', default = 400.0, type = float, \
                        help = 'altitude :  alt in km (closest)') 

    # variable to plot as a number
    parser.add_argument('-var',  \
                        default = '3', \
                        help = 'variable to plot (number or variable name)')

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
    varToPlot = args.var

    allData = util.read_all_files(filelist, varToPlot)

    altGoal = args.alt
    alts1d = allData['alts'][0, 0, :]
    diff = np.abs(alts1d - altGoal)
    iAlt = np.argmin(diff)
    realAlt = alts1d[iAlt]

    lons1d = allData['lons'][:, 0, iAlt]
    lats1d = allData['lats'][0, :, iAlt]
    lonsEdge = utils.move_centers_to_edges(lons1d)
    latsEdge = utils.move_centers_to_edges(lats1d)

    allSlices = utils.data_slice(allData, iAlt = iAlt)
    varName = allData['longname'][0]
    sVarNum = allData['shortname'][0] + '_'
    sAltNum = 'alt%04d_' % int(realAlt)

    allTimes = allData['times']

    # get min and max values, plus color table:
    dataMinMax = plotutils.get_min_max_data(allSlices, lats1d, \
                     yMin = args.latmin, yMax = args.latmax, \
                     color = 'red', \
                     minVal = 1e32, maxVal = -1e32)

    sFilePre = sVarNum + sAltNum
    sTitleAdd = '; Alt: %.0f km' % realAlt
    plotters.plot_series_of_slices(allSlices,
                                   allTimes,
                                   lonsEdge,
                                   latsEdge,
                                   dataMinMax,
                                   xLabel = 'Longitude (deg)',
                                   yLabel = 'Latitude (deg)',
                                   varName = varName,
                                   titleAddOn = sTitleAdd,
                                   filenamePrefix = sFilePre,
                                   xLimits = [args.lonmin, args.lonmax],
                                   yLimits = [args.latmin, args.latmax])   
