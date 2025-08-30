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

# ----------------------------------------------------------------------------
# Get arguments as inputs into the code
#-----------------------------------------------------------------------------

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
    
    parser.add_argument('-latmin',  default = 0, type = float, \
                        help = 'manually set the min latitude for the plots')
    parser.add_argument('-latmax',  default = 0, type = float, \
                        help = 'manually set the max latitude for the plots')

    parser.add_argument('-list',  \
                        action='store_true', default = False, \
                        help = 'list variables in file')
    
    parser.add_argument('-north',  \
                        action='store_true', default = False, \
                        help = 'Plot Northern Hemisphere')
    parser.add_argument('-south',  \
                        action='store_true', default = False, \
                        help = 'Plot Southern Hemisphere')
    
    # Get the files to plot:
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------------
# This assumes you have no blocks, and an easily defined grid
#-----------------------------------------------------------------------------

def plot_sphere(args, allData, plotNorth = True, plotSouth = True):

    alts1d = allData['alts'][0, 0, :]
    diff = np.abs(alts1d - altGoal)
    iAlt = np.argmin(diff)
    realAlt = alts1d[iAlt]

    lons2d = allData['lons'][:, :, iAlt]
    lats2d = allData['lats'][:, :, iAlt]

    allSlices = utils.data_slice(allData, iAlt = iAlt)
    varName = allData['longname'][0]
    sVarNum = allData['shortname'][0] + '_'
    sAltNum = 'alt%04d_' % int(realAlt)

    allTimes = allData['times']

    # get min and max values, plus color table:
    dataMinMax = plotutils.get_min_max_data(allSlices, None, \
                     color = 'red', \
                     minVal = 1e32, maxVal = -1e32)

    sFilePre = sVarNum + sAltNum
    sTitleAdd = '; Alt: %.0f km' % realAlt
    plotters.plot_slices_polar_only(allSlices,
                                    allTimes,
                                    lons2d,
                                    lats2d,
                                    dataMinMax,
                                    plotNorth = plotNorth,
                                    plotSouth = plotSouth,
                                    xLabel = 'Longitude (deg)',
                                    yLabel = 'Latitude (deg)',
                                    varName = varName,
                                    titleAddOn = sTitleAdd,
                                    filenamePrefix = sFilePre,
                                    yLimits = [args.latmin, args.latmax])   

    return

# ----------------------------------------------------------------------------
# This assumes you have blocks, like a cubesphere grid
#-----------------------------------------------------------------------------

def plot_cubesphere(args, allData):

    alts1d = allData['alts'][0, 0, 0, :]
    diff = np.abs(alts1d - altGoal)
    iAlt = np.argmin(diff)
    realAlt = alts1d[iAlt]

    allSlices = utils.data_slice(allData, iAlt = iAlt)
    varName = allData['longname'][0]
    sVarNum = variables.strip_varname(allData['shortname'][0]) + '_'
    sAltNum = 'alt%04d_' % int(realAlt)

    allTimes = allData['times']

    # get min and max values, plus color table:
    dataMinMax = plotutils.get_min_max_data(allSlices, None, \
                     color = 'red', \
                     minVal = 1e32, maxVal = -1e32)
    
    lons3d = allData['lons'][:, :, :, iAlt]
    lats3d = allData['lats'][:, :, :, iAlt]
    
    sFilePre = sVarNum + sAltNum
    sTitleAdd = '; Alt: %.0f km' % realAlt
    plotters.plot_series_of_slices_wblocks(allSlices,
                                           allTimes,
                                           lons3d,
                                           lats3d,
                                           dataMinMax,
                                           xLabel = 'Longitude (deg)',
                                           yLabel = 'Latitude (deg)',
                                           varName = varName,
                                           titleAddOn = sTitleAdd,
                                           filenamePrefix = sFilePre,
                                           yLimits = [args.latmin, args.latmax])

# Needed to run main script as the default executable from the command line
if __name__ == '__main__':

    # Get the input arguments
    args = get_args()
    filelist = args.filelist
    varToPlot = args.var

    plotNorth = args.north 
    plotSouth = args.south

    if ((not plotNorth) and (not plotSouth)):
        print('You need to select at least one of -north or -south')
        exit()
        
    if (args.list):
        util.list_file_info(filelist)
        exit()
    
    allData = util.read_all_files(filelist, varToPlot)

    if (not allData):
        util.list_file_info(filelist)
        exit()

    altGoal = args.alt
    
    if (allData['nblocks'] == 0):
        plot_sphere(args, allData, plotNorth, plotSouth)
        
    if (allData['nblocks'] > 0):
        plot_cubesphere(args, allData, plotNorth, plotSouth)
        
