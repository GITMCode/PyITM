#!/usr/bin/env python3
""" This is a simple lat/lon plotter with polar dials
"""

import numpy as np
import argparse

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
    
    parser.add_argument('-v',  \
                        action='store_true', default = False, \
                        help = 'set verbose')
    
    # select altitude to plot:
    parser.add_argument('-alt', metavar = 'alt', \
                        default = 400.0, type = float, \
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
                        help = 'set the minimum latitude for the plots')
    parser.add_argument('-latmax',  default = 90, type = float, \
                        help = 'set the maxiumum latitude for the plots')
    parser.add_argument('-lonmin',  default = 0, type = float, \
                        help = 'set the minimum longitude for the plots')
    parser.add_argument('-lonmax',  default = 360, type = float, \
                        help = 'set the maxiumum latitude for the plots')

    # directory to use as a background, so you can subtract one run from another
    parser.add_argument('-backdir',  \
                        default = '', \
                        help = 'background directory to subtract off')
    parser.add_argument('-percent',  \
                        action='store_true', default = False, \
                        help = 'plot percent different (if backdir specified)')

    parser.add_argument('-list',  \
                        action='store_true', default = False, \
                        help = 'list variables in file')
    
    # Get the files to plot:
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')

    parser.add_argument('-iStart',  default = 0, type = int, \
                        help = 'start index for time in multi-time file')
    parser.add_argument('-iEnd',  default = 0, type = int, \
                        help = 'end index for time in multi-time file')
    parser.add_argument('-iStep',  default = 1, type = int, \
                        help = 'step size for time in multi-time file')

    parser.add_argument('-start',  default = None, \
                        help = 'first time to plot (index or date string)')
    parser.add_argument('-stop',  default = None, \
                        help = 'last time to plot (index or date string)')
    parser.add_argument('-time',  default = None, \
                        help = 'only plot the time nearest this (index or date string)')

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------------
# This assumes you have no blocks, and an easily defined grid
#-----------------------------------------------------------------------------

def plot_sphere(args, allData):

    altGoal = args.alt
    sliceData = utils.data_slice(allData, targetAlt = altGoal)
    allSlices = sliceData['slices']
    realAlt = sliceData['realAlt']
    lons2d = sliceData['lons2d']
    lats2d = sliceData['lats2d']

    varName = allData['longname'][0]
    sVarNum = allData['shortname'][0] + '_'
    sAltNum = 'alt%04d_' % int(realAlt)

    if (len(args.backdir) > 0):
        sAltNum = sAltNum + 'diff_'
        if (args.percent):
            varName = varName + ' (Per. Diff.)'
        else:
            varName = varName + ' (Diff.)'

    allTimes = allData['times']

    # get min and max values, plus color table:
    dataMinMax = plotutils.get_min_max_data(allSlices, None, \
                                            color = 'red', \
                                            minVal = args.mini, \
                                            maxVal = args.maxi)

    sFilePre = sVarNum + sAltNum
    sTitleAdd = '; Alt: %.0f km' % realAlt
    plotters.plot_slices_wpolar(allSlices,
                                allTimes,
                                lons2d,
                                lats2d,
                                dataMinMax,
                                xLabel = 'Longitude (deg)',
                                yLabel = 'Latitude (deg)',
                                varName = varName,
                                titleAddOn = sTitleAdd,
                                filenamePrefix = sFilePre,
                                xLimits = [args.lonmin, args.lonmax],
                                yLimits = [args.latmin, args.latmax])   

    return

# ----------------------------------------------------------------------------
# This assumes you have blocks, like a cubesphere grid
#-----------------------------------------------------------------------------

def plot_cubesphere(args, allData):

    altGoal = args.alt
    sliceData = utils.data_slice(allData, targetAlt = altGoal)
    allSlices = sliceData['slices']
    realAlt = sliceData['realAlt']
    lons3d = sliceData['lons2d']
    lats3d = sliceData['lats2d']

    varName = allData['longname'][0]
    sVarNum = allData['shortname'][0] + '_'
    sAltNum = 'alt%04d_' % int(realAlt)

    if (len(args.backdir) > 0):
        sAltNum = sAltNum + 'diff_'
        if (args.percent):
            varName = varName + ' (Per. Diff.)'
        else:
            varName = varName + ' (Diff.)'

    allTimes = allData['times']

    # get min and max values, plus color table:
    dataMinMax = plotutils.get_min_max_data(allSlices, None, \
                                            color = 'red', \
                                            minVal = args.mini, \
                                            maxVal = args.maxi)
    
    #lons3d = allData['lons'][:, :, :, iAlt]
    #lats3d = allData['lats'][:, :, :, iAlt]
    
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
                                           xLimits = [args.lonmin, args.lonmax],
                                           yLimits = [args.latmin, args.latmax],
                                           forcePolar = True)

# Needed to run main script as the default executable from the command line
if __name__ == '__main__':

    # Get the input arguments
    args = get_args()
    filelist = args.filelist
    varToPlot = args.var
    isVerbose = args.v

    if (args.list):
        util.list_file_info(filelist)
        exit()

    allData = util.read_all_files(filelist, varToPlot, \
                                  verbose = isVerbose, \
                                  iStart = args.iStart, \
                                  iEnd = args.iEnd, \
                                  iStep = args.iStep, \
                                  start = args.start, \
                                  stop = args.stop, \
                                  time = args.time)

    if (not allData):
        util.list_file_info(filelist)
        exit()

    if (len(args.backdir) > 0):
        backfiles = util.find_files_in_different_directory(filelist, \
                                                           dir = args.backdir)
        allBackground = util.read_all_files(backfiles, \
                                            varToPlot, \
                                            verbose = isVerbose,
                                            iStart = args.iStart, \
                                            iEnd = args.iEnd, \
                                            iStep = args.iStep, \
                                            start = args.start, \
                                            stop = args.stop, \
                                            time = args.time)
        allData = utils.subtract_all_slices(allData, \
                                            allBackground, \
                                            percent = args.percent)

    altGoal = args.alt
    
    if (allData['nblocks'] == 0):
        plot_sphere(args, allData)
        
    if (allData['nblocks'] > 0):
        plot_cubesphere(args, allData)
        
