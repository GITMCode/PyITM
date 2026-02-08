#!/usr/bin/env python3
""" This is a simple lat/lon plotter
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse

from pyitm.fileio import util
from pyitm.modeldata import utils
from pyitm.plotting import plotutils
from pyitm.plotting import plotters, contour_wtime

# ----------------------------------------------------------------------------
# Get arguments as inputs into the code
#-----------------------------------------------------------------------------

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Plot Aether / GITM model results')
    
    # select latitude to plot:
    parser.add_argument('-lat', metavar = 'lat',
                        default = -1e32, type = float, \
                        help = 'latitude to plot in deg (closest)') 
    # select longitude to plot:
    parser.add_argument('-lon', metavar = 'lon',
                        default = -1e32, type = float, \
                        help = 'longitude to plot in deg (closest)') 

    # variable to plot as a number
    parser.add_argument('-var',  \
                        default = '3', \
                        help = 'variable to plot (number or variable name)')

    # User can set max and min of the colorbar:
    parser.add_argument('-mini',  default = 1e32, type = float, \
                        help = 'manually set the minimum value for the plots')
    parser.add_argument('-maxi',  default = -1e32, type = float, \
                        help = 'manually set the maxiumum value for the plots')
    
    parser.add_argument('-altmin',  default = -1e32, type = float, \
                        help = 'manually set the minimum altitude for the plots')
    parser.add_argument('-altmax',  default = -1e32, type = float, \
                        help = 'manually set the maxiumum altitude for the plots')

    # directory to use as a background, so you can subtract one run from another
    parser.add_argument('-backdir',  \
                        default = '', \
                        help = 'background directory to subtract off')
    parser.add_argument('-percent',  \
                        action='store_true', default = False, \
                        help = 'plot percent different (if backdir specified)')

    parser.add_argument('-log',  \
                        action='store_true', default = False, \
                        help = 'plot log of variable')
    
    parser.add_argument('-list',  \
                        action='store_true', default = False, \
                        help = 'list variables in file')
    
    # Get the files to plot:
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')

    args = parser.parse_args()

    return args


# ----------------------------------------------------------------------------
# Should write an interpolator for lat / lon.
# ----------------------------------------------------------------------------



# ----------------------------------------------------------------------------
# Needed to run main script as the default executable from the command line
# ----------------------------------------------------------------------------

if __name__ == '__main__':

    # Get the input arguments
    args = get_args()
    filelist = args.filelist
    varToPlot = args.var

    if (args.list):
        util.list_file_info(filelist)
        exit()
    
    allData = util.read_all_files(filelist, varToPlot)

    if (not allData):
        util.list_file_info(filelist)
        exit()

    if (len(args.backdir) > 0):
        backfiles = util.find_files_in_different_directory(filelist, \
                                                           dir = args.backdir)
        allBackground = util.read_all_files(backfiles, \
                                            varToPlot, \
                                            verbose = True)
        allBackSlices = allBackground['data'][:, 0, 0, :]
        useBack = True
    else:
        useBack = False

    if ((allData['nlons'] == 1) and (allData['nlats'] == 1)):
        alts = allData['alts'][0,0,:]
        altmin = args.altmin
        altmax = args.altmax
        if (altmin == -1e32):
            altmin = np.min(alts)
        if (altmax == -1e32):
            altmax = np.max(alts)

        varName = allData['longname'][0]
        sVarNum = allData['shortname'][0] + '_'

        allSlices = allData['data'][:, 0, 0, :]

        if (useBack):
            if (args.percent):
                allSlices = (allSlices - allBackSlices) / allBackSlices * 100.0
                varName = varName + '(% diff)'
                sVarNum = sVarNum + 'perdiff_'
            else:
                allSlices = (allSlices - allBackSlices)
                varName = varName + '(diff)'
                sVarNum = sVarNum + 'diff_'

        if (args.log):
            allSlices = np.log(allSlices)
            varName = 'log(' + varName + ')'
            sVarNum = 'log_' + sVarNum
        
        # get min and max values, plus color table:
        dataMinMax = plotutils.get_min_max_data(allSlices, alts, \
                                                yMin = altmin, \
                                                yMax = altmax, \
                                                color = 'red', \
                                                minVal = args.mini, \
                                                maxVal = args.maxi)
        
        contour_wtime.plot_vs_time(allSlices,
                                   allData['times'],
                                   alts,
                                   dataMinMax,
                                   varName = varName,
                                   titleAddOn = '',
                                   yLabel = 'Altitude (km)',
                                   filenamePrefix = sVarNum,
                                   yLimits = [altmin, altmax])
