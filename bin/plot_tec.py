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
from pyitm.plotting import plotters

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Plot Aether / GITM model results')
    
    # average low files together
    parser.add_argument('-low', default = -1, type = int, \
                        help = 'low pass filter (average low files together)') 
    # average high files together and subtract them for each file
    parser.add_argument('-high', default = -1, type = int, \
                        help = 'high pass filter (ave files together and subtract)') 

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
    varToPlot = [34]

    allData3D = gitmio.read_gitm_all_files(filelist, varToPlot)

    lons1d = allData3D['lons'][:, 0, 0]
    lats1d = allData3D['lats'][0, :, 0]
    lonsEdge = utils.move_centers_to_edges(lons1d)
    latsEdge = utils.move_centers_to_edges(lats1d)

    allSlicesRaw = utils.calc_tec(allData3D)
    sVarName = 'TEC'
    sFilePre = 'varTEC_'
    allTimes = allData3D['times']
    dt = (allTimes[1] - allTimes[0]).total_seconds()

    sTitleAdd = ' TEC '
    if (args.low > 0):
        sVarName = 'TEC (low pass)'
        sFilePre = 'varTECl_'
        dtFilter = dt * args.low/60.0
        sTitleAdd = sTitleAdd + '(%.0fm Low Pass Filtered)' % dtFilter
    if (args.high > 0):
        sVarName = 'TEC (high pass)'
        sFilePre = 'varTECh_'
        dtFilter = dt * args.high/60.0
        sTitleAdd = sTitleAdd + '(%.0fm High Pass Filtered)' % dtFilter

    allSlices = utils.filter_slices(allSlicesRaw, \
                                    iLowPass = args.low, iHighPass = args.high)        

    # get min and max values, plus color table:
    dataMinMax = plotutils.get_min_max_data(allSlices, lats1d, \
                                            yMin = args.latmin, yMax = args.latmax, \
                                            color = 'red', \
                                            minVal = args.mini, maxVal = args.maxi)

    plotters.plot_series_of_slices(allSlices,
                                   allTimes,
                                   lonsEdge,
                                   latsEdge,
                                   dataMinMax,
                                   xLabel = 'Longitude (deg)',
                                   yLabel = 'Latitude (deg)',
                                   varName = sVarName,
                                   titleAddOn = sTitleAdd,
                                   filenamePrefix = sFilePre,
                                   xLimits = [args.lonmin, args.lonmax],
                                   yLimits = [args.latmin, args.latmax])                   

    
