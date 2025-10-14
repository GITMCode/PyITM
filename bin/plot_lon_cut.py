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

    # select altitude to plot:
    parser.add_argument('-ilon', metavar = 'ilon', default = -1, type = int, \
                        help = 'longitude index (> -1 to use this instead of lon)') 

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

    parser.add_argument('-log',  \
                        action='store_true', default = False, \
                        help = 'plot the log10() of the data')
    
    parser.add_argument('-list',  \
                        action='store_true', default = False, \
                        help = 'list variables in file')

    # directory to use as a background, so you can subtract one run from another
    parser.add_argument('-backdir',  \
                        default = '', \
                        help = 'background directory to subtract off')
    parser.add_argument('-percent',  \
                        action='store_true', default = False, \
                        help = 'plot percent different (if backdir specified)')

    # Get the files to plot:
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')

    args = parser.parse_args()

    return args

# -----------------------------------------------------------------------------
# Plot, assuming no blocks, just a nice sphere
# -----------------------------------------------------------------------------

def plot_lon_cut_noblocks(args, allData):

    if (args.ilon > -1):
        iLon = args.ilon
        realLon = 0.0
    else:
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
        sLonNum = 'lonZave_'
    else:
        sLonNum = 'lon%04d_' % int(realLon)

    if (len(args.backdir) > 0):
        sLonNum = sLonNum + 'diff_'
        if (args.percent):
            varName = varName + ' (Per. Diff.)'
        else:
            varName = varName + ' (Diff.)'

    if (args.log):
        allSlices = np.log10(allSlices)
        varName = 'log10(' + varName + ')'
        sLonNum = sLonNum + 'log_'

    allTimes = allData['times']

    # get min and max values, plus color table:
    dataMinMax = plotutils.get_min_max_data(allSlices, alts1d, \
                     yMin = args.altmin, yMax = args.altmax, \
                     color = 'red', \
                     minVal = args.mini, maxVal = args.maxi,
                     isLog = args.log)

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
        
    return

# -----------------------------------------------------------------------------
# plot with a bunch of different blocks
# -----------------------------------------------------------------------------

def plot_lon_cut_wblocks(args, allData):

    realLon = args.lon
    nBlocks = allData['nblocks']
    
    if (args.ilon > -1):
        iLon = args.ilon
        realLon = iLon
    else:
        if (args.lon >= 0):
            lonGoal = args.lon
            lons1d = allData['lons'][:, 0, 0]
            diff = np.abs(lons1d - lonGoal)
            iLon = np.argmin(diff)
            realLon = lons1d[iLon]
        else:
            realLon = -1.0
            iLon = -1

    allSlices = utils.data_slice(allData, iLon = iLon)

    varName = allData['longname'][0]
    sVarNum = allData['shortname'][0] + '_'
        
    if (args.lon < 0):
        sLonNum = 'lonZave_'
    else:
        sLonNum = 'lon%04d_' % int(realLon)

    if (args.log):
        allSlices = np.log10(allSlices)
        varName = 'log10(' + varName + ')'
        sLonNum = sLonNum + 'log_'

    if (len(args.backdir) > 0):
        sLonNum = sLonNum + 'diff_'

    allTimes = allData['times']

    # get min and max values, plus color table:
    dataMinMax = plotutils.get_min_max_data(allSlices, None, \
                    color = 'red', \
                    minVal = args.mini, maxVal = args.maxi)

    dpi = 120

    for iTime, uTime in enumerate(allTimes):

        fig = plt.figure(figsize=(10, 5.5), dpi = dpi)
        ax = fig.add_axes([0.07, 0.06, 0.97, 0.9])

        title = uTime.strftime("%d %b %Y %H:%M:%S UT")
        if (iLon > -1):
            title = title + '; Lon: %.0f deg' % realLon
        else:
            title = title + '; Zonal Average'

        for iBlock in range(nBlocks):
            alts2d = allData['alts'][iBlock, iLon, :, :]
            lats2d = allData['lats'][iBlock, iLon, :, :]
            #altsEdge = utils.move_centers_to_edges(alts1d)
            #latsEdge = utils.move_centers_to_edges(lats1d)
            value2d = allSlices[iTime, iBlock, :, :]
            con = ax.scatter(lats2d, alts2d, c = value2d, \
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
    
    return


# -----------------------------------------------------------------------------
# Needed to run main script as the default executable from the command line
# -----------------------------------------------------------------------------

if __name__ == '__main__':

    # Get the input arguments
    args = get_args()
    filelist = args.filelist
    varToPlot = args.var

    allData = util.read_all_files(filelist, varToPlot, verbose = True)

    if (len(args.backdir) > 0):
        backfiles = util.find_files_in_different_directory(filelist, dir = args.backdir)
        allBackground = util.read_all_files(backfiles, varToPlot, verbose = True)
        allData = utils.subtract_all_slices(allData, allBackground, percent = args.percent)

    if (allData['nblocks'] == 0):
        plot_lon_cut_noblocks(args, allData)
    else:
        plot_lon_cut_wblocks(args, allData)

