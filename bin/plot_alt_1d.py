#!/usr/bin/env python3
""" This is a simple lat/lon plotter
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse

from pyitm.fileio import util, controlfile
from pyitm.modeldata import utils
from pyitm.plotting import plotutils, plotters, line_plots

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

    parser.add_argument('-csv',
                        help = 'file to use to describe files, variables and styles',
                        default = 'none')

    # variable to plot as a number
    parser.add_argument('-var',  \
                        default = '3', \
                        help = 'variable to plot (number or variable name)')
    # variables to plot as numbers (can choose multiple variables or one)
    parser.add_argument('-vars',  \
                        default = '', \
                        help = 'variables to plot (numbers or names)')
    parser.add_argument('-xlabel', \
                        default = 'none', \
                        help = 'label for the x axis (default = none)')

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
    
    parser.add_argument('-plotfile',
                        help = 'output file for plot',
                        default = 'alt_plot.png')

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
    plt.rcParams.update({'font.size': 14})

    useController = False
    if (args.csv != 'none'):
        controller = controlfile.read_logfile_styles(args.csv)
        if (len(controller['files']) > 0):
            useController = True
        else:
            print('Controller file provided, but something went wrong! Stopping!')
            exit()

    if (useController):
        filelist = controller['files']
        varsToPlot = controller['vars']
        # stuff is returned in a weird way from this.
        # find all of the unique variable names
        print(filelist)
        print(varsToPlot)
        uniqueVars = []
        uniqueFiles = []
        nRows = len(varsToPlot)
        for iRow in range(nRows):
            if (not (varsToPlot[iRow][0] in uniqueVars)):
                uniqueVars.append(varsToPlot[iRow][0])
            if (not (filelist[iRow] in uniqueFiles)):
                uniqueFiles.append(filelist[iRow])
        print(uniqueFiles)
        print(uniqueVars)
        nVars = len(uniqueVars)
        nFiles = len(uniqueFiles)
        if (nRows != nVars * nFiles):
            print('Controller file mismatch!')
            print('  -> This code only works if each nFiles * nVars = number of rows in controller file')
            print('  -> Each file to read needs to have all of the variables listed')
            exit()
        colors = controller['colors'][0]
        linestyles = controller['styles'][0]
        linewidths = controller['widths'][0]
        labels = controller['labels']
        filelist = uniqueFiles
        varsToPlot = uniqueVars
        print(colors)
    else:
        filelist = args.filelist
        varsToPlot = args.var
        labels = None
        if (len(args.vars) > 0):
            if (',' in args.vars):
                varsToPlot = args.vars.split(',')
            else:
                varsToPlot = [args.vars]
    
    if (args.list):
        util.list_file_info(filelist)
        exit()

    print('reading files...')
    allData = util.read_all_files(filelist, varsToPlot)

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

        varNames = allData['longname']

        s = np.shape(allData['data'])
        print('Shape of allData : ', s)
        if (len(s) <= 4):
            # only one variable given
            allSlices = np.zeros((s[0], 1, s[-1]))
            allSlices[:,0,:] = allData['data'][:, 0, 0, :]
        if (len(s) == 5):
            # multiple variables given
            allSlices = allData['data'][:, :, 0, 0, :]

        if (useBack):
            if (args.percent):
                allSlices = (allSlices - allBackSlices) / allBackSlices * 100.0
                varName = varName + '(% diff)'
            else:
                allSlices = (allSlices - allBackSlices)
                varName = varName + '(diff)'

        if (args.log):
            allSlices = np.log(allSlices)
            varName = 'log(' + varName + ')'
        
        # get min and max values, plus color table:
        dataMinMax = plotutils.get_min_max_data(allSlices, alts, \
                                                yMin = altmin, \
                                                yMax = altmax, \
                                                minVal = args.mini, \
                                                maxVal = args.maxi)
        
        line_plots.lineplot_alt(allSlices,
                                allData['times'],
                                alts,
                                dataMinMax,
                                varNames = varNames,
                                timeLabels = labels,
                                titleAddOn = '',
                                xLabel = args.xlabel,
                                outfile = args.plotfile,
                                yLimits = [altmin, altmax])

