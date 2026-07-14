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
        description = 'Make altitude plot, reading in CSV controller file')
    
    parser.add_argument('-xlabel', \
                        default = 'none', \
                        help = 'label for the x axis (default = none)')
    parser.add_argument('-title', \
                        default = 'none', \
                        help = 'title of plot (default = none)')

    # User can set max and min of the colorbar:
    parser.add_argument('-mini',  default = 1e32, type = float, \
                        help = 'manually set the minimum value for the plots')
    parser.add_argument('-maxi',  default = -1e32, type = float, \
                        help = 'manually set the maxiumum value for the plots')
    
    parser.add_argument('-altmin',  default = -1e32, type = float, \
                        help = 'manually set the minimum altitude for the plots')
    parser.add_argument('-altmax',  default = -1e32, type = float, \
                        help = 'manually set the maxiumum altitude for the plots')

    parser.add_argument('-dpi',  default = 120, type = int, \
                        help = 'dots per inch for plot file')
    
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
    parser.add_argument('csvfile', nargs = 1, \
                        help = 'csv controller file')

    args = parser.parse_args()

    return args


# ----------------------------------------------------------------------------
# Needed to run main script as the default executable from the command line
# ----------------------------------------------------------------------------

if __name__ == '__main__':

    # Get the input arguments
    args = get_args()
    plt.rcParams.update({'font.size': 14})
    nGCs = 2
    
    controller = controlfile.read_logfile_styles(args.csvfile[0])
    nRows = len(controller['files'])
    
    if (nRows == 0):
        print('Controller file provided, but something went wrong! Stopping!')
        exit()

    filelist = controller['files']
    if (args.list):
        util.list_file_info(filelist)
        exit()
        
    # we are going to make one single plot, reading in each file one by one and plotting it:

    fig = plt.figure(figsize=(10, 6), dpi = args.dpi)
    ax = fig.add_axes([0.08, 0.1, 0.9, 0.85])
     
    vars = controller['vars']

    for iRow in range(nRows):
        allData = util.read_all_files(filelist[iRow], vars[iRow])

        if (not allData):
            util.list_file_info(filelist)
            exit()

        alts = allData['alts'][0,0,nGCs:-nGCs]
        varNames = allData['longname']

        vals = allData['data'][0, 0, 0, nGCs:-nGCs]
        ax.plot(vals, alts, label = controller['labels'][iRow], \
                color = controller['colors'][iRow], \
                linestyle = controller['styles'][iRow], \
                linewidth = controller['widths'][iRow])
    
    ax.axvline(0.0, linestyle = 'dotted', color = 'grey')
    ax.set_ylabel('Altitude (km)')
    yLimits = [np.min(alts), np.max(alts)]
    if (args.altmin > -1e32):
        yLimits[0] = args.altmin 
    if (args.altmax > -1e32):
        yLimits[1] = args.altmax 
    ax.set_ylim(yLimits[0], yLimits[1])
    if (args.xlabel != 'none'):
        ax.set_xlabel(args.xlabel)

    if ((args.maxi > -1e32) and (args.mini < 1e32)):
        ax.set_xlim(args.mini, args.maxi)

    if (args.log):
        ax.set_xscale('log')

    if (args.title != 'none'):
        ax.set_title(args.title)
        
    ax.legend()
    outfile = args.plotfile
    print(" ==> Writing file : ", outfile)
    fig.savefig(outfile, dpi = args.dpi)
    plt.close(fig)

