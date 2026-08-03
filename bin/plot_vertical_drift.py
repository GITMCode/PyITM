#!/usr/bin/env python3
""" Extract vertical ion drift and plot as a function of local time
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse

from pyitm.fileio import util, logfile
from pyitm.modeldata import satellite
from pyitm.plotting import axes

# ----------------------------------------------------------------------------
# Get arguments as inputs into the code
#-----------------------------------------------------------------------------

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Extract vertical ion drifts and plot vs local time')
    
    # select latitude to plot:
    parser.add_argument('-lat', metavar = 'lat',
                        default = -12.0, type = float, \
                        help = 'latitude to plot in deg (closest)') 
    # select longitude to plot:
    parser.add_argument('-lon', metavar = 'lon',
                        default = 283.0, type = float, \
                        help = 'longitude to plot in deg (closest)') 
    parser.add_argument('-alt',  default = 175, type = float, \
                        help = 'altitude to plot in km')

    # User can set max and min of the plot:
    parser.add_argument('-mini',  default = 1e32, type = float, \
                        help = 'manually set the minimum value for the plots')
    parser.add_argument('-maxi',  default = -1e32, type = float, \
                        help = 'manually set the maxiumum value for the plots')
    
    parser.add_argument('-plotfile',
                        help = 'output file for plot',
                        default = 'vertical_ion_drift.png')

    parser.add_argument('-var',
                        help = 'Variable that is the vertical ion drift',
                        default = 'Viv')

    parser.add_argument('-v',  \
                        action='store_true', default = False, \
                        help = 'set verbose to true')
    
    # Get the files to plot:
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------------
# Needed to run main script as the default executable from the command line
# ----------------------------------------------------------------------------

if __name__ == '__main__':

    # Get the input arguments
    args = get_args()
    verbose = args.v
    plt.rcParams.update({'font.size': 14})

    modelData = util.read_all_files(args.filelist, \
                                   varsToRead = args.var, \
                                   verbose = verbose)
    locInfo = {'times': [],
               'lons': [],
               'lats': [],
               'alts': []}
    
    localtimes = []
    for time in modelData['times']:
        locInfo['times'].append(time)
        locInfo['lons'].append(args.lon)
        locInfo['lats'].append(args.lat)
        locInfo['alts'].append(args.alt)
        ut = time.hour + time.minute/60.0 + time.second/60.0
        localtimes.append((ut + args.lon) % 24.0)

    modelViv = satellite.extract_1d(locInfo, modelData, extrapolate=False, \
                                    interpVar=[0], \
                                    skipTimeCheck=True, verbose=verbose)

    dpi = 150
    fig = plt.figure(figsize=(10, 10), dpi = dpi)
    nPlots = 1
    yBuffer = 1e-6
    yTop = 0.05
    yBottom = 0.075
    xRight = 0.05
    ax = axes.get_axes_one_column(fig,
                                  nPlots,
                                  yBottom,
                                  yTop,
                                  yBuffer,
                                  xRight = xRight)

    nTimes = len(localtimes)
    localtimesNew = []
    dataNew = []
    timeNew = []
    for i in range(nTimes-1):
        localtimesNew.append(localtimes[i])
        timeNew.append( (locInfo['times'][i] - \
                         locInfo['times'][0]).total_seconds()/3600.0)
        dataNew.append(modelViv['model_'+args.var][i])
        if (np.abs(localtimes[i+1] - localtimes[i]) > 18.0):
            localtimesNew.append(None)
            dataNew.append(None)
            timeNew.append(None)
                    
    localtimesNew.append(localtimes[-1])
    dataNew.append(modelViv['model_'+args.var][-1])
    timeNew.append( (locInfo['times'][-1] - \
                     locInfo['times'][0]).total_seconds()/3600.0)
    ax[0].plot(localtimesNew, \
               dataNew,\
               color = 'b', label = 'Vertical Ion Drift')

    con1 = ax[0].scatter(localtimesNew, \
                         dataNew,\
                         c = timeNew) 
    cbar1 = fig.colorbar(con1, ax = ax[0], shrink = 0.5, pad = 0.02)
    sTime = locInfo['times'][0].strftime("%d %b %Y %H:%M UT")
    eTime = locInfo['times'][-1].strftime("%d %b %Y %H:%M UT")
    
    cbar1.set_label(sTime + ' to ' + eTime + ' (Hours)', rotation=90)
    
    ax[0].set_ylabel('Vertical Ion Drift (m/s)')
    ax[0].set_xlabel('Local Time (hrs)')
    ax[0].axhline(0.0, color = 'grey', linestyle = '--')
    ax[0].axvline(6.0, color = 'orange', linestyle = '--')
    ax[0].axvline(18.0, color = 'orange', linestyle = '--')
    outFile = args.plotfile
    print(" ==> Writing file : ", outFile)
    fig.savefig(outFile, dpi = dpi)
    plt.close(fig)
    

    
