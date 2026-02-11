#!/usr/bin/env python3
""" This is a simple lat/lon plotter
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from pyitm.modeldata import utils
from pyitm.general import time_conversion

# ----------------------------------------------------------------------------
# Plot 1d variable vs time
#   allSlices: 2D - [time, ypositions]
#   allTimes: 1D - [times] 
#   yPos: 1D - [ypositions] 
#   dataMinMax: dictionary that contains info on min, max, and color map
#   varName: puts this string on the colorbar
#   titleAddOn: puts this string after the title
#   yLabel: string label for the y axis
#   filenamePrefix: string to add before time for filename (like varX_)
#   yLimits: limits for the y axis - (array of 2 elements)
#   dpi: dots per inch for file
# ----------------------------------------------------------------------------

def plot_vs_time(allSlices,
                 allTimes,
                 yPosCenters,
                 dataMinMax,
                 varName = '',
                 titleAddOn = '',
                 yLabel = '',
                 filenamePrefix = '',
                 yLimits = [0, 0],
                 dpi = 120):

    fig = plt.figure(figsize=(10, 5.5), dpi = dpi)
    ax = fig.add_axes([0.07, 0.08, 0.97, 0.88])

    value2d = allSlices.transpose()
    con = ax.pcolor(allTimes, yPosCenters, value2d, \
                    cmap = dataMinMax['cmap'], \
                    vmin = dataMinMax['mini'], \
                    vmax = dataMinMax['maxi'])

    sTime = allTimes[0].strftime("%d %b %Y %H:%M:%S UT")
    eTime = allTimes[-1].strftime("%d %b %Y %H:%M:%S UT")

    ax.set_xlabel('Time from ' + sTime + ' to ' + eTime)
    ax.set_title(varName)
    ax.set_ylabel(yLabel)
    ax.set_ylim(yLimits[0], yLimits[1])

    cbar = fig.colorbar(con, ax = ax, shrink = 0.75, pad = 0.02)
    cbar.set_label(varName, rotation=90)

    outFile = filenamePrefix + 'alt_vs_time.png'
    print(" ==> Writing file : ", outFile)
    fig.savefig(outFile, dpi = dpi)
    plt.close(fig)

    return

    
