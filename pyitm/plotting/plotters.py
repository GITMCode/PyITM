#!/usr/bin/env python3
""" This is a simple lat/lon plotter
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------------
# Plot a series of slices
#   allSlices: 3D - [time, xpositions, ypositions]
#   xPosEdges: 1D - [xpositions] (edges of cells, not centers)
#   yPosEdges: 1D - [ypositions] (edges of cells, not centers)
#   dataMinMax: dictionary that contains info on min, max, and color map
#   varName: puts this string on the colorbar
#   titleAddOn: puts this string after the time on the title
#   xLabel: string label for the x axis
#   yLabel: string label for the y axis
#   filenamePrefix: string to add before time for filename (like varX_)
#   xLimits: limits for the x axis - (array of 2 elements)
#   yLimits: limits for the y axis - (array of 2 elements)
#   dpi: dots per inch for file
# ----------------------------------------------------------------------------

def plot_series_of_slices(allSlices,
                          allTimes,
                          xPosEdges,
                          yPosEdges,
                          dataMinMax,
                          varName = '',
                          titleAddOn = '',
                          xLabel = '',
                          yLabel = '',
                          filenamePrefix = '',
                          xLimits = [0, 0],
                          yLimits = [0, 0],
                          dpi = 120):

    for iTime, uTime in enumerate(allTimes):

        fig = plt.figure(figsize=(10, 5.5), dpi = dpi)
        ax = fig.add_axes([0.07, 0.06, 0.97, 0.9])

        title = uTime.strftime("%d %b %Y %H:%M:%S UT") + titleAddOn
        value2d = allSlices[iTime, :, :].transpose()
        con = ax.pcolormesh(xPosEdges, yPosEdges, value2d, \
                            cmap = dataMinMax['cmap'], \
                            vmin = dataMinMax['mini'], \
                            vmax = dataMinMax['maxi'])
        if (xLimits[1] > xLimits[0]):
            ax.set_xlim(xLimits)
        if (yLimits[1] > yLimits[0]):
            ax.set_ylim(yLimits)
        ax.set_aspect(1.0)
        ax.set_title(title)
        ax.set_xlabel(xLabel)
        ax.set_ylabel(yLabel)

        cbar = fig.colorbar(con, ax = ax, shrink = 0.75, pad = 0.02)
        cbar.set_label(varName, rotation=90)

        sTimeOut = uTime.strftime('%y%m%d_%H%M%S')
        outFile = filenamePrefix + sTimeOut + '.png'
        print(" ==> Writing file : ", outFile)
        fig.savefig(outFile, dpi = dpi)
        plt.close(fig)
    
