#!/usr/bin/env python3
""" This is a simple lat/lon plotter
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from pyitm.modeldata import utils
from pyitm.general import time_conversion
from pyitm.plotting import polar

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

    return

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

def plot_slices_wpolar(allSlices,
                       allTimes,
                       lons2d,
                       lats2d,
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

        fig = plt.figure(constrained_layout=False, figsize=(8.1, 7.5))
        ax = fig.add_axes([0.07, 0.06, 0.97, 0.48])
        # Top Left Graph Northern Hemisphere
        ax2 = fig.add_axes([0.06, 0.55, 0.425, 0.43])
        # Top Right Graph Southern Hemisphere
        ax3 = fig.add_axes([0.535, 0.55, 0.425, 0.43])

        title = uTime.strftime("%d %b %Y %H:%M:%S UT") + titleAddOn
        lon2dedges = utils.move_centers_to_corners(lons2d)
        lat2dedges = utils.move_centers_to_corners(lats2d)
        values2d = allSlices[iTime, :, :]
        con = ax.pcolormesh(lon2dedges, lat2dedges, values2d, \
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

        whichPole = 'North'
        cbar_label = ''
        lats1d = lats2d[0,:]
        mask = lats1d >= 40.0
        polar.plot_polar_region(fig, ax2,
                                lons2d[:,mask],
                                lats2d[:,mask],
                                values2d[:,mask],
                                uTime, whichPole,
                                dataMinMax['mini'],
                                dataMinMax['maxi'],
                                dataMinMax['cmap'],
                                title = '', cbar_label = '')

        whichPole = 'South'
        cbar_label = ''
        lats1d = lats2d[0,:]
        mask = lats1d <= -40.0
        polar.plot_polar_region(fig, ax3,
                                lons2d[:,mask],
                                lats2d[:,mask],
                                values2d[:,mask],
                                uTime, whichPole,
                                dataMinMax['mini'],
                                dataMinMax['maxi'],
                                dataMinMax['cmap'],
                                title = '', cbar_label = '')
        
        cbar = fig.colorbar(con, ax = ax, shrink = 0.75, pad = 0.02)
        cbar.set_label(varName, rotation=90)

        sTimeOut = uTime.strftime('%y%m%d_%H%M%S')
        outFile = filenamePrefix + sTimeOut + '.png'
        print(" ==> Writing file : ", outFile)
        fig.savefig(outFile, dpi = dpi)
        plt.close(fig)

    return

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

def plot_series_of_slices_wblocks(allSlices,
                                  allTimes,
                                  lonPos3d,
                                  latPos3d,
                                  dataMinMax,
                                  varName = '',
                                  titleAddOn = '',
                                  xLabel = '',
                                  yLabel = '',
                                  filenamePrefix = '',
                                  xLimits = [0, 0],
                                  yLimits = [0, 0],
                                  dpi = 120):

    nBlocks, nX, nY = np.shape(lonPos3d)
    
    for iTime, uTime in enumerate(allTimes):

        fig = plt.figure(constrained_layout=False, figsize=(8.1, 7.5))
        ax = fig.add_axes([0.07, 0.06, 0.97, 0.48])
        # Top Left Graph Northern Hemisphere
        ax2 = fig.add_axes([0.06, 0.55, 0.425, 0.43])
        # Top Right Graph Southern Hemisphere
        ax3 = fig.add_axes([0.535, 0.55, 0.425, 0.43])
        
        title = uTime.strftime("%d %b %Y %H:%M:%S UT") + titleAddOn

        for iBlock in range(nBlocks):
            lon2d = lonPos3d[iBlock, 2:-2, 2:-2]
            lat2d = latPos3d[iBlock, 2:-2, 2:-2]
            lon2dedges = utils.move_centers_to_corners(lon2d)
            lat2dedges = utils.move_centers_to_corners(lat2d)
            values2d = allSlices[iTime, iBlock, 2:-2, 2:-2]

            if (np.abs(np.mean(lat2d)) < 45.0):
                con = ax.pcolormesh(lon2dedges, lat2dedges, values2d, \
                                    cmap = dataMinMax['cmap'], \
                                    vmin = dataMinMax['mini'], \
                                    vmax = dataMinMax['maxi'])
            else:
                con = ax.scatter(lon2d, lat2d, c = values2d, \
                                 cmap = dataMinMax['cmap'], \
                                 vmin = dataMinMax['mini'], \
                                 vmax = dataMinMax['maxi'])

                if (np.mean(lat2d) > 45.0):
                    whichPole = 'North'
                    cbar_label = ''
                    polar.plot_polar_region(fig, ax2, lon2d, lat2d, values2d,
                                            uTime, whichPole,
                                            dataMinMax['mini'],
                                            dataMinMax['maxi'],
                                            dataMinMax['cmap'],
                                            title = '', cbar_label = '')
                if (np.mean(lat2d) < -45.0):
                    whichPole = 'South'
                    cbar_label = ''
                    polar.plot_polar_region(fig, ax3, lon2d, lat2d, values2d,
                                            uTime, whichPole,
                                            dataMinMax['mini'],
                                            dataMinMax['maxi'],
                                            dataMinMax['cmap'],
                                            title = '', cbar_label = '')
            
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
    
    return
