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
# Plot polar region
# ----------------------------------------------------------------------------

def plot_polar_region(fig, axplot, lon_pos2d, lat_pos2d, values2d,
                      utime, whichPole,
                      miniPole, maxiPole, cmap,
                      title = '', cbar_label = ''):

    maxR = 50.0
    
    if (whichPole =='North'):
        ylabels = [r'80$^\circ$', r'70$^\circ$', r'60$^\circ$',
                   r'50$^\circ$']
        fac = 1.0
    else:
        ylabels = [r'-80$^\circ$', r'-70$^\circ$', r'-60$^\circ$',
                   r'-50$^\circ$']
        fac = -1.0
        
    # Find rotation for convertion to local time
    shift = time_conversion.calc_time_shift(utime)

    xlabels = []
    xlabelpos = []

    ylabelpos = [10.0, 20.0, 30.0, 40.0]
    xticks = np.arange(0, 2 * np.pi, np.pi / 2.0)
    yticks = np.arange(10, 50, 10)

    rp2d = 90.0 - fac * lat_pos2d
    tp2d = np.radians(lon_pos2d + shift - 90.0)

    x2d = rp2d * np.cos(tp2d)
    y2d = rp2d * np.sin(tp2d)

    x2de = utils.move_centers_to_corners(x2d)
    y2de = utils.move_centers_to_corners(y2d)
    
    print(x2de)
    axplot.grid(False)
    conn = axplot.pcolor(x2de, y2de, values2d,
                             shading = 'auto',
                             vmin = miniPole, vmax = maxiPole,
                             cmap = cmap)
    axplot.set_aspect(1.0)

    def rt_to_xy(r, t):
        return r * np.cos(t), r * np.sin(t)
    
    t = np.radians(np.arange(0, 361, 1))
    for r in np.arange(10, maxR + 10, 10):
        x, y = rt_to_xy(r, t)
        axplot.plot(x, y, linestyle = ':', color = 'k', linewidth = 0.5)
        x, y = rt_to_xy( r, 3 * np.pi/4)
        axplot.text(x, y, '%d' % int(90.0 - r), fontsize = 10)

    axplot.plot([-maxR, maxR], [0.0, 0.0], \
                linestyle = ':', color = 'k', linewidth = 0.5)
    axplot.plot([0.0, 0.0], [-maxR, maxR], \
                linestyle = ':', color = 'k', linewidth = 0.5)
    axplot.set_axis_off()

    x, y = rt_to_xy( maxR * 1.5, 132.0/180. * np.pi)
    axplot.text(x, y, title,
                verticalalignment = 'top',
                horizontalalignment = 'left',
                fontsize = 14)
    x, y = rt_to_xy( maxR + 1, 55/180. * np.pi)
    axplot.text(x, y, cbar_label,
                verticalalignment = 'bottom',
                horizontalalignment = 'left',
                fontsize = 14)
    x, y = rt_to_xy( maxR + 1, -np.pi/2)
    axplot.text(x, y, '00 LT',
             verticalalignment='top', horizontalalignment='center')
    x, y = rt_to_xy( maxR + 1, np.pi/2)
    axplot.text(x, y, '12 LT',
             verticalalignment='bottom', horizontalalignment='center')
    x, y = rt_to_xy( maxR + 2, -np.pi)
    axplot.text(x, y, '18 LT',
             verticalalignment='center', horizontalalignment='center',
             rotation = 90)
    x, y = rt_to_xy(maxR, 3*np.pi/4)
    axplot.text(x, y, whichPole,
             verticalalignment='bottom',
             horizontalalignment='center', rotation = 45)
    axplot.set_xlim([-maxR, maxR])
    axplot.set_ylim([-maxR, maxR])

    cbar = fig.colorbar(conn, ax = axplot, shrink=0.5, pad=0.005)

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
                    
                    plot_polar_region(fig, ax2, lon2d, lat2d, values2d,
                                      uTime, whichPole,
                                      dataMinMax['mini'],
                                      dataMinMax['maxi'],
                                      dataMinMax['cmap'],
                                      title = '', cbar_label = '')
                if (np.mean(lat2d) < -45.0):
                    whichPole = 'South'
                    cbar_label = ''
                    
                    plot_polar_region(fig, ax3, lon2d, lat2d, values2d,
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
