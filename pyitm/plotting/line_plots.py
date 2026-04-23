#!/usr/bin/env python3
""" This is a simple lat/lon plotter
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from pyitm.plotting import axes

# ----------------------------------------------------------------------------
# Plot 1d variable vs altitue
#   allSlices: 2D - [time, ypositions]
#   yPos: 1D - [ypositions] 
#   dataMinMax: dictionary that contains info on min, max, and color map
#   varName: puts this string on the colorbar
#   titleAddOn: puts this string after the title
#   yLabel: string label for the y axis
#   filenamePrefix: string to add before time for filename (like varX_)
#   yLimits: limits for the y axis - (array of 2 elements)
#   dpi: dots per inch for file
# ----------------------------------------------------------------------------

def lineplot_alt(allSlices,
                 allTimes,
                 yPosCenters,
                 dataMinMax,
                 varNames = [''],
                 timeLabels = None,
                 titleAddOn = '',
                 xLabel = '',
                 outfile = 'alt_plot.png',
                 linestyles = ['solid', 'dashed', 'dotted', 'dash-dotted'],
                 colors = ['k','b','r','c','m','darkgreen'],
                 yLimits = [0, 0],
                 nGCs = 2, 
                 dpi = 120):

    fig = plt.figure(figsize=(10, 6), dpi = dpi)
    ax = fig.add_axes([0.08, 0.1, 0.9, 0.85])
    
    nTimes, nVars, nAlts = np.shape(allSlices)
    print(nTimes, nVars, nAlts)
    alts = yPosCenters[nGCs:-nGCs]
    for iVar in range(nVars):
        print(iVar, nVars)
        for iTime in range(nTimes):
            varName = varNames[iVar]
            if (nTimes > 1):
                if (timeLabels == None):
                    t = allTimes[iTime].strftime('%H%M UT')
                    varName = varName + ' ' + t
                else:
                    varName = varName + ' ' + timeLabels[iTime]
            linestyle = linestyles[iTime % len(linestyles)]
            color = colors[iVar % len(colors)]
            ax.plot(allSlices[iTime, iVar, nGCs:-nGCs], \
                    alts, \
                    label = varName,
                    linestyle = linestyle, \
                    color = color)

    ax.axvline(0.0, linestyle = 'dotted', color = 'grey')
    ax.set_ylabel('Altitude (km)')
    if (yLimits[0] != yLimits[1]):
        ax.set_ylim(yLimits[0], yLimits[1])
    if (xLabel != 'none'):
        ax.set_xlabel(xLabel)
    ax.legend()
    print(" ==> Writing file : ", outfile)
    fig.savefig(outfile, dpi = dpi)
    plt.close(fig)

    return
    

# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------

def lineplot_data(data, outFile=None, fig=None, ax=None, vars = None,
                  logScale = False,
                  linewidth = 1.0, color = None, linestyle = None,
                  xStart = None, xEnd = None,
                  yMin = None, yMax = None,
                  medianTime = -1,
                  dpi = 120,
                  ylabel = None, label = None, nLabels = None, title = None):
    """
    Make lineplot out of data

    if outFile is None, the figure & axes are returned, otherwise figure is saved

    if fig is None, figure is created. Otherwise, lines are drawn on existing plot.

    if vars is None, it will plot out all of the variables, else only listed vars
    
    """

    if (vars == None):
        vars = data['vars']
    nVars = len(vars)

    xSize = 10.0
    ySize = 6.0 * nVars

    if fig is None:
        fig = plt.figure(figsize=(xSize, ySize), dpi = dpi)

        yBot = 0.11
        yTop = 0.05
        yBuf = 0.05
        ax = axes.get_axes_one_column(fig,
                                      nVars,
                                      yBot,
                                      yTop,
                                      yBuf)

    times = data['times']
    if (xStart == None):
        xStart = times[0]
    if (xEnd == None):
        xEnd = times[-1]
    for iPlot, var in enumerate(vars):
        vals = data[var]
        if (medianTime > 0):
            iT = 0
            totalVals = 0.0
            dt = (times[iT] - xStart).total_seconds()/3600
            print(dt)
            while (dt < 0):
                iT = iT + 1
                dt = (times[iT] - xStart).total_seconds()/3600
            nPts = 0
            while ((dt < medianTime) and (dt >= 0)):
                totalVals = totalVals + vals[iT]
                iT = iT + 1
                dt = (times[iT] - xStart).total_seconds()/3600
                nPts = nPts + 1
            meanVal = totalVals / nPts
            print(' -> subtracting : ', totalVals / nPts, iT, nPts, np.shape(vals), xStart)
            vals = (vals - meanVal) / meanVal * 100.0
        ax[iPlot].plot(times, vals, \
                       linewidth = linewidth, \
                       linestyle = linestyle, \
                       label = label, \
                       color = color)
        if (ylabel == None):
            ylabel = var
        ax[iPlot].set_ylabel(ylabel)
        ax[iPlot].set_xlim(xStart, xEnd)
        if (title != None):
            ax[iPlot].set_title(title)
        if ((yMin != None) and (yMax != None)):
            ax[iPlot].set_ylim(yMin, yMax)

    sTime = xStart.strftime("%d %b %Y %H:%M:%S UT")
    eTime = xEnd.strftime("%d %b %Y %H:%M:%S UT")
    ax[-1].set_xlabel('Time from ' + sTime + ' to ' + eTime)
    if (label != None):
        if (nLabels != None):
            ax[-1].legend(ncol = int(np.ceil(nLabels/5)))
        else:
            ax[-1].legend()

    if ('title' in data.keys()):
        ax[0].set_title(data['title'])

    if (logScale):
        ax[0].set_yscale('log')
        
    if outFile:
        print(" ==> Writing file : ", outFile)
        fig.savefig(outFile, dpi = dpi)
        plt.close(fig)
    else:
        return fig, ax
