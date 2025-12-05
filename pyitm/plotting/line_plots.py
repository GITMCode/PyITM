#!/usr/bin/env python3
""" This is a simple lat/lon plotter
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from pyitm.plotting import axes

# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------

def lineplot_data(data, outFile=None, fig=None, ax=None, vars = None,
                  linewidth = 1.0, color = None, linestyle = None,
                  xStart = None, xEnd = None,
                  yMin = None, yMax = None,
                  xlabel = None, label = None, nLabels = None, title = None):
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
    ySize = 5.0 * nVars

    if fig is None:
        fig = plt.figure(figsize=(xSize, ySize))

        yBot = 0.1
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
        ax[iPlot].plot(times, data[var], \
                       linewidth = linewidth, \
                       linestyle = linestyle, \
                       label = label, \
                       color = color)
        if (xlabel == None):
            xlabel = var
        ax[iPlot].set_ylabel(xlabel)
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

    if outFile:
        print(" ==> Writing file : ", outFile)
        fig.savefig(outFile)
        plt.close(fig)
    else:
        return fig, ax
