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

def lineplot_data(data, outFile, vars = None):

    if (vars == None):
        vars = data['vars']
    nVars = len(vars)

    xSize = 10.0
    ySize = 5.0 * nVars

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
    for iPlot, var in enumerate(vars):
        ax[iPlot].plot(times, data[var])
        ax[iPlot].set_ylabel(var)
    
    sTime = times[0].strftime("%d %b %Y %H:%M:%S UT")
    eTime = times[-1].strftime("%d %b %Y %H:%M:%S UT")
    ax[-1].set_xlabel('Time from ' + sTime + ' to ' + eTime)

    if ('title' in data.keys()):
        ax[0].set_title(data['title'])

    print(" ==> Writing file : ", outFile)
    fig.savefig(outFile)
    plt.close(fig)
    
