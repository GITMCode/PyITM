#!/usr/bin/env python

# Functions to manually set the sizes of plots in python, since the
# python defaults suck

import matplotlib
matplotlib.use('AGG')
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
# This gives the axes for multi-row, but single column plots.
# Inputs:
#    - figIn: this is the figure (e.g., fig = plt.figure(figsize=(10, 10)))
#    - nPlots: number of plots in total
#    - yBottom: space desired below the bottom plot
#    - yTop: space desired above the top plot
#    - yBuffer: space desired between each plot
# Optional Inputs:
#    - xLeft: space desired to the left of the plots
#    - xRight: space desired to the right of the plots
# Outputs:
#    - ax: array of axes
# -------------------------------------------------------------------

def get_axes_one_column(figIn,
                        nPlots,
                        yBottom,
                        yTop,
                        yBuffer,
                        xLeft = 0.11,
                        xRight = 0.02):
    
    # plot size in y direction:
    ySize = (1.0 - yBottom - yTop) / nPlots - yBuffer * (nPlots - 1) / nPlots
    xSize = 1.0 - xLeft - xRight
    
    ax = []
    for iPlot in range(nPlots):
        # I want to reverse this, so that 0 is the top plot:
        i = nPlots - iPlot - 1
        ax.append(figIn.add_axes([xLeft,
                                  yBottom + i * (ySize + yBuffer),
                                  xSize,
                                  ySize]))
    
    return ax

# -------------------------------------------------------------------
# This gives the axes for multi-row, multi-column plots.
# Inputs:
#    - figIn: this is the figure (e.g., fig = plt.figure(figsize=(10, 10)))
#    - nPlots: number of plots in total
#    - yBottom: space desired below the bottom plot
#    - yTop: space desired above the top plot
#    - yBuffer: space desired between each plot
# Optional Inputs:
#    - xLeft: space desired to the left of the plots
#    - xRight: space desired to the right of the plots
# Outputs:
#    - ax: array of axes
# -------------------------------------------------------------------

def get_axes(figIn,
             nPlots,
             yBottom,
             yTop,
             yBuffer,
             xLeft,
             xRight,
             xBuffer):
    
    n = np.sqrt(nPlots)
    nX = int(np.ceil(n))
    nY = int(np.ceil(nPlots / nX))

    # plot size in y direction:
    ySize = (1.0 - yBottom - yTop) / nY - yBuffer * (nY - 1) / nY
    xSize = (1.0 - xLeft - xRight) / nX - xBuffer * (nX - 1) / nX
    
    ax = []
    for iPlot in range(nPlots):
        iX = int(iPlot % nX)
        iY = int(iPlot / nX)
        # I want to reverse y, so that 0 is the top left plot:
        iY = nY - iY - 1
        ax.append(figIn.add_axes([xLeft + iX * (xSize + xBuffer),
                                  yBottom + iY * (ySize + yBuffer),
                                  xSize,
                                  ySize]))
    
    return ax
