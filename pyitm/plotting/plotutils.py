#!/usr/bin/env python3
""" This is a simple lat/lon plotter
"""

import matplotlib as mpl
import numpy as np

# ----------------------------------------------------------------------------
# Get max and min values
# model_data is [time, xpositions, ypositions]
# ----------------------------------------------------------------------------

def get_min_max_data(allSlices, yPos, \
                     yMin = -1e32, yMax = 1e32, \
                     color = 'default', \
                     minVal = 1e32, maxVal = -1e32,
                     isLog = False):
    
    symmetric = False
    if (color == 'default'):
        cmap = mpl.cm.plasma
    if (color == 'red'):
        cmap = mpl.cm.YlOrRd

    mask = ((yPos >= yMin) & (yPos <= yMax))
    if (mask.max()):    
        doPlot = True
        maxi = allSlices[:, :, mask].max() * 1.01
        mini = allSlices[:, :, mask].min() * 0.99
    
        if ((mini < 0.0) and (not isLog)):
            symmetric = True
            cmap = mpl.cm.bwr
            maxi = abs(allSlices[:, :, mask]).max() * 1.05
            mini = -maxi

    else:
        doPlot = False
        mini = -1.0
        maxi = 1.0

    if (minVal < 1e31):
        mini = minVal
    if (maxVal > -1e31):
        maxi = maxVal
        
    min_max_data = {'mini' : mini,
                    'maxi' : maxi,
                    'cmap' : cmap,
                    'symmetric' : symmetric,
                    'doPlot': doPlot,
                    'mask': mask}
    
    return min_max_data
