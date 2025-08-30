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

    axplot.text(0.0, maxR * 1.1, title,
                verticalalignment = 'bottom',
                horizontalalignment = 'center',
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
    cbar.set_label(cbar_label, rotation=90)

    return

