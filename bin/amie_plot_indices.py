#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from pylab import cm
import sys
import argparse

from pyitm.fileio import amieio
from pyitm.modeldata import utils
from pyitm.plotting import line_plots

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(description = 'Plot AMIE files')
    parser.add_argument('files', metavar = 'file', nargs = '+', \
                        help = 'Files to process')

    parser.add_argument('-plotfile',
                        help = 'output file for plot',
                        default = 'amie_timeline.png')    

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# Main Code
# ----------------------------------------------------------------------

args = parse_args()
plt.rcParams.update({'font.size': 14})

plotfile = args.plotfile

file = args.files[0]
data = amieio.amie_read_binary(file)

iPot_ = -1
iEflux_ = -1

vars = data["Vars"]
for i, v in enumerate(vars):
    print(i, v)
    if ((iPot_ < 0) and ('potential' in v.lower())):
        iPot_ = i
        print('  --> found potential')
    if ((iEflux_ < 0) and ('electron energy flux' in v.lower())):
        iEflux_ = i
        print('  --> found electron eflux')
    if ((iEflux_ < 0) and ('auroral energy flux' in v.lower())):
        iEflux_ = i
        print('  --> found electron eflux')

if (iPot_ < 0):
    print("Can't find potential!")
    exit()

pot_ = vars[iPot_]
eflux_ = vars[iEflux_]

# convert to M/m2
eflux3d = np.array(data[eflux_]) / 1000.0
potential3d = np.array(data[pot_]) / 1000.0

nTimes = len(data["times"])
lats = data["lats"]
mlts = data["mlts"]

lats_edges = utils.move_centers_to_edges(lats)
mlts_edges = utils.move_centers_to_edges(mlts)

dtom = (6372.0 + 120.0) * 1000.0 * 2 * np.pi / 360.0
dlats = np.abs(lats_edges[1:] - lats_edges[:-1]) * dtom
dmlts = (mlts_edges[1:] - mlts_edges[:-1]) * 15.0 * dtom

dmlts2d, dlats2d = np.meshgrid(dmlts, dlats)
mlts2d, lats2d = np.meshgrid(mlts, lats * np.pi / 180.0)

area2d = np.cos(lats2d) * dmlts2d * dlats2d
areasum = np.sum(area2d)

cpcp = np.zeros(nTimes)
hp = np.zeros(nTimes)

for iT in range(nTimes):
    cpcp[iT] = np.max(potential3d[iT]) - np.min(potential3d[iT])
    hp[iT] = np.sum(eflux3d[iT] * area2d) / 1e9

fig = plt.figure(figsize = (11,6))
ax = fig.add_axes([0.07, 0.09, 0.9, 0.86])
ax = fig.add_axes([0.07, 0.09, 0.9, 0.86])
ax.plot(data["times"], cpcp)
#color = color, linestyle = line, linewidth = linewidth)

data_to_plot = {'CPCP (kv)': cpcp,
                'HP (GW)': hp,
                'times': data['times'],
                'vars': ['CPCP (kv)', 'HP (GW)'],
                'title': 'Cross Polar Cap Potential and Hemispheric Power',
                'note': file}

line_plots.lineplot_data(data_to_plot, plotfile)

