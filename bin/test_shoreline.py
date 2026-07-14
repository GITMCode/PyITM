#!/usr/bin/env python3
""" This is a simple lat/lon plotter
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from pyitm.modeldata import shoreline

fig = plt.figure(figsize=(10, 10))
ax = fig.add_axes([0.07, 0.07, 0.9, 0.9])

lon, lat = shoreline.shoreline()

ax.plot(lon,lat)
ax.set_aspect(1.0)

outFile = 'test.png'
print(" ==> Writing file : ", outFile)
fig.savefig(outFile)
plt.close(fig)


