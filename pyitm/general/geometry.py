#!/usr/bin/env python
# Copyright 2025, the PyITM Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md

"""Routines to do stuff with geometry."""

import numpy as np

#-----------------------------------------------------------------------------
# alt and rPlanet are in km
# lons and lats are in degrees
#-----------------------------------------------------------------------------

def calc_areas(lons2d, lats2d, alt, rPlanet = 6372.0):

    dlons2d = lons2d * np.pi / 360.0
    dlons2d[1:-1, :] = lons2d[2:, :] - lons2d[0:-2, :]
    dlons2d[0, :] = dlons2d[1, :]
    dlons2d[-1, :] = dlons2d[-2, :]
    dlats2d = lats2d * np.pi / 360.0
    dlats2d[:, 1:-1] = lats2d[:, 2:] - lats2d[:, 0:-2]
    dlats2d[:, 0] = dlats2d[:, 1]
    dlats2d[:, -1] = dlats2d[:, -2]

    r = (rPlanet + alt) * 1000.0
    area = r * r * dlons2d * dlats2d * np.cos(lats2d * np.pi / 360.0)

    return area
    
