#!/usr/bin/env python3

import numpy as np

# ----------------------------------------------------------------------------
# take a dictionary containing all of the model data and
# return slices. The data should have the shape:
# [nTimes, nVars, nLons, nLats, nAlts]
# or
# [nTimes, nLons, nLats, nAlts]
# ----------------------------------------------------------------------------

def data_slice(allData3D, iLon = -1, iLat = -1, iAlt = -1):

    nTimes = allData3D['nTimes']
    nVars = allData3D['nVars']
    nLons = allData3D['nLons']
    nLats = allData3D['nLats']
    nAlts = allData3D['nAlts']

    if (nVars > 1):
        if (iAlt > -1):
            slices = np.zeros((nTimes, nVars, nLons, nLats))
            slices[:, :, :, :] = allData3D['data'][:, :, :, :, iAlt]
        elif (iLat > -1):
            slices = np.zeros((nTimes, nVars, nLons, nAlts))
            slices[:, :, :, :] = allData3D['data'][:, :, :, iLat, :]
        else:
            slices = np.zeros((nTimes, nVars, nLats, nAlts))
            slices[:, :, :, :] = allData3D['data'][:, :, iLon, :, :]
    else:
        if (iAlt > -1):
            slices = np.zeros((nTimes, nLons, nLats))
            slices[:, :, :] = allData3D['data'][:, :, :, iAlt]
        elif (iLat > -1):
            slices = np.zeros((nTimes, nLons, nAlts))
            slices[:, :, :] = allData3D['data'][:, :, iLat, :]
        else:
            slices = np.zeros((nTimes, nLats, nAlts))
            slices[:, :, :] = allData3D['data'][:, iLon, :, :]

    return slices

# ----------------------------------------------------------------------------
# This function calculates the edges of cells based on the centers of the cells
# it assumes a 1D array.
# ----------------------------------------------------------------------------

def move_centers_to_edges(pos):
    edges = (pos[1:] + pos[:-1])/2
    dpLeft = pos[1] - pos[0]
    dpRight = pos[-1] - pos[-2]
    edges = np.append(edges[0] - dpLeft, edges)
    edges = np.append(edges, dpRight + edges[-1])
    return edges


#-----------------------------------------------------------------------------
# vertically integrate the 3D data given the altitudes.
#-----------------------------------------------------------------------------

def vertically_integrate(value, alts, calc3D = False):
    [nLons, nLats, nAlts] = value.shape
    integrated = np.zeros((nLons, nLats, nAlts))
    descending = np.arange(nAlts-2, -1, -1)
    dz = alts[:,:,-1] - alts[:,:,-2]
    integrated[:,:,-1] = value[:,:,-1] * dz
    for i in descending:
        dz = alts[:,:,i+1] - alts[:,:,i]
        integrated[:,:,i] = integrated[:,:,i+1] + value[:,:,i] * dz
    if (not calc3D):
        integrated = integrated[:,:,0]
    return integrated

#-----------------------------------------------------------------------------
# Calculate TEC using the integration function
#  --> Assume that we are giving it the [e-] as a variable!
#-----------------------------------------------------------------------------

def calc_tec(allData3D):
    nTimes = allData3D['nTimes']
    nVars = allData3D['nVars']
    nLons = allData3D['nLons']
    nLats = allData3D['nLats']
    alts1d = allData3D['alts']

    slices = np.zeros((nTimes, nLons, nLats))
    for iTime in range(nTimes):
        tec2d = vertically_integrate(allData3D['data'][iTime, :, :, :], alts1d)
        # Convert km->m and /m3 to TECU
        slices[iTime, :, :] = tec2d[:, :] * 1000.0 / 1e16
    return slices



