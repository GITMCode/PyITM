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


#-----------------------------------------------------------------------------
# This function takes either a low-pass filter or a high-pass filter
# of the slices. It assumes that:
#  allSlices = np.array of size (nTimes, nX, nY)
#  iLow or iHigh is set to the number of times to average (centered on current time)
#  if iLow is set, then it calculates a center-averaged valaue (low pass)
#  if iHigh is set, then it subtracts the center-averaged value (high pass)
#-----------------------------------------------------------------------------

def filter_slices(allSlices, iLowPass = -1, iHighPass = -1):
    
    nTimes, nXs, nYs = np.shape(allSlices)
    nAve = 0
    if (iLowPass > 0):
        nAve = iLowPass
        doSubtract = False
        print(' --> Low Pass Filtering the Slices!')
    else:
        nAve = iHighPass
        doSubtract = True
        print(' --> High Pass Filtering the Slices!')
    if (nAve % 2 == 0):
        print('high/low is set to even number, which is not great.')
        print('Consider setting it to odd number!')
    iHalf = int(nAve/2)
    allSlicesAve = np.zeros((nTimes, nXs, nYs))
    for iTime in range(nTimes):
        iLow = iTime - iHalf
        iHigh = iTime + iHalf + 1
        if (iLow < 0):
            iLow = 0
        if (iHigh > nTimes):
            iHigh = nTimes
        nT = 0
        sumSlice = np.zeros((nXs, nYs))
        for iTimeSub in np.arange(iLow, iHigh):
            sumSlice = sumSlice + allSlices[iTimeSub, :, :]
            nT += 1
        if (iLowPass > 0):
            allSlicesAve[iTime, :, :] = sumSlice / nT
        elif (iHighPass > 0):
            if (nT == nAve):
                allSlicesAve[iTime, :, :] = allSlices[iTime, :, :] - (sumSlice / nT)
        else:
            # no filter, just move the data!
            allSlicesAve[iTime, :, :] = allSlices[iTime, :, :]
            
    return allSlicesAve
    

