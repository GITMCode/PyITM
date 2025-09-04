#!/usr/bin/env python3

import numpy as np

#----------------------------------------------------------------------------
# Test to see if altitude is changing as a function of lat/lon/block
#  - this is really to test for the dipole grid
#----------------------------------------------------------------------------

def calc_if_same_alts(altsWblock3d):

    isSame = True
    nBlocks = len(altsWblock3d[:, 0, 0, 0])

    nAlts = len(altsWblock3d[0, 0, 0, :])

    for iAlt in range(nAlts):
        reference = altsWblock3d[0, 0, 0, iAlt]
        small = 1e-6 * reference

        for iBlock in range(nBlocks):
            alts2d = altsWblock3d[iBlock, :, :, iAlt]
            if (np.abs(alts2d[0,0] - reference) > small):
                isSame = False
            if (np.abs(alts2d[-1,-1] - reference) > small):
                isSame = False
    return isSame

#----------------------------------------------------------------------------
# Test to see if grid is uniform:
#----------------------------------------------------------------------------

def calc_if_uniform_grid(lonsWblock2d, latsWblock2d):

    nBlocks = len(lonsWblock2d[:, 0, 0])

    # Let's figure out if we have a uniform horizontal grid:
    isUniform = True

    for iBlock in range(nBlocks):
        # Assume first 3 variables are lon, lat, alt:
        longitude = lonsWblock2d[iBlock, :, :]
        latitude = latsWblock2d[iBlock, :, :]

        if (iBlock == 0):
            dLon = longitude[1, 0] - longitude[0, 0]
            dLat = latitude[0, 1] - latitude[0, 0]
        else:
            dLonT = longitude[1, 0] - longitude[0, 0]
            dLatT = latitude[0, 1] - latitude[0, 0]
            if (np.abs(dLat - dLatT) > dLat/1000.0):
                isUniform = False
            if (np.abs(dLon - dLonT) > dLon/1000.0):
                isUniform = False

    return isUniform

#-----------------------------------------------------------------------------
# find cut in altitude, given that the altitude may vary from point to point
#-----------------------------------------------------------------------------

def find_alts_oneblock(alts3d, goalAlt):

    nLons, nLats, nAlts = np.shape(alts3d)

    iAlts = np.zeros((nLons, nLats)).astype('int')

    for iLon in range(nLons):
        for iLat in range(nLats):
            diff = np.abs(goalAlt - alts3d[iLon, iLat, :])
            iAlts[iLon, iLat] = np.argmin(diff)

    return iAlts


def find_alts(alts, goalAlt, blocks = False):

    sizes = np.shape(alts)
    nDims = len(sizes)

    if (nDims == 3):
        iAlts = find_alts_oneblock(alts, goalAlt)
        nBlocks = 0
    elif (nDims == 4):
        nBlocks, nLons, nLats, nAlts = sizes
        iAlts = np.zeros((nBlocks, nLons, nLats)).astype('int')
    else:
        print('Dont understand dimensions of alts in find_alts')
        iAlt = -1
        return iAlt

    for iBlock in range(nBlocks):
        iAlts[iBlock, :, :] = find_alts_oneblock(alts[iBlock, :, :, :], goalAlt)

    return iAlts

#-----------------------------------------------------------------------------
# find which cut direction to make, then extract a bunch of information
# based on that decision.
#-----------------------------------------------------------------------------

def find_cut(args, allData):

    cutValue = -1e32
    cutString = ''
    cutShort = ''
    iLon = -1
    iLat = -1
    iAlt = -1
    xLabel = ''
    yLabel = ''
    lons1d = allData['lons'][:, 0, 0]
    lats1d = allData['lats'][0, :, 0]
    alts1d = allData['alts'][0, 0, :]
    xRange = [0,0]
    yRange = [0,0]

    if (args.cut == 'alt'):
        altGoal = float(args.alt)
        diff = np.abs(alts1d - altGoal)
        iAlt = np.argmin(diff)
        cutValue = alts1d[iAlt]
        cutString = 'Alt : %d km' % int(cutValue)
        cutShort = 'alt%04d_' % iAlt
        xLabel = 'Longitude (deg)'
        yLabel = 'Latitude (deg)'
        xPos1d = lons1d
        yPos1d = lats1d
        xRange = [0, 360]
        yRange = [-90.0, 90.0]
    if (args.cut == 'lon'):
        lonGoal = float(args.lon)
        diff = np.abs(lons1d - lonGoal)
        iLon = np.argmin(diff)
        cutValue = lons1d[iLon]
        cutString = 'Lon : %d deg' % int(cutValue)
        cutShort = 'lon%04d_' % iLon
        xLabel = 'Latitude (deg)'
        yLabel = 'Altitude (km)'
        xPos1d = lats1d
        yPos1d = alts1d
        xRange = [-90.0, 90.0]
        yRange = [alts1d[2], alts1d[-3]]
    if (args.cut == 'lat'):
        latGoal = float(args.lat)
        diff = np.abs(lats1d - latGoal)
        iLat = np.argmin(diff)
        cutValue = lats1d[iLat]
        cutString = 'Lat : %d deg' % int(cutValue)
        cutShort = 'lat%04d_' % iLat
        xLabel = 'Longitude (deg)'
        yLabel = 'Altitude (km)'
        xPos1d = lons1d
        yPos1d = alts1d
        xRange = [0.0, 360.0]
        yRange = [alts1d[2], alts1d[-3]]

    cut = {'iLon': iLon,
           'iLat': iLat,
           'iAlt': iAlt,
           'cutValue': cutValue,
           'cutString': cutString,
           'cutShort': cutShort,
           'xLabel': xLabel,
           'yLabel': yLabel,
           'xRange': xRange,
           'yRange': yRange,
           'xPos': xPos1d,
           'yPos': yPos1d}

    return cut

# ----------------------------------------------------------------------------
# This takes a 3D array (lon, lat, alt) and returns a (lon, lat), where
# the altitude is (lon, lat) dependent
# ----------------------------------------------------------------------------

def slice_alt_with_array(var3d, iAlt):

    nLons, nLats, nAlts = np.shape(var3d)
    slice2d = np.zeros((nLons, nLats))

    for iLon in range(nLons):
        for iLat in range(nLats):
            iAlt_ = iAlt[iLon, iLat]
            slice2d[iLon, iLat] = var3d[iLon, iLat, iAlt_]
    return slice2d

# ----------------------------------------------------------------------------
# This takes a 4D array (block, lon, lat, alt) and returns a
# (block, lon, lat), where the altitude is (lon, lat) dependent
# ----------------------------------------------------------------------------

def slice_alt_with_array_block(var4d, iAlt):

    nBlocks, nLons, nLats, nAlts = np.shape(var4d)
    slice3d = np.zeros((nBlocks, nLons, nLats))

    for iBlock in range(nBlocks):
        for iLon in range(nLons):
            for iLat in range(nLats):
                iAlt_ = iAlt[iBlock, iLon, iLat]
                slice3d[iBlock, iLon, iLat] = var4d[iBlock, iLon, iLat, iAlt_]
    return slice3d

# ----------------------------------------------------------------------------
# This takes a 3D array (lon, lat, alt) and returns a (lat, alt) slice
# if iLon < 0, then it returns a zonal average
# ----------------------------------------------------------------------------

def slice_lon_3d(var3d, iLon):

    nLons, nLats, nAlts = np.shape(var3d)
    slice2d = np.zeros((nLats, nAlts))

    if (iLon >= 0):
        slice2d = var3d[iLon, :, :]
    else:
        for iL in range(nLons-4):
            slice2d = slice2d + \
                var3d[iLon + 2, :, :]
        slice2d = slices / (nLons-4)

    return slice2d

# ----------------------------------------------------------------------------
# This takes a 4D array (n, lon, lat, alt) and returns a (n, lat, alt) slice
# ----------------------------------------------------------------------------

def slice_lon_4d(var4d, iLon):

    nVals, nLons, nLats, nAlts = np.shape(var4d)
    slice3d = np.zeros((nVals, nLats, nAlts))

    for iVal in range(nVals):
        slice3d[iVal, :, :] = slice_lon_3d(var4d[iVal, :, :, :], iLon)
    
    return slice3d

# ----------------------------------------------------------------------------
# This takes a 5D array (n1, n2, lon, lat, alt) and returns a
# (n1, n2, lat, alt) slice
#  - if iLon < 0, then returns zonal average
# ----------------------------------------------------------------------------

def slice_lon_5d(var5d, iLon):

    nVal1, nVal2, nLons, nLats, nAlts = np.shape(var5d)
    slice4d = np.zeros((nVal1, nVal2, nLats, nAlts))

    for iVal1 in range(nVal1):
        for iVal2 in range(nVal2):
            slice4d[iVal1, iVal2, :, :] = slice_lon_3d(var5d[iVal1, iVal2, :, :, :], iLon)
    
    return slice4d

# ----------------------------------------------------------------------------
# take a dictionary containing all of the model data and
# return slices. The data should have the shape:
# [nTimes, nVars, nLons, nLats, nAlts]
# or
# [nTimes, nLons, nLats, nAlts]
# ----------------------------------------------------------------------------

def data_slice(allData3D, iLon = -1, iLat = -1, iAlt = -1):

    nTimes = allData3D['ntimes']
    nVars = allData3D['nvars']
    nLons = allData3D['nlons']
    nLats = allData3D['nlats']
    nAlts = allData3D['nalts']
    nBlocks = allData3D['nblocks']

    print(nTimes, nVars, nBlocks, nLons, nLats, nAlts)
    
    doAltCut = False
    altArray = False
    if (np.isscalar(iAlt)):
        if (iAlt > -1):
            doAltCut = True
    else:
        doAltCut = True
        altArray = True
    
    # This has become somewhat complicated to handle all of the different cases:
    # nVars == 1 vs nVars > 1
    # nBlocks = 0 vs block-based arrays
    # iAlt = scalar or iAlt is an array - We can now slice with different 
    #                  iAlts for each lat/lon point

    if (nVars > 1):
        if (doAltCut):
            if (nBlocks == 0):
                slices = np.zeros((nTimes, nVars, nLons, nLats))
            else:
                slices = np.zeros((nTimes, nVars, nBlocks, nLons, nLats))
            if (not altArray):
                if (nBlocks == 0):
                    slices[:, :, :, :] = allData3D['data'][:, :, :, :, iAlt]
                else:
                    slices[:, :, :, :, : ] = allData3D['data'][:, :, :, :, :, iAlt]
            else:
                if (nBlocks == 0):
                    for iLon in range(nLons):
                        for iLat in range(nLats):
                            iAlt_ = iAlt[iLon, iLat]
                            slices[:, :, iLon, iLat] = allData3D['data'][:, :, iLon, iLat, iAlt_]
                else:    
                    for iBlock in range(nBlocks):
                        for iLon in range(nLons):
                            for iLat in range(nLats):
                                iAlt_ = iAlt[iBlock, iLon, iLat]
                                slices[:, :, iBlock, iLon, iLat] = \
                                    allData3D['data'][:, :, iBlock, iLon, iLat, iAlt_]
        elif (iLat > -1):
            slices = np.zeros((nTimes, nVars, nLons, nAlts))
            slices[:, :, :, :] = allData3D['data'][:, :, :, iLat, :]
        else:
            slices = np.zeros((nTimes, nVars, nLats, nAlts))
            if (iLon > 0):
                slices[:, :, :, :] = allData3D['data'][:, :, iLon, :, :]
            else:
                for iL in range(nLons-4):
                    slices[:, :, :, :] = slices[:, :, :, :] + \
                        allData3D['data'][:, :, iLon, :, :]
                slices[:, :, :, :] = slices[:, :, :, :] / (nLons-4)
    else:
        if (doAltCut):
            if (nBlocks == 0):
                slices = np.zeros((nTimes, nLons, nLats))
            else:
                slices = np.zeros((nTimes, nBlocks, nLons, nLats))

            if (not altArray):
                if (nBlocks == 0):
                    slices[:, :, :] = allData3D['data'][:, :, :, iAlt]
                else:
                    slices[:, :, :, : ] = allData3D['data'][:, :, :, :, iAlt]
            else:
                if (nBlocks == 0):
                    for iTime in range(nTimes):
                        slices[iTime, :, :] = \
                            slice_alt_with_array(allData3D['data'][iTime, :, :, :], iAlt)
                else:
                    for iTime in range(nTimes):
                        slices[iTime, :, :, :] = \
                            slice_alt_with_array_block(allData3D['data'][iTime, :, :, :], iAlt)

        elif (iLat > -1):
            slices = np.zeros((nTimes, nLons, nAlts))
            slices[:, :, :] = allData3D['data'][:, :, iLat, :]
        else:
            if (nBlocks == 0):
                slices = slice_lon_4d(allData3D['data'], iLon)
            else:
                slices = slice_lon_5d(allData3D['data'], iLon)

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

# ----------------------------------------------------------------------------
# This function calculates the edges of cells based on the centers of the cells
# it assumes a 2D array.
# - first is finds the edges in the X direction
# - then it finds the edges in the Y direction
# ----------------------------------------------------------------------------

def move_centers_to_corners(pos2d):

    nX, nY = np.shape(pos2d)
    intermediate = np.zeros((nX + 1, nY))
    final = np.zeros((nX + 1, nY + 1))
    for iY in range(nY):
        intermediate[:, iY] = move_centers_to_edges(pos2d[:,iY])
    for iX in range(nX + 1):
        final[iX, :] = move_centers_to_edges(intermediate[iX, :])
    return final


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
    nTimes = allData3D['ntimes']
    nVars = allData3D['nvars']
    nLons = allData3D['nlons']
    nLats = allData3D['nlats']
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
    

