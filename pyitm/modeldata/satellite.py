#!/usr/bin/env python3

import numpy as np
from pyitm.fileio import variables

def extract_1d(sat_locations, model_data, interpVar=None,
               extrapolate=False, verbose=False):
    """ Fly a satellite through model results, returning a timeseries

    Parameters
    ----------
    sat_locations (dict): dictionary containing (at minimum) keys (time, lon, lat, [alt]).
        If alt is not present, model data will be altitude integrated (TEC). All other
        variables present will be returned in the output dict.

    model_data (dict): dictionary returned from one of the read routines.
        Should contain info on time, grid, and any variables we want to interpolate.
    
    interpVar (int or list of int's): which variable (indices) to interpolate. default=None, so 
        interpolate every variable (except lon, lat, [alt]) in the model_data file.

    extrapolate (bool): whether to use data outside the time range covered by both data.

    verbose (bool): print debugging info? default=False

    Returns
    -------
    (dict): dictionary with keys (time, lon, lat, alt, [etc.]), where etc. is any other 
        key that was present in the original satellite data (lst, quality flags, etc.),
        and the variable(s) that have been interpolated.

    """

    sat_locations['times'] = np.array(sat_locations['times'])
    model_data['time'] = np.array(model_data['times'])

    t_min = max(min(sat_locations['times']), min(model_data['times']))
    t_max = min(max(sat_locations['times']), max(model_data['times']))

    timesliceModel = False
    timesliceSat = False
    if min(model_data['times']) < t_min or max(model_data['times']) > t_max:
        timesliceModel = True
    if min(sat_locations['times']) < t_min or max(sat_locations['times']) > t_max:
        timesliceSat = True

    if verbose:
        print(f" -> found (tmin, tmax): ({t_min}, {t_max})")
        print(f" --> timesliceSat: {timesliceSat}, timesliceModel: {timesliceModel}")

    if timesliceModel:
        t_ma_model = np.where((model_data['time'] >= t_min)
                             & (model_data['time'] <= t_max))[0]
        if len(t_ma_model) == 0:
            raise ValueError("None of the satellite data and model outputs overlap!!"
                            f"min/max sat: {sat_locations['times'][0]} / {sat_locations['times'][-1]}"
                            f"min/max model: {model_data['times'][0]} / {model_data['times'][-1]}")
        
        model_data['data'] = model_data['data'][t_ma_model, ...]
        model_data['times'] = np.array(model_data['times'])[t_ma_model]
        
    if timesliceSat:
        t_ma_sat = np.where((sat_locations['times'] >= t_min) 
                            & (sat_locations['times'] < t_max))[0]
        if len(t_ma_sat) == 0:
            raise ValueError("None of the satellite data and model outputs overlap!!"
                            f"min/max sat: {sat_locations['times'][0]} / {sat_locations['times'][-1]}"
                            f"min/max model: {model_data['times'][0]} / {model_data['times'][-1]}")
        
        for varname in sat_locations.keys():
            # might not be a numpy array
            if 'name' not in varname:
                sat_locations[varname] = np.array(sat_locations[varname])[t_ma_sat]

    # setup for loop 
    itb4 = 0 # iTime before sat time

    lons = np.sort(np.unique(model_data['lons']))
    lats = np.sort(np.unique(model_data['lats']))
    # Will be present in 2D files, but only one value
    alts = np.sort(np.unique(model_data['alts']))

    # rad to degrees if needed
    if all(np.abs(lons) < 7):
        lons = np.rad2deg(lons)
    if all(np.abs(lats) < 7):
        lats = np.rad2deg(lats)

    dLon = lons[1] - lons[0]
    dLat = lats[1] - lats[0]

    if 'alts' in sat_locations.keys():
        nAlts = len(alts)
    else:
        nAlts = 1

    if interpVar is None:
        # skip first 3 vars, usually lon, lat, alt
        interpVar = []
        for iVar, varname in enumerate(model_data['vars']):
            if varname not in ['Longitude', 'Latitude', 'Altitude']:
                interpVar.append(iVar)

    else:
        if isinstance(interpVar, int): #single int
            interpVar = [interpVar]
            
    # check if we only read in one variable. Reshape model_data if so
    if len(model_data['data'].shape) == 4:
        t, x, y, z = model_data['data'].shape
        model_data['data'] = model_data['data'].reshape((t, 1, x, y, z))

    # outData is a dict to simplify lookups.
    outVals = {iVar: [] for iVar in interpVar}

    if verbose:
        print(f" -> Interpolating variables: {[model_data['vars'][i] for i in interpVar]}.")

    for i, time in enumerate(sat_locations['times']):
        
        if model_data['times'][itb4 + 1] <= time:
            itb4 += 1
            if itb4+1 == len(model_data['times']):
                itb4 -= 1
        
        dt = (model_data["times"][itb4] - \
              model_data["times"][itb4 + 1]).total_seconds()
        xt = (time - model_data["times"][itb4]).total_seconds() / dt

        xLon = (sat_locations['lons'][i] - lons[0])/dLon
        iLon = int(xLon)
        xLon = xLon - iLon
        
        yLat = (sat_locations['lats'][i] - lats[0])/dLat
        jLat = int(yLat)
        yLat = yLat - jLat

        kAlt = 0
        zAlt = 0.0
        if ((sat_locations['alts'][i] > alts[0]) and (nAlts > 1)):
            if (sat_locations['alts'][i] > alts[nAlts-1]):
                if extrapolate:
                    # above domain:
                    kAlt = nAlts-2
                    zAlt = 1.0
                else:
                    raise ValueError(
                        "Above model altitude domain!\n"
                        f" Max 3 modeled altitudes:{np.round(alts[-3:], 2)}\n"
                        f" Min/Max satellite altitude: "
                        f"({round(np.min(sat_locations['alts']), 1)}/{round(np.max(sat_locations['alts']),1)})")
            else:
                while (alts[kAlt] < sat_locations['alts'][i]):
                    kAlt = kAlt + 1
                kAlt = kAlt - 1
                zAlt = (sat_locations['alts'][i] - alts[kAlt]) / (alts[kAlt+1] - alts[kAlt])
            kAltp1 = kAlt + 1
        else:
            kAltp1 = kAlt

        for iVar in interpVar:
            BeforeVal = \
                (1-xLon)*(1-yLat)*(1-zAlt)*model_data['data'][itb4, iVar, iLon,   jLat,   kAlt]+\
                (  xLon)*(1-yLat)*(1-zAlt)*model_data['data'][itb4, iVar, iLon+1, jLat,   kAlt]+\
                (1-xLon)*(  yLat)*(1-zAlt)*model_data['data'][itb4, iVar, iLon,   jLat+1, kAlt]+\
                (  xLon)*(  yLat)*(1-zAlt)*model_data['data'][itb4, iVar, iLon+1, jLat+1, kAlt]+\
                (1-xLon)*(1-yLat)*(  zAlt)*model_data['data'][itb4, iVar, iLon,   jLat,   kAltp1]+\
                (  xLon)*(1-yLat)*(  zAlt)*model_data['data'][itb4, iVar, iLon+1, jLat,   kAltp1]+\
                (1-xLon)*(  yLat)*(  zAlt)*model_data['data'][itb4, iVar, iLon,   jLat+1, kAltp1]+\
                (  xLon)*(  yLat)*(  zAlt)*model_data['data'][itb4, iVar, iLon+1, jLat+1, kAltp1]
            
            AfterVal = \
                (1-xLon)*(1-yLat)*(1-zAlt)*model_data['data'][itb4+1, iVar, iLon,   jLat,   kAlt]+\
                (  xLon)*(1-yLat)*(1-zAlt)*model_data['data'][itb4+1, iVar, iLon+1, jLat,   kAlt]+\
                (1-xLon)*(  yLat)*(1-zAlt)*model_data['data'][itb4+1, iVar, iLon,   jLat+1, kAlt]+\
                (  xLon)*(  yLat)*(1-zAlt)*model_data['data'][itb4+1, iVar, iLon+1, jLat+1, kAlt]+\
                (1-xLon)*(1-yLat)*(  zAlt)*model_data['data'][itb4+1, iVar, iLon,   jLat,   kAltp1]+\
                (  xLon)*(1-yLat)*(  zAlt)*model_data['data'][itb4+1, iVar, iLon+1, jLat,   kAltp1]+\
                (1-xLon)*(  yLat)*(  zAlt)*model_data['data'][itb4+1, iVar, iLon,   jLat+1, kAltp1]+\
                (  xLon)*(  yLat)*(  zAlt)*model_data['data'][itb4+1, iVar, iLon+1, jLat+1, kAltp1]

            outVals[iVar].append((1-xt) * BeforeVal + xt * AfterVal)

    # make & fill output dictionary:
    outData = {}

    # First, a list of some expected data that we want. some may exist, some may not.
    wantAlways = ['times', 'lons', 'lats', 'alts', 'lst', 'sat_name']
    for key in wantAlways:
        if key in sat_locations.keys():
            outData[key] = sat_locations[key]

    # Put any interpolated outputs into the file!
    for iVar in interpVar:
        dat_name = 'model_' + variables.get_short_names(model_data['vars'][iVar])[0]
        outData[dat_name] = np.array(outVals[iVar])
    
    # This adds density, wind, flags, etc. present in original satellite file
    for inKey in sat_locations.keys():
        if inKey not in wantAlways:
            outData['sat_' + inKey] = np.array(sat_locations[inKey])

    if verbose:
        print(f" -> Interpolation done! Returning dict with keys:\n\t{outData.keys()}")

    return outData

#-----------------------------------------------------------------------------
# finds the index of the time within a time array 
#-----------------------------------------------------------------------------

def find_index(time, t):

    if (t < time[0]):
        return 0
    if (t > time[-1]):
        return len(time)-1

    iLow = 0
    iHigh = len(time)
    iMid = int((iHigh + iLow)/2)
    while (iHigh - iLow > 1):
        if (time[iMid] == t):
            iHigh = iMid
            iLow = iMid
        else:
            if (t > time[iMid]):
                iLow = iMid
            else:
                iHigh = iMid
            iMid = int((iHigh + iLow)/2)
    return iMid


def calc_period(satData, verbose=False):
    """
    Determine the orbital period of the satellite, and add it to satData dict

    Inputs
    ------
        satData (dict) - satellite data read from satelliteio
        verbose (bool=False) - print extra info?

    Returns
    -------
        (dict) same as input but has 'orbital_period' key added

    Notes
    -----
        This is not the most robust, however will work for most satellites/planets.
    """

    # Find where the satellite crosses the north pole (latitudes go from + to -)
    # Then take the median time between those crossings as the period
    signLats = np.sign(np.diff(satData['lats']))
    
    iAtNorthPole = []
    tAtNorthPole = []
    for i in range(2, len(satData['lats']) - 2):
        if signLats[i] > signLats[i + 1]:
            iAtNorthPole.append(i)
            tAtNorthPole.append(satData['times'][i])
    try:
        iPeriod = int(np.round(np.median(np.diff(iAtNorthPole))))
        tPeriod = np.median(np.diff(tAtNorthPole))
    except ValueError:
        # This is mostly for the test cases... May be useful for synthetic data too
        print("Times too short! Cannot orbit average. Sorry.")
        return satData

    satData['orbital_period'] = tPeriod
    if verbose: 
        print(f"Found satellite orbital period of {tPeriod}, or {iPeriod} indices")
    return satData


def orbit_average(satData, varlist=None, verbose=False):
    """
    Attempt to smooth satData using a rolling average over each orbital period

    Inputs
    ------
        satData (dict) - satellite data read from satelliteio. If 'orbital_period'
            key is not present, will attempt to calculate it.
        varlist (list-like) - strings corresponding to keys in satData to smooth.
            If None, smooth anything with sat_* or model_* (default=None)
        verbose (bool=False) - print extra info?

    Returns
    -------
        (dict) same as input but has smoothed values of each variable requested
    
    """

    if 'orbital_period' not in satData.keys():
        satData = calc_period(satData, verbose=verbose)
        if 'orbital_period' not in satData.keys():
            if verbose:
                print("No orbital period found. Cannot smooth.")
            return satData
    tPeriod = satData['orbital_period']

    if varlist is None:
        varlist = []
        for i in satData.keys():
            if ('sat_' in i or 'model_' in i) and 'name' not in i:
                varlist.append(i)

    strs2ignore_smooth = ['flag','arg']
    for var in varlist:
        dosmooth = True
        for badvar in strs2ignore_smooth:
            # If we have a variable we want to ignore...
            if badvar in var.lower():
                dosmooth = False
        
        if dosmooth:
            if verbose:
                print("Smoothing ", var)
            smoothed = np.zeros(len(satData['times']))
            for i, t in enumerate(satData['times']):
                iMin = find_index(satData['times'], t-tPeriod/2)
                iMax = find_index(satData['times'], t+tPeriod/2)
                s = np.mean(satData[var][iMin : iMax+1])
                smoothed[i] = s
            satData['smoothed_' + var] = smoothed

    return satData

def calc_wind_dir(lons, lats):

    dlats = lats[1:] - lats[:-1]
    dlats = np.concatenate((dlats, [dlats[-1]]))

    dlons = lons[1:] - lons[:-1]
    dlons = np.concatenate((dlons, [dlons[-1]]))

    # Longitude can go across the 0 - 360 or 360 - 0 boundary, so
    # we need to correct for this possibility:
    
    dlons[dlons > 180.0] = dlons[dlons > 180.0] - 360.0
    dlons[dlons < -180.0] = dlons[dlons < -180.0] + 360.0

    # Longitudes get closer together near the poles, so we need to
    # correct for that also:
    dlons = dlons * np.cos(lats * np.pi / 180.0)

    # Make a unit vector of the direction of travel:
    mag = np.sqrt(dlats**2 + dlons**2)
    unitn = dlats / mag 
    unite = dlons / mag

    # Rotate the unit vector, so that it points orthogonal to the
    # orbit plane.  This will be the actual wind vector direction:
    uniteRot = unitn
    unitnRot = -unite

    return uniteRot, unitnRot

def calc_zon_merid_wind(satDataDict, verbose=False):
    """
    Takes satDataDict and transforms ve/vn columns to zonal & meridional winds. Assumes
    satellite's ve/vn are already aligned to the satellite trajectory (cross-track)
    
    """
    
    unitE, unitN = calc_wind_dir(satDataDict['lons'], satDataDict['lats'])

    if 'model_Ve' and 'model_Vn' in satDataDict:
        satDataDict['model_Ve'] = unitE * satDataDict['model_Ve']
        satDataDict['model_Vn'] = unitN * satDataDict['model_Vn']
    elif 'model_Vie' and 'model_Vin' in satDataDict:
        satDataDict['model_Vie'] = unitE * satDataDict['model_Vie']
        satDataDict['model_Vin'] = unitN * satDataDict['model_Vin']
    else:
        raise KeyError(
            f"Horizontal winds not found! Found keys:\n\t{satDataDict.keys()}")

    return satDataDict

