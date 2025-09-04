#!/usr/bin/env python3

import numpy as np
from pyitm.fileio import variables

def extract_1d(sat_locations, model_data, extrapolate=False, verbose=False, interpVar=None):
    """ Fly a satellite through model results, returning a timeseries

    Parameters
    ----------
    sat_locations (dict): dictionary containing (at minimum) keys (time, lon, lat, alt)

    model_data (dict): dictionary returned from one of the read routines.
        Should contain info on time, grid, and any variables we want to interpolate.
    
    extrapolate (bool): whether to use data outside the time range covered by both data.

    verbose (bool): print debugging info? default=False

    interpVar (int or list): which variable (indices) to interpolate. default=None, so 
        interpolate every variable (except lon, lat, alt) in the model_data file

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
            raise ValueError("None of the satellite data and model outputs overlap!!")
        
        model_data['data'] = model_data['data'][t_ma_model, ...]
        model_data['times'] = np.array(model_data['times'])[t_ma_model]
        
    if timesliceSat:
        t_ma_sat = np.where((sat_locations['times'] >= t_min) 
                            & (sat_locations['times'] < t_max))[0]
        if len(t_ma_sat) == 0:
            raise ValueError("None of the satellite data and model outputs overlap!!")
        
        for varname in sat_locations.keys():
            # might not be a numpy array
            sat_locations[varname] = np.array(sat_locations[varname])[t_ma_sat]

    # setup for loop 
    itb4 = 0 # iTime before sat time

    lons = np.rad2deg(np.unique(model_data['lons']))
    lats = np.rad2deg(np.unique(model_data['lats']))
    alts = np.unique(model_data['alts'])

    dLon = lons[1] - lons[0]
    dLat = lats[1] - lats[0]
    nAlts = len(alts)

    if interpVar is None:
        # skip first 3 vars, usually lon,lat,alt
        interpVar = range(3, model_data['nvars'])
    if isinstance(interpVar, int):
        interpVar = [interpVar]

    outTimes = []
    # outData is a dict to simplify lookups.
    outVals = {iVar: [] for iVar in interpVar}

    for i, time in enumerate(sat_locations['times']):
        
        if model_data['times'][itb4 + 1] < time:
            itb4 += 1
            if itb4+1 == len(model_data['times']):
                itb4 -= 1

        if verbose:
            print(f"->> in loop. sat time: {time}, gitmtimes: {model_data['times'][itb4]},"
                  f"{model_data['times'][itb4+1]}")
        
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
                    raise ValueError("Above model altitude domain!")
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

    return outData, outTimes




