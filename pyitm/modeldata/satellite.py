#!/usr/bin/env python3

import numpy as np
import pandas as pd # for dealing with numpy vs python datetimes
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

    iVar (int or list): which variable (indices) to interpolate. default=None, so 
        interpolate everything.

    Returns
    -------
    (dict): dictionary with keys (time, lon, lat, alt, [etc.]), where etc. is any other 
        key that was present in the original satellite data (lst, quality flags, etc.),
        and the variable(s) that have been interpolated.

    """

    sat_locations['time'] = pd.to_datetime(sat_locations['time'])
    model_data['time'] = pd.to_datetime(model_data['time'])

    t_min = max(min(sat_locations['time']), min(model_data['time']))
    t_max = min(max(sat_locations['time']), max(model_data['time']))

    timesliceModel = False
    timesliceSat = False
    if min(model_data['time']) < t_min or max(model_data['time']) > t_max:
        timesliceModel = True
    if min(sat_locations['time']) < t_min or max(sat_locations['time']) > t_max:
        timesliceSat = True

    if verbose:
        print(f" -> found (tmin, tmax): ({t_min}, {t_max})")
        print(f" --> timesliceSat: {timesliceSat}, timesliceModel: {timesliceModel}")

    if timesliceModel:
        t_ma_model = np.where((model_data['time'] >= t_min)
                             & (model_data['time'] <= t_max))[0]
        for varname in model_data.keys():
            # hopefully already a numpy array...
            model_data[varname] = model_data[varname][t_ma_model, ...]
        
    if timesliceSat:
        t_ma_sat = np.where((sat_locations['time'] >= t_min) 
                               & (sat_locations['time'] <= t_max))
        print(t_ma_sat)
        for varname in sat_locations.keys():
            # might not be a numpy array
            sat_locations[varname] = np.array(sat_locations[varname])[t_ma_sat]

    # setup for loop 
    itb4 = 0 # iTime before sat time
    lons = np.rad2deg(np.unique(model_data['data'][0, 0, :,0,0]))
    lats = np.rad2deg(np.unique(model_data['data'][0, 1, 0,:,0]))
    alts = model_data['data'][0, 2, 0, 0, :]/1000.0
    print(lons, lats, alts)
    dLon = lons[1] - lons[0]
    dLat = lats[1] - lats[0]
    nAlts = len(alts)

    if interpVar is None:
        interpVar = range(3, model_data['nVars'])
    if isinstance(interpVar, int):
        interpVar = [interpVar]

    outTimes = []
    outVals = [[] for iVar in interpVar]

    for i, time in enumerate(sat_locations['time']):
        
        if model_data['time'][itb4 + 1] < time:
            itb4 += 1

        if verbose:
            print(f"->> in loop. sat time: {time}, gitmtimes: {model_data['time'][itb4]},{model_data['time'][itb4+1]}")
        
        dt = (model_data["time"][itb4] - \
            model_data["time"][itb4 + 1]).total_seconds()
        xt = (time - model_data["time"][itb4]).total_seconds() / dt

        xLon = (sat_locations['lon'][i] - lons[0])/dLon
        iLon = int(xLon)
        xLon = xLon - iLon
        
        yLat = (sat_locations['lat'][i] - lats[0])/dLat
        jLat = int(yLat)
        yLat = yLat - jLat

        kAlt = 0
        zAlt = 0.0
        if ((sat_locations['alt'][i] > alts[0]) and (nAlts > 1)):
            if (sat_locations['alt'][i] > alts[nAlts-1]):
                # above domain:
                kAlt = nAlts-2
                zAlt = 1.0
            else:
                while (alts[kAlt] < sat_locations['alt'][i]):
                    kAlt = kAlt + 1
                kAlt = kAlt - 1
                zAlt = (sat_locations['alt'][i] - alts[kAlt]) / (alts[kAlt+1] - alts[kAlt])
            kAltp1 = kAlt + 1
        else:
            kAltp1 = kAlt

        for iVar in interpVar:
            BeforeVal = \
                (1-xLon)*(1-yLat)*(1-zAlt)*model_data['data'][itb4, iVar, iLon, jLat, kAlt]+\
                (  xLon)*(1-yLat)*(1-zAlt)*model_data['data'][itb4, iVar, iLon+1, jLat, kAlt]+\
                (1-xLon)*(  yLat)*(1-zAlt)*model_data['data'][itb4, iVar, iLon, jLat+1, kAlt]+\
                (  xLon)*(  yLat)*(1-zAlt)*model_data['data'][itb4, iVar, iLon+1, jLat+1, kAlt]+\
                (1-xLon)*(1-yLat)*(  zAlt)*model_data['data'][itb4, iVar, iLon, jLat, kAltp1]+\
                (  xLon)*(1-yLat)*(  zAlt)*model_data['data'][itb4, iVar, iLon+1, jLat, kAltp1]+\
                (1-xLon)*(  yLat)*(  zAlt)*model_data['data'][itb4, iVar, iLon, jLat+1, kAltp1]+\
                (  xLon)*(  yLat)*(  zAlt)*model_data['data'][itb4, iVar, iLon+1, jLat+1, kAltp1]
            
            AfterVal = \
                (1-xLon)*(1-yLat)*(1-zAlt)*model_data['data'][itb4+1, iVar, iLon, jLat, kAlt]+\
                (  xLon)*(1-yLat)*(1-zAlt)*model_data['data'][itb4+1, iVar, iLon+1, jLat, kAlt]+\
                (1-xLon)*(  yLat)*(1-zAlt)*model_data['data'][itb4+1, iVar, iLon, jLat+1, kAlt]+\
                (  xLon)*(  yLat)*(1-zAlt)*model_data['data'][itb4+1, iVar, iLon+1, jLat+1, kAlt]+\
                (1-xLon)*(1-yLat)*(  zAlt)*model_data['data'][itb4+1, iVar, iLon, jLat, kAltp1]+\
                (  xLon)*(1-yLat)*(  zAlt)*model_data['data'][itb4+1, iVar, iLon+1, jLat, kAltp1]+\
                (1-xLon)*(  yLat)*(  zAlt)*model_data['data'][itb4+1, iVar, iLon, jLat+1, kAltp1]+\
                (  xLon)*(  yLat)*(  zAlt)*model_data['data'][itb4+1, iVar, iLon+1, jLat+1, kAltp1]
            
            outVals[iVar].append((1-xt) * BeforeVal + xt * AfterVal)

        outTimes.append(time)

    # make & fill output dictionary:
    outData = {'time': np.array(outTimes)}

    for iVar in interpVar:
        outData[variables.get_short_names(model_data['vars'][iVar])] = np.array(outVals[iVar])
    
    # in case any other data was present in original satellite file
    for inKey in sat_locations.keys():
        if inKey != time:
            outData[inKey] = np.array(sat_locations[inKey])

    return outData




