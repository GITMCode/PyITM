#!/usr/bin/env python3
""" This is for satellite plots.
"""


import os
import numpy as np
from pyitm.modeldata import satellite
from pyitm.general import time_conversion
import matplotlib.pyplot as plt
import scipy

def makesatplot(satDataDict, savepath, verbose=False, saveName=None):

    fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    didSmooth = 'smoothed_sat_rho' in satDataDict.keys() 
    
    if didSmooth:
        axs[0].plot(satDataDict['times'], satDataDict['model_rho'], 'b', 
                    alpha = 0.1)
        axs[0].plot(satDataDict['times'], satDataDict['sat_rho'], 'r', 
                    alpha = 0.1)
        axs[0].plot(satDataDict['times'], satDataDict['smoothed_model_rho'], 'b', 
                    label="Modeled Density")
        axs[0].plot(satDataDict['times'], satDataDict['smoothed_sat_rho'], 'r', 
                    label=satDataDict['sat_name'].upper() + " Density")
    else:
        axs[0].plot(satDataDict['times'], satDataDict['model_rho'], 'b', 
                    label="Modeled Density")
        axs[0].plot(satDataDict['times'], satDataDict['sat_rho'], 'r', 
                    label=satDataDict['sat_name'].upper() + " Density")
    axs[0].legend()

    
    if didSmooth:
        axs[1].plot(satDataDict['times'], 
                    100*(satDataDict['sat_rho']- satDataDict['model_rho'])
                    /satDataDict['sat_rho'],
                    alpha = 0.1, color='k')
        axs[1].plot(satDataDict['times'], 
                    100*(satDataDict['smoothed_sat_rho']-satDataDict['smoothed_model_rho'])
                    /satDataDict['smoothed_sat_rho'],
                    color='k')
    else:
        axs[1].plot(satDataDict['times'], 
                    100*(satDataDict['sat_rho']- satDataDict['model_rho'])
                    /satDataDict['sat_rho'],
                    color='k')
    
    axs[1].hlines(0, min(satDataDict['times']), max(satDataDict['times']),
                  linestyle='--', color='k')

    axs[1].set_ylabel(' Density (kg/m3)')
    axs[1].set_ylabel(' Diff (%)')
    fig.suptitle(satDataDict['sat_name'].upper())

    i = 0
    outdate = satDataDict['times'][0].strftime('%Y%m%d')
    savename = os.path.join(savepath, f"{satDataDict['sat_name']}_density_{outdate}_000.png")
    
    old_end = str(i).rjust(3, '0') + '.png'
    while os.path.exists(savename):
        # start adding zeros!
        i += 1
        new_end = str(i).rjust(3, '0') + '.png'
        savename = savename.replace(old_end, new_end)
        old_end = new_end
    
    print('--> Saving plot as: ' + savename)
    plt.savefig(savename)

    return



#-----------------------------------------------------------------------------  
def gridData4map(satDataDict, var='sat_rho',
                 lat_lim=None, dlat=2.5,
                 verbose=False):
    """
    Interpolate satellite data variable to a time/latitude grid for mapping.

    Separates satellite data into ascending and descending orbital passes at
    mid-latitudes, then interpolates the specified variable onto a regular
    time-latitude grid.

    Parameters
    ----------
    satDataDict : dict
        Dictionary containing satellite data with keys:
        - 'times' : array-like
            Datetime objects representing observation times
        - 'lats' : array-like
            Latitude values in degrees
        - **var** : array-like
            Variable data to interpolate (key name specified by var parameter)
        - 'orbital_period' : float
            Satellite orbital period in seconds
    var : str, optional
        Variable name to interpolate from satDataDict (default: 'sat_rho')
    lat_lim : float, optional
        Latitude limits for output grid in degrees. Set to the min/max observed by
        satellite if None (default: None)
    dlat : float, optional
        Latitude step for output grid in degrees (default: 2.5)
    verbose : bool, optional
        If True, print debug information (default: False)

    Returns
    -------
    dict
        Dictionary containing gridded data with keys:
        - 'desc' : numpy.ndarray
            Gridded data for descending orbital passes
        - 'asc' : numpy.ndarray
            Gridded data for ascending orbital passes
        - 'times' : numpy.ndarray
            Time coordinate array (hours from start)
        - 'lats' : numpy.ndarray
            Latitude coordinate array (degrees)
        - 'asc_n' : float
            Local time of ascending node in hours (0 if not found, orbit too short)
        - 'desc_n' : float
            Local time of descending node in hours (0 if not found)

    """

    # Make sure we have the orbital period
    if 'orbital_period' not in satDataDict.keys():
        if verbose:
            print("-> Griddata4map: No orbital period found. Calling calc_period()")
        satDataDict = satellite.calc_period(satDataDict, verbose=verbose)
        if 'orbital_period' not in satDataDict.keys():
            raise ValueError("No orbital period found. Cannot grid data for map.")

    # find where sat is descending & at midlatitudes

    if lat_lim is None:
        low_lat_lim = np.min(satDataDict['lats'])
        high_lat_lim = np.max(satDataDict['lats'])
        if verbose:
            print(f"-> Griddata4map: No lat_lim given. Setting to ({high_lat_lim:.1f},{low_lat_lim:.1f}) deg")
    else:
        low_lat_lim = -lat_lim
        high_lat_lim = lat_lim
        if verbose:
            print(f"-> Griddata4map: Setting lat_lim to +/- {lat_lim:.1f} deg")
            
    desc_indices = np.where((np.diff(satDataDict['lats']) < 0) & (
                            (satDataDict['lats'] < high_lat_lim)[1:]) & (
                            (satDataDict['lats'] > low_lat_lim)[1:]))[0]

    asc_indices = np.where((np.diff(satDataDict['lats']) > 0) & (
                           (satDataDict['lats'] < high_lat_lim)[1:]) & (
                           (satDataDict['lats'] > low_lat_lim)[1:]))[0]

    t0, t1 = min(satDataDict['times']), max(satDataDict['times'])

    Xs = np.arange(0, (t1 - t0).total_seconds()/3600.0, 
                   satDataDict['orbital_period'].total_seconds()/3600.0)
    Ys = np.arange(low_lat_lim, high_lat_lim+dlat, dlat)
    
    gridX, gridY = np.meshgrid(Xs, Ys, indexing='ij')

    ts = []
    for dt in satDataDict['times'] - t0:
        ts.append(dt.total_seconds()/3600.0)
    ts = np.array(ts)

    desc = scipy.interpolate.griddata(
        ((ts[desc_indices]),
         satDataDict['lats'][desc_indices]),
        satDataDict[var][desc_indices],
        (gridX, gridY),
        method='linear')

    asc = scipy.interpolate.griddata(
        ((ts[asc_indices]),
         satDataDict['lats'][asc_indices]),
        satDataDict[var][asc_indices],
        (gridX, gridY),
        method='linear')
    
    # get local time
    if 'lst' in satDataDict.keys():
        lts = satDataDict['lst']
    elif 'lsts' in satDataDict.keys():
        lts = satDataDict['lsts']
    else:
        lts = time_conversion.ut_to_lt(satDataDict['times'],
                                       satDataDict['lons'])
        
    # find where lat crosses 0 deg north going north (ascending node)
    an_idx = np.where((np.diff(satDataDict['lats']) > 0) 
                      & ((satDataDict['lats'][:-1] < 0)))[0]
    if len(an_idx) > 0:
        an = np.median(lts[an_idx])
    else:
        an = 0
        if verbose:
            print("-> Griddata4map: No ascending node found.")
    
    dn_idx = np.where((np.diff(satDataDict['lats']) < 0) 
                      & ((satDataDict['lats'][:-1] > 0)))[0]
    if len(dn_idx) > 0:
        dn = np.median(lts[dn_idx])
    else:
        dn = 0
        if verbose:
            print("-> Griddata4map: No descending node found.")

    return dict(desc =desc,
                asc = asc,
                times = Xs,
                lats = Ys,
                asc_n = an,
                desc_n = dn)
