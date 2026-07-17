#!/usr/bin/env python3

import numpy as np

from netCDF4 import Dataset

from pyitm.fileio import util
from pyitm.fileio import netcdfio
from pyitm.fileio import variables
from pyitm.modeldata import utils as mutils

# TIEGCM writes missing_value = 1e36 as a float64 attribute on float32 data,
# which netCDF4 refuses to apply - fills are handled by value instead.
cFillValue = 1e35

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def check_whether_tiegcm(filename):
    state = False
    with Dataset(filename, 'r') as ncfile:
        if (('lev' in ncfile.dimensions) and \
            ('ilev' in ncfile.dimensions) and \
            ('mlev' in ncfile.dimensions)):
            state = True
    return state

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def fills_to_nan(val):
    val = val.astype(float)
    val[val >= cFillValue] = np.nan
    return val

#------------------------------------------------------------------------------
# (lev, lat, lon) or (lat, lon) -> (lon, lat[, lev]), with the longitudes
# rolled from -180..180 to 0..360 order
#------------------------------------------------------------------------------

def orient(val, iRoll):
    return np.roll(val.T, -iRoll, axis = 0)

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def read_tiegcm_one_header(filename):
    """Read variable and time info from a TIEGCM secondary/primary history
    file. Only variables on the geographic grid are listed.
    """

    inventory = netcdfio.read_netcdf_inventory(filename)

    data = {'filename': filename}
    data['vars'] = []
    data['longname'] = []
    data['units'] = []
    for sig, names in inventory['groups'].items():
        if not (('lat' in sig) and ('lon' in sig) and
                set(sig) <= {'lev', 'ilev', 'lat', 'lon'}):
            continue
        for name in names:
            data['vars'].append(name)
            data['longname'].append(inventory['varinfo'][name]['longname'])
            data['units'].append(inventory['varinfo'][name]['units'])
    data['nvars'] = len(data['vars'])
    data['shortname'] = variables.get_short_names(data['vars'])

    with Dataset(filename, 'r') as ncfile:
        data['nlons'] = len(ncfile.dimensions['lon'])
        data['nlats'] = len(ncfile.dimensions['lat'])
        data['nalts'] = len(ncfile.dimensions['lev'])
    data['nblocks'] = 0
    data['times'] = inventory['times']

    return data

#------------------------------------------------------------------------------
# Read a series of vars / files, subset to start/stop/time at read time
#------------------------------------------------------------------------------

def read_tiegcm_all_files(filelist, varlist = [-1], verbose = False,
                          start = None, stop = None, time = None):

    filelist = util.any_to_filelist(filelist)

    inventory = netcdfio.read_netcdf_inventory(filelist, verbose = verbose)
    header = read_tiegcm_one_header(filelist[0])

    if varlist == [-1]:
        # reading everything spans multiple grids -> refused below with
        # a message listing the groups; vars must be named explicitly
        varlist = header['vars']
    nVars = len(varlist)
    netcdfio.check_signature_groups(inventory, varlist)

    sigs = set()
    for var in varlist:
        sig = tuple(d for d in inventory['varinfo'][var]['dims']
                    if d != 'time')
        if not (('lat' in sig) and ('lon' in sig) and
                set(sig) <= {'lev', 'ilev', 'lat', 'lon'}):
            raise ValueError(f"'{var}' is not on the geographic grid "
                             f"(dims {inventory['varinfo'][var]['dims']}), "
                             "cannot read it with this reader")
        sigs.add(sig)

    # midpoint (lev) or interface (ilev) grid; None when all vars are 2D
    grid = None
    for sig in sigs:
        if ('lev' in sig):
            grid = 'lev'
        elif ('ilev' in sig) and (grid is None):
            grid = 'ilev'

    if start is None and stop is None and time is None:
        iSelected = np.arange(inventory['ntimes'])
    else:
        iSelected = mutils.resolve_time_indices(inventory['times'],
                                                start, stop, time)
    nTimes = len(iSelected)
    allTimes = [inventory['times'][i] for i in iSelected]

    if verbose:
        est = netcdfio.estimate_read_size(inventory, varlist, nTimes)
        print(' -> Reading %d vars x %d times, ~%.1f MB' %
              (nVars, nTimes, est / 1e6))

    nLons = header['nlons']
    nLats = header['nlats']
    nAlts = header['nalts'] if grid is not None else 1

    zName = None
    if grid is not None:
        # altitude comes from geometric height ZG (cm); primary files only
        # carry geopotential Z, which is used as-is
        zName = 'ZG' if 'ZG' in inventory['varinfo'] else 'Z'
        if verbose and zName == 'Z':
            print(' -> no ZG in file, using geopotential Z for altitudes')

    out_shape = [nTimes]
    if nVars > 1:
        out_shape.append(nVars)
    out_shape.extend([nLons, nLats, nAlts])
    allData = np.zeros(out_shape)
    allAlts = np.zeros((nTimes, nLons, nLats, nAlts)) \
        if grid is not None else None

    ncfile = None
    iOpen = -1
    lons1d = None
    for iOut, i in enumerate(iSelected):
        iFile = inventory['ifile'][i]
        if iFile != iOpen:
            if ncfile is not None:
                ncfile.close()
            if verbose:
                print('-> Reading tiegcm : ', filelist[iFile],
                      ' --> Vars : ', varlist)
            ncfile = Dataset(filelist[iFile], 'r')
            # fills are caught by value in fills_to_nan; auto-masking only
            # half works on these files (missing_value is never applied)
            ncfile.set_auto_mask(False)
            iOpen = iFile
            if lons1d is None:
                lon = np.asarray(ncfile.variables['lon'][:])
                lats1d = np.asarray(ncfile.variables['lat'][:])
                iRoll = int(np.searchsorted(lon, 0.0))
                lons1d = np.roll(lon, -iRoll) % 360
        iLocal = inventory['ilocal'][i]
        for iVar, var in enumerate(varlist):
            val = orient(fills_to_nan(ncfile.variables[var][iLocal, ...]),
                         iRoll)
            # 2D vars fill all levels, like the generic netcdf reader
            if val.ndim == 2:
                val = val[..., np.newaxis]
            if (nVars == 1):
                allData[iOut, ...] = val
            else:
                allData[iOut, iVar, ...] = val
        if grid is not None:
            z = orient(fills_to_nan(ncfile.variables[zName][iLocal, ...]),
                       iRoll) / 1e5
            if grid == 'lev':
                # midpoint altitudes = adjacent-interface mean; the top
                # midpoint level is half-spacing extrapolated (its data is
                # all fill anyway)
                zm = np.empty_like(z)
                zm[..., :-1] = 0.5 * (z[..., :-1] + z[..., 1:])
                zm[..., -1] = z[..., -1] + 0.5 * (z[..., -1] - z[..., -2])
                z = zm
            allAlts[iOut] = z
    ncfile.close()

    # plots are on model levels: each level gets its mean altitude over
    # the selected times, and allalts keeps the full per-point values
    if grid is not None:
        alts1d = np.nanmean(allAlts, axis = (0, 1, 2))
        if np.any(np.isnan(alts1d)):
            raise ValueError("no valid altitudes at the selected times "
                             f"({zName} is all fill values there)")
    else:
        alts1d = np.array([0.0])
    lons, lats, alts = np.meshgrid(lons1d, lats1d, alts1d, indexing = 'ij')

    vars = []
    for var in varlist:
        vars.append(var)

    data = {'times': allTimes,
            'data': allData,
            'vars': vars,
            'shortname': variables.get_short_names(vars), \
            'longname': variables.get_long_names(vars), \
            'lons': lons,
            'lats': lats,
            'alts': alts,
            'ntimes': nTimes,
            'nvars': nVars,
            'nblocks': 0,
            'nlons': nLons,
            'nlats': nLats,
            'nalts': nAlts}
    if allAlts is not None:
        data['allalts'] = allAlts

    return data
