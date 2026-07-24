#!/usr/bin/env python3

import os
from datetime import datetime
import numpy as np
from pyitm.fileio import util
from pyitm.fileio import variables
from pyitm.general import time_conversion as tc

from netCDF4 import Dataset

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def check_whether_wam(filename):
    state = False
    with Dataset(filename, 'r') as ncfile:
        # Process header information: nlons, nlats, nalts, nblocks
        if (('lat' in ncfile.dimensions) and \
            ('lon' in ncfile.dimensions) and \
            ('hlevs' in ncfile.dimensions) and \
            ('time' in ncfile.dimensions)):
            state = True
        if (('x01' in ncfile.dimensions) and \
            ('x02' in ncfile.dimensions) and \
            ('x03' in ncfile.dimensions)):
            state = True
    ncfile.close
    return state

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def read_wam_one_file(filename, \
                      file_vars = None, \
                      verbose = False):
    """Read all data from a WAM netcdf file.

    Parameters
    ----------
    filename : str
        An WAM netCDF filename
    file_vars : list or NoneType
        List of desired variable neames to read, or None to read all
        (default=None)

    Returns
    -------
    data : dict
        A dictionary containing all data from the netCDF file, including:
        filename - filename of file containing header data
        nlons - number of longitude grids per block
        nlats - number of latitude grids per block
        nalts - number of altitude grids per block
        nblocks - number of blocks in file
        vars - list of data variable names
        time - datetime for time of file
        isEnsemble - if true, stores ensembleNumber and ensembleMembers
        The dictionary also contains a read_routines.DataArray keyed to the
        corresponding variable name. Each DataArray carries both the variable's
        data from the netCDF file and the variable's corresponding attributes.

    Raises
    --------
    IOError
        If the input file does not exist
    KeyError
        If any expected dimensions of the input netCDF file are not present

    Notes
    -----
    This routine only works with WAM netCDF files.

    """

    # Checks for file existence
    if not os.path.isfile(filename):
        raise IOError(f"unknown WAM netCDF blocked file: {filename}")

    # NOTE: Includes header information for easy access until
    #       updated package structure is confirmed
    # Initialize data dict with defaults (will remove these defaults later)
    data = {'filename': filename,
            'units': {},
            'long_name': None,
            'isAltVarying': False}
    
    if verbose:
        print(' -> Reading WAM netcdf: ', filename, ' --> Vars : ', file_vars)
    else:
        print(' -> Reading WAM netcdf: ', filename)
        
    with Dataset(filename, 'r') as ncfile:
        # Process header information: nlons, nlats, nalts, nblocks
        if ('lon' in ncfile.dimensions):
            data['nlons'] = len(ncfile.dimensions['lon'])+1
        else:
            if ('x01' in ncfile.dimensions):
                data['nlons'] = len(ncfile.dimensions['x01'])+1
            else:
                data['nlons'] = 0
        if ('lat' in ncfile.dimensions):
            data['nlats'] = len(ncfile.dimensions['lat'])
        else:
            if ('x02' in ncfile.dimensions):
                data['nlats'] = len(ncfile.dimensions['x02'])
            else:
                data['nlats'] = 0
        if ('hlevs' in ncfile.dimensions):
            data['nalts'] = len(ncfile.dimensions['hlevs'])
        else:
            if ('x03' in ncfile.dimensions):
                data['nalts'] = len(ncfile.dimensions['x03'])
            else:
                data['nalts'] = 0
        data['nblocks'] = 0

        if (file_vars):
            if ('height' in ncfile.variables.keys()):
                # We are going to assume that the height is varying here!
                data['isAltVarying'] = True
                if (verbose):
                    print('  --> WAM file has height as a variable!')
            if ('hlevs' in file_vars):
                if (verbose):
                    print('  --> Requesting hlevs, subing height')
                file_vars[file_vars.index('hlevs')] = 'height'
        
        # Get all of the variable names:
        allVars = [var for var in ncfile.variables.keys()
                   if file_vars is None or var in file_vars]
        nLons = data['nlons']
        nLats = data['nlats']
        nAlts = data['nalts']

        data['vars'] = []
        # Fetch requested variable data

        for varName in allVars:
            temp = np.array(ncfile.variables[varName])
            # we really only want the variables that are 4d:
            isFound = False
            if (len(np.shape(temp)) == 4):
                # Some WAM files have time as an index...
                data['vars'].append(varName)
                temp3d = np.zeros((nLons, nLats, nAlts))
                for iAlt in range(nAlts):
                    for iLat in range(nLats):
                        temp3d[:-1, iLat, iAlt] = temp[0, iAlt, iLat, :]
                        temp3d[-1, iLat, iAlt] = temp3d[0, iLat, iAlt]
                data[varName] = temp3d
                data['units'][varName] = ncfile.variables[varName].units
                isFound = True
            if ((len(np.shape(temp)) == 3) and (varName != 'height')):
                # Some don't have time as an index
                data['vars'].append(varName)
                temp3d = np.zeros((nLons, nLats, nAlts))
                for iAlt in range(nAlts):
                    for iLat in range(nLats):
                        temp3d[:-1, iLat, iAlt] = temp[iAlt, iLat, :]
                        temp3d[-1, iLat, iAlt] = temp3d[0, iLat, iAlt]
                # it seems like the poles have bad values, so lets overwrite them:
                for iAlt in range(nAlts):
                    temp3d[:, 0, iAlt] = np.mean(temp3d[:, 1, iAlt])
                    temp3d[:, -1, iAlt] = np.mean(temp3d[:, -2, iAlt])
                data[varName] = temp3d
                data['units'][varName] = ncfile.variables[varName].units
                isFound = True
            if (not isFound):
                # want to include other variables:
                if (varName == 'lon'):
                    lons3d = np.zeros((nLons, nLats, nAlts))
                    for iAlt in range(nAlts):
                        for iLat in range(nLats):
                            lons3d[:-1, iLat, iAlt] = temp
                            lons3d[-1, iLat, iAlt] = temp[0] + 360.0
                    data['Longitude'] = lons3d
                if (varName == 'lat'):
                    lats3d = np.zeros((nLons, nLats, nAlts))
                    for iAlt in range(nAlts):
                        for iLon in range(nLons):
                            lats3d[iLon, :, iAlt] = temp
                    data['Latitude'] = lats3d
                if (varName == 'hlevs'):
                    alts3d = np.zeros((nLons, nLats, nAlts))
                    for iLat in range(nLats):
                        for iLon in range(nLons):
                            alts3d[iLon, iLat, :] = temp
                    data['Altitude'] = alts3d
                if (varName == 'height'):
                    # The heights are 3D variables in some WAM files.
                    alts3d = np.zeros((nLons, nLats, nAlts))
                    for iLat in range(nLats):
                        for iLon in range(nLons-1):
                            alts3d[iLon, iLat, :] = temp[:, iLat, iLon]
                    alts3d[alts3d < 0] = 0.0
                    alts3d[-1, :, :] = alts3d[0, :, :]
                    # it seems like the poles have bad altitudes,
                    # so lets overwrite them:
                    for iAlt in range(nAlts):
                        alts3d[:, 0, iAlt] = np.mean(alts3d[:, 1, iAlt])
                        alts3d[:, -1, iAlt] = np.mean(alts3d[:, -2, iAlt])
                    
                    data['Altitude'] = alts3d/1000.0
                    data['height'] = alts3d/1000.0
                    data['vars'].append('height')
                    
        if ('time' in ncfile.variables.keys()):
            temp = np.array(ncfile.variables['time'])
            time = tc.epoch_to_datetime(temp[0]*86400.0, \
                                        t0 = datetime(1970, 1, 1))
        else:
            time = tc.convert_string_to_datetime(ncfile.fcst_date, \
                                                 format = 'ymd_hms')
        data['times'] = time
        data['isEnsemble'] = False

    return data


def read_wam_one_header(filename):
    """Read all keys and such from netcdf file

    Parameters
    ----------
    filename : str
        An WAM netCDF filename
    file_vars : list or NoneType
        List of desired variable neames to read, or None to read all
        (default=None)

    Returns
    -------
    data : dict
        A dictionary containing all data from the netCDF file, including:
        filename - filename of file containing header data
        nlons - number of longitude grids per block
        nlats - number of latitude grids per block
        nalts - number of altitude grids per block
        nblocks - number of blocks in file
        vars - list of data variable names
        times - datetime for time of file
        isEnsemble - if true, stores ensembleNumber and ensembleMembers
        The dictionary also contains a read_routines.DataArray keyed to the
        corresponding variable name. Each DataArray carries both the variable's
        data from the netCDF file and the variable's corresponding attributes.

    Raises
    --------
    IOError
        If the input file does not exist
    KeyError
        If any expected dimensions of the input netCDF file are not present

    Notes
    -----
    This routine only works with blocked Aether netCDF files.

    """

    # Checks for file existence
    if not os.path.isfile(filename):
        raise IOError(f"unknown aether netCDF blocked file: {filename}")

    # NOTE: Includes header information for easy access until
    #       updated package structure is confirmed
    # Initialize data dict with defaults (will remove these defaults later)
    data = {'filename': filename}

    with Dataset(filename, 'r') as ncfile:
        # Process header information: nlons, nlats, nalts, nblocks
        if ('lon' in ncfile.dimensions):
            data['nlons'] = len(ncfile.dimensions['lon'])
            data['nlats'] = len(ncfile.dimensions['lat'])
            data['nalts'] = len(ncfile.dimensions['hlevs'])
        else:
            data['nlons'] = len(ncfile.dimensions['x01'])
            data['nlats'] = len(ncfile.dimensions['x02'])
            data['nalts'] = len(ncfile.dimensions['x03'])
        data['nblocks'] = 0

        # Included for compatibility
        data['vars'] = []
        data['longname'] = []
        data['units'] = []
        for key in ncfile.variables.keys():
            # Only store variables that have at least 2 dimensions
            if (len(ncfile.variables[key].shape) > 1): 
                data['vars'].append(key)
                if hasattr(ncfile.variables[key], 'units'):
                    data['units'].append(getattr(ncfile.variables[key], 'units'))
                else:
                    data['units'].append('')
                if hasattr(ncfile.variables[key], 'long_name'):
                    data['longname'].append(getattr(ncfile.variables[key],
                                                'long_name'))
                else:
                    data['longname'].append(key)

        if ('time' in ncfile.variables.keys()):
            temp = np.array(ncfile.variables['time'])
            time = tc.epoch_to_datetime(temp[0]*86400.0, \
                                        t0 = datetime(1970, 1, 1))
        else:
            time = tc.convert_string_to_datetime(ncfile.fcst_date, format = 'ymd_hms')
        data['times'] = time
        data['isEnsemble'] = False

    data['shortname'] = variables.get_short_names(data['vars'])
    
    return data


# -----------------------------------------------------------------------------
# This reads in a series of vars / files and returns the 3D information
#-----------------------------------------------------------------------------

def read_wam_all_files(filelist, varlist = [-1], verbose = False):

    filelist = util.any_to_filelist(filelist)

    # Get the prefixes for all entries in filelist;
    prefixes = np.unique([file.split('/')[-1].split('_')[0] \
                          for file in filelist])
    if len(prefixes) > 1: # make sure there is only one output type.
        raise ValueError("Multiple output types cannot be read by this " +
                         "function.\n\tProvided: " + str(prefixes))

    # first read in spatial information:
    vars = ['lon', 'lat', 'hlevs']
    spatialData = read_wam_one_file(filelist[0], vars, verbose=False)
    if (verbose):
        print('  -> Spatial Data Alt Varying : ', spatialData['isAltVarying'])

    nTimes = len(filelist)
    if (varlist != [-1]):
        nVars = len(varlist)
    else: # varlist=[-1] means we read in all variables
        header = read_wam_one_header(filelist[0], verbose=False)
        varlist = header['vars']
        nVars = len(varlist)

    allTimes = []
        
    if (spatialData['nblocks'] == 0):
        # This assumes we have 3D arrays for the coord info.
        lons = spatialData['Longitude']
        lats = spatialData['Latitude']
        alts = spatialData['Altitude']  # Convert from m to km
        nDims = len(np.shape(lons))
        nLons = len(lons[:, 0, 0])
        nLats = len(lats[0, :, 0])
        nAlts = len(alts[0, 0, :])
        nBlocks = 0
        
        if (nVars == 1):
            allData = np.zeros((nTimes, nLons, nLats, nAlts))
        else:
            allData = np.zeros((nTimes, nVars, nLons, nLats, nAlts))
        allAlts = np.zeros((nTimes, nLons, nLats, nAlts))

    for iTime, filename in enumerate(filelist):
        data = read_wam_one_file(filename, varlist, verbose=verbose)
        allTimes.append(data["times"])
        if (spatialData['isAltVarying']):
            altData = read_wam_one_file(filename, ['height'], verbose=verbose)
            allAlts[iTime, :, :, :] = altData['Altitude'][:, :, :]
        for iVar, var in enumerate(varlist):
            if (nVars == 1):
                allData[iTime, :, :, :] = data[var][:, :, :]
            else:
                allData[iTime, iVar, :, :, :] = data[var][:, :, :]
                
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
            'isAltVarying': spatialData['isAltVarying'],
            'ntimes': nTimes,
            'nvars': nVars,
            'nblocks' : nBlocks,
            'nlons' : nLons,
            'nlats': nLats,
            'nalts': nAlts}
    if (spatialData['isAltVarying']):
        data['allalts'] = allAlts
    
    return data

