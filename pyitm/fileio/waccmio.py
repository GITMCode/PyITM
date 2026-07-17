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

def check_whether_waccm(filename):
    state = False
    with Dataset(filename, 'r') as ncfile:
        # Process header information: nlons, nlats, nalts, nblocks
        if (('lat' in ncfile.dimensions) and \
            ('lon' in ncfile.dimensions) and \
            ('lev' in ncfile.dimensions) and \
            ('time' in ncfile.dimensions)):
            state = True
    ncfile.close
    return state

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def read_waccm_one_file(filename, \
                        file_vars = None, \
                        iStart = 0, \
                        iEnd = 0, \
                        iStep = 1, \
                        verbose = False):
    """Read all data from a WACCM netcdf file.

    Parameters
    ----------
    filename : str
        An WACCM netCDF filename
    file_vars : list or NoneType
        List of desired variable neames to read, or None to read all
        (default=None)
    iStart : starting time to read
    iEnd : ending time to read

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
    This routine only works with WACCM netCDF files.

    """

    # Checks for file existence
    if not os.path.isfile(filename):
        raise IOError(f"unknown WACCM netCDF blocked file: {filename}")

    # NOTE: Includes header information for easy access until
    #       updated package structure is confirmed
    # Initialize data dict with defaults (will remove these defaults later)
    data = {'filename': filename,
            'units': {},
            'long_name': None,
            'isAltVarying': False}
    
    if verbose:
        print(' -> Reading WACCM netcdf: ', filename, ' --> Vars : ', file_vars)
    else:
        print(' -> Reading WACCM netcdf: ', filename)
        
    with Dataset(filename, 'r') as ncfile:
        # Process header information: nlons, nlats, nalts, nblocks
        if ('lon' in ncfile.dimensions):
            data['nlons'] = len(ncfile.dimensions['lon'])+1
        else:
            data['nlons'] = 0
        if ('lat' in ncfile.dimensions):
            data['nlats'] = len(ncfile.dimensions['lat'])
        else:
            data['nlats'] = 0
        if ('lev' in ncfile.dimensions):
            data['nalts'] = len(ncfile.dimensions['lev'])
        else:
            data['nalts'] = 0
        data['nblocks'] = 0

        if (file_vars):
            if ('Z3GM' in ncfile.variables.keys()):
                # We are going to assume that the height is varying here!
                data['isAltVarying'] = True
                if (verbose):
                    print('  --> WACCM file has Z3GM as a variable!')
        
        # Get all of the variable names:
        allVars = [var for var in ncfile.variables.keys()
                   if file_vars is None or var in file_vars]
        nLons = data['nlons']
        nLats = data['nlats']
        nAlts = data['nalts']

        data['vars'] = []

        iTimes = list(range(iStart, iEnd+1, iStep))
        nTimes = len(iTimes)

        for varName in allVars:
            isFound = False
            if (varName == 'lon'):
                temp = np.array(ncfile.variables[varName])
                isFound = True
                lons3d = np.zeros((nLons, nLats, nAlts))
                for iAlt in range(nAlts):
                    for iLat in range(nLats):
                        lons3d[:-1, iLat, iAlt] = temp
                        lons3d[-1, iLat, iAlt] = temp[0] + 360.0
                    data['Longitude'] = lons3d
            if (varName == 'lat'):
                temp = np.array(ncfile.variables[varName])
                isFound = True
                lats3d = np.zeros((nLons, nLats, nAlts))
                for iAlt in range(nAlts):
                    for iLon in range(nLons):
                        lats3d[iLon, :, iAlt] = temp
                data['Latitude'] = lats3d
            if (varName == 'lev'):
                temp = np.array(ncfile.variables[varName])
                isFound = True
                pressure3d = np.zeros((nLons, nLats, nAlts))
                for iLat in range(nLats):
                    for iLon in range(nLons):
                        pressure3d[iLon, iLat, :] = temp[::-1]
                data['Pressure'] = pressure3d * 100.0
            if (not isFound):
                dims = ncfile.variables[varName].shape
                nDims = len(dims)
                
                # we really only want the variables that are 4d:
                if (nDims == 4):
                    if (verbose):
                        print('  -> Reading in var : ', varName)
                    # Some WACCM files have time as an index...
                    temp = np.array(ncfile.variables[varName][iTimes,:,:,:])
                    temp4d = np.zeros((nTimes, nLons, nLats, nAlts))
                    for iTime in range(nTimes):
                        for iAlt in range(nAlts):
                            for iLat in range(nLats):
                                # WACCM data starts at the top:
                                temp4d[iTime, :-1, iLat, iAlt] = \
                                    temp[iTime, nAlts -1 - iAlt, iLat, :]
                                temp4d[iTime, -1, iLat, iAlt] = \
                                    temp4d[iTime, 0, iLat, iAlt]
                    if (varName == 'Z3GM'):
                        data['Altitude'] = temp4d / 1000.0
                        data[varName] = temp4d
                        data['vars'].append(varName)
                    else:
                        data[varName] = temp4d
                        data['vars'].append(varName)
                    data['units'][varName] = ncfile.variables[varName].units
                    isFound = True
                    #for iAlt in range(nAlts):
                    #    temp3d[:, 0, iAlt] = np.mean(temp3d[:, 1, iAlt])
                    #    temp3d[:, -1, iAlt] = np.mean(temp3d[:, -2, iAlt])
                    
        if ('time' in ncfile.variables.keys()):
            temp = np.array(ncfile.variables['time'][iTimes])
            timestamp = ncfile.variables['time'].units[11:] 
            t0 = tc.convert_string_to_datetime(timestamp)
            time = tc.epoch_to_datetime(temp*86400.0, t0 = t0)
        data['times'] = time
        data['isEnsemble'] = False
        
    return data

def read_waccm_one_header(filename, verbose = False):
    """Read all keys and such from netcdf file

    Parameters
    ----------
    filename : str
        An WACCM netCDF filename
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
        data['nlons'] = len(ncfile.dimensions['lon'])
        data['nlats'] = len(ncfile.dimensions['lat'])
        data['nalts'] = len(ncfile.dimensions['lev'])
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
            timestamp = ncfile.variables['time'].units[11:] 
            t0 = tc.convert_string_to_datetime(timestamp)
            time = tc.epoch_to_datetime(temp*86400.0, t0 = t0)
        data['times'] = time
        data['isEnsemble'] = False

    data['shortname'] = variables.get_short_names(data['vars'])
    
    return data


# -----------------------------------------------------------------------------
# This reads in a series of vars / files and returns the 3D information
#-----------------------------------------------------------------------------

def read_waccm_all_files(filelist, varlist = [-1], verbose = False, \
                         iStart = 0, iEnd = 0, iStep = 1):

    filelist = util.any_to_filelist(filelist)

    # Get the prefixes for all entries in filelist;
    prefixes = np.unique([file.split('/')[-1].split('_')[0] \
                          for file in filelist])
    if len(prefixes) > 1: # make sure there is only one output type.
        raise ValueError("Multiple output types cannot be read by this " +
                         "function.\n\tProvided: " + str(prefixes))

    # first read in spatial information:
    cAltName = 'Z3GM'
    vars = ['lon', 'lat', cAltName]
    spatialData = read_waccm_one_file(filelist[0], vars, verbose = verbose)
    if (verbose):
        print('  -> Spatial Data Alt Varying : ', spatialData['isAltVarying'])
    nTimes = len(filelist)
    if (varlist != [-1]):
        nVars = len(varlist)
    else: # varlist=[-1] means we read in all variables
        header = read_waccm_one_header(filelist[0], verbose = verbose)
        varlist = header['vars']
        nVars = len(varlist)

    allTimes = []
    
    if (iStart < 0):
        iStart = 0 
    if (iEnd < 0):
        header = read_waccm_one_header(filelist[0], verbose = verbose)
        iEnd = len(header['times'])-1
    if (iStep < 0):        
        iStep = 1 
    iTimes = list(range(iStart, iEnd+1, iStep))
    nTimes = len(iTimes)
    
    # This assumes we have 3D arrays for the coord info.
    lons = spatialData['Longitude']
    lats = spatialData['Latitude']
    alts = spatialData['Altitude']  # Convert from m to km
    # It is hard coded to return pressure when asking for lev:
    vars = ['lev']
    pressureData = read_waccm_one_file(filelist[0], vars, verbose = verbose)
    pressure = pressureData['Pressure']
    nDims = len(np.shape(lons))
    nLons = len(lons[:, 0, 0])
    nLats = len(lats[0, :, 0])
    nAlts = spatialData['nalts']
    nBlocks = 0
        
    if (nVars == 1):
        allData = np.zeros((nTimes, nLons, nLats, nAlts))
    else:
        allData = np.zeros((nTimes, nVars, nLons, nLats, nAlts))
    allAlts = np.zeros((nTimes, nLons, nLats, nAlts))

    if (len(filelist) > 1):
        print('WACCM reader is not smart enough to read in multiple files')
        print(' -> only reading in first file!')
    filename = filelist[0]
    data = read_waccm_one_file(filename, varlist, verbose = verbose,
                               iStart = iStart, \
                               iEnd = iEnd, \
                               iStep = iStep)
    allTimes = data["times"]
    if (spatialData['isAltVarying']):
        altData = read_waccm_one_file(filename, [cAltName], \
                                      verbose=verbose,
                                      iStart = iStart, \
                                      iEnd = iEnd, \
                                      iStep = iStep)
        allAlts = altData['Altitude']
    for iVar, var in enumerate(varlist):
        if (nVars == 1):
            allData = data[var]
        else:
            for iTime in range(nTimes):
                allData[iTime, iVar, :, :, :] = data[var][iTime, :, :, :]
                
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
            'allalts': allAlts,
            'pressure': pressure,
            'isAltVarying': spatialData['isAltVarying'],
            'ntimes': nTimes,
            'nvars': nVars,
            'nblocks' : nBlocks,
            'nlons' : nLons,
            'nlats': nLats,
            'nalts': nAlts}
    
    return data

