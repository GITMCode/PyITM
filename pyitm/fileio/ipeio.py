#!/usr/bin/env python3

import os
from datetime import datetime
from struct import unpack
import numpy as np
from pyitm.fileio import util
from pyitm.fileio import variables
from pyitm.general import time_conversion as tc

from netCDF4 import Dataset

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

def check_whether_ipe(filename):
    state = False
    with Dataset(filename, 'r') as ncfile:
        # Process header information: nlons, nlats, nalts, nblocks
        if (('x01' in ncfile.dimensions) and \
            ('x02' in ncfile.dimensions) and \
            ('x03' in ncfile.dimensions)):
            state = True
        if (('phony_dim_0' in ncfile.dimensions) and \
            ('phony_dim_1' in ncfile.dimensions) and \
            ('phony_dim_2' in ncfile.dimensions)):
            state = True
    ncfile.close
    return state

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

def is_grid_file(filename):
    state = False
    with Dataset(filename, 'r') as ncfile:
        if (('phony_dim_0' in ncfile.dimensions) and \
            ('phony_dim_1' in ncfile.dimensions) and \
            ('phony_dim_2' in ncfile.dimensions) and \
            ('phony_dim_3' in ncfile.dimensions)):
            state = True
    ncfile.close
    return state

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

def get_date_from_filename(filename):
    year = int(filename[-15:-11])
    month = int(filename[-11:-9])
    day = int(filename[-9:-7])
    hour = int(filename[-7:-5])
    minute = int(filename[-5:-3])
    return datetime(year, month, day, hour, minute, 0)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

def read_ipe_grid_header(filename, verbose = False):
    """Read all keys and such from netcdf file

    Parameters
    ----------
    filename : str
        An IPE GRID netCDF filename
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
        data['nlons'] = len(ncfile.dimensions['phony_dim_2'])
        data['nlats'] = len(ncfile.dimensions['phony_dim_0'])
        data['nalts'] = len(ncfile.dimensions['phony_dim_1'])
        data['nblocks'] = 0

        # Included for compatibility
        data['vars'] = []
        data['longname'] = []
        data['units'] = []
        for key in ncfile.variables.keys():
            # Only store variables that have at least 2 dimensions (exclude time!)
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
        
        data['isEnsemble'] = False

    data['shortname'] = variables.get_short_names(data['vars'])
    data['times'] = datetime(1965,1,1,0,0,0)

    print(data['vars'])
    return data

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

def read_ipe_one_header(filename):
    """Read all keys and such from netcdf file

    Parameters
    ----------phony_dim_0
    filename : str
        An IPE netCDF filename
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
        raise IOError(f"unknown IPE netCDF blocked file: {filename}")

    # NOTE: Includes header information for easy access until
    #       updated package structure is confirmed
    # Initialize data dict with defaults (will remove these defaults later)
    data = {'filename': filename}

    with Dataset(filename, 'r') as ncfile:
        # Process header information: nlons, nlats, nalts, nblocks
        if ('x01' in ncfile.dimensions):
            data['nlons'] = len(ncfile.dimensions['x01'])
            data['nlats'] = len(ncfile.dimensions['x02'])
            data['nalts'] = len(ncfile.dimensions['x03'])
        else:
            data['nlons'] = len(ncfile.dimensions['phony_dim_0'])
            data['nlats'] = len(ncfile.dimensions['phony_dim_1'])
            data['nalts'] = len(ncfile.dimensions['phony_dim_2'])

        data['nblocks'] = 0

        # Included for compatibility
        data['vars'] = []
        data['longname'] = []
        data['units'] = []
        for key in ncfile.variables.keys():
            # Only store variables that have at least 2 dimensions (exclude time!)
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
        
        data['times'] = get_date_from_filename(filename)
        data['isEnsemble'] = False

    data['shortname'] = variables.get_short_names(data['vars'])
    data['longname'] = variables.get_long_names(data['vars'])
    
    return data


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

def reshape_ipe_array(varIn, ipeGridShape, verbose = False):

    iNt_ = ipeGridShape['iNorthTop']
    iSt_ = ipeGridShape['iSouthTop']-1
    iM_ = ipeGridShape['iMax']

    nLons = ipeGridShape['nLons']
    nLats = ipeGridShape['nLats'] * 2
    nAlts = ipeGridShape['nAlts']

    if (np.max(iNt_) < nAlts):
        nAlts = np.max(iNt_)

    if (nLons > 0):
        varOut = np.empty((nLons, nLats, nAlts))
    else:
        varOut = np.empty((nLats, nAlts))

    varOut.fill(np.nan)

    for i_ in range(ipeGridShape['nLats']):
        if (nLons > 0):
            # grab south, put it in first, but reverse it, since it is high alt first:
            sub = varIn[:, i_, iSt_[i_]:iM_[i_]]
            n = len(sub[0,:])
            varOut[:, i_, 0:n] = np.flip(sub, axis = 1)
            if (n < nAlts):
                for j in range(n,nAlts):
                    varOut[:, i_, j] = varOut[:, i_, n-1] 
            # grab north, put it in the opposite side of the array:
            sub = varIn[:, i_, 0:iNt_[i_]]
            n = len(sub[0,:])
            varOut[:, -(i_+1), 0:n] = sub
            if (n < nAlts):
                for j in range(n,nAlts):
                    varOut[:, -(i_+1), j] = varOut[:, -(i_+1), n-1] 
        else:
            # grab south, put it in first, but reverse it, since it is high alt first:
            sub = varIn[i_, iSt_[i_]:iM_[i_]]
            n = len(sub)
            varOut[i_, 0:n] = np.flip(sub)
            if (verbose):
                print('south : ', n, varOut[i_, n-1])
            if (n < nAlts):
                for j in range(n, nAlts):
                    varOut[i_, j] = varOut[i_, n-1] 
            # grab north, put it in the opposite side of the array:
            sub = varIn[i_, 0:iNt_[i_]]
            n = len(sub)
            varOut[-(i_+1), 0:n] = sub
            if (verbose):
                print('north : ', n, varOut[-(i_+1), n-1])
            if (n < nAlts):
                for j in range(n,nAlts):
                    varOut[-(i_+1), j] = varOut[-(i_+1), n-1] 

    return varOut

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

def read_ipe_grid_file(filename, verbose = False):

    # Checks for file existence
    if not os.path.isfile(filename):
        raise IOError(f"unknown IPE netCDF GRID file: {filename}")

    # NOTE: Includes header information for easy access until
    #       updated package structure is confirmed
    # Initialize data dict with defaults (will remove these defaults later)

    header = read_ipe_grid_header(filename, verbose = verbose)
    if (verbose):
        print('Variables from IPE grid (header) file: ', header['vars'])

    with Dataset(filename, 'r') as ncfile:
        alts2draw = np.array(ncfile.variables['altitude']) / 1000.0
        geolons3draw = np.array(ncfile.variables['longitude']) * 180.0 / np.pi
        lats2draw = 90.0 - np.array(ncfile.variables['m_colat']) * 180.0 / np.pi
        geolats3draw = 90.0 - np.array(ncfile.variables['colatitude']) * 180.0 / np.pi
        iNorthTop_ = np.array(ncfile.variables['northern_top']).astype(np.int_)
        iSouthTop_ = np.array(ncfile.variables['southern_top']).astype(np.int_)
        iMax_ = np.array(ncfile.variables['tube_max']).astype(np.int_)
        ipeGridShape = {
            'iNorthTop': iNorthTop_,
            'iSouthTop': iSouthTop_,
            'iMax': iMax_
        }
        nLons, nLats, nAlts = np.shape(geolons3draw)
        ipeGridShape['nLons'] = 0
        ipeGridShape['nLats'] = nLats
        ipeGridShape['nAlts'] = nAlts

        alts2d = reshape_ipe_array(alts2draw, ipeGridShape, verbose = verbose)
        lats2d = reshape_ipe_array(lats2draw, ipeGridShape)
        ipeGridShape['nLons'] = nLons
        geolons3d = reshape_ipe_array(geolons3draw, ipeGridShape)
        geolats3d = reshape_ipe_array(geolats3draw, ipeGridShape)

    nLons, nLats, nAlts = np.shape(geolons3d)
    # lons is actually the geographic longitude instead of the magnetic
    # longitude.  Need to make up the magnetic longitude:
    dLon = 360.0/(nLons)
    lons1d = np.arange(0, 360, dLon)
    lons3d = np.zeros((nLons, nLats, nAlts))
    for iLat in range(nLats):
        for iAlt in range(nAlts):
            lons3d[:,iLat,iAlt] = lons1d
    lats3d = np.zeros((nLons, nLats, nAlts))
    alts3d = np.zeros((nLons, nLats, nAlts))
    for iLon in range(nLons):
        lats3d[iLon, :, :] = lats2d
        alts3d[iLon, :, :] = alts2d

    vars = ['lons', 'lats', 'alts']
    data = {'filename': filename,
            'nblocks': 0,
            'times': datetime(1965,1,1,0,0,0),
            'vars': vars,
            'shortname': variables.get_short_names(vars), \
            'longname': variables.get_long_names(vars), \
            'lon': lons3d,
            'geolon': geolons3d,
            'geolat': geolats3d,
            'lat': lats3d,
            'z': alts3d,
            'ipeGridShape': ipeGridShape,
            'nlons' : nLons,
            'nlats': nLats,
            'nalts': nAlts}

    return data

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def read_ipe_one_file(filename, \
                      file_vars = None, \
                      verbose = False):
    """Read all data from an IPE netcdf file.

    Parameters
    ----------
    filename : str
        An IPE netCDF filename
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
    This routine only works with blocked Aether netCDF files.

    """

    # Checks for file existence
    if not os.path.isfile(filename):
        raise IOError(f"unknown IPE netCDF blocked file: {filename}")

    # NOTE: Includes header information for easy access until
    #       updated package structure is confirmed
    # Initialize data dict with defaults (will remove these defaults later)
    data = {'filename': filename,
            'units': {},
            'long_name': None}
    
    if verbose:
        print('-> Reading IPE netcdf : ', filename, ' --> Vars : ', file_vars)

    with Dataset(filename, 'r') as ncfile:
        # Process header information: nlons, nlats, nalts, nblocks
        if ('x01' in ncfile.dimensions):
            data['nlons'] = len(ncfile.dimensions['x01'])
            data['nlats'] = len(ncfile.dimensions['x02'])
            data['nalts'] = len(ncfile.dimensions['x03'])
        else:
            data['nlons'] = len(ncfile.dimensions['phony_dim_0'])
            data['nlats'] = len(ncfile.dimensions['phony_dim_1'])
            data['nalts'] = len(ncfile.dimensions['phony_dim_2'])
        data['nblocks'] = 0
            
        # Included for compatibility
        data['vars'] = [var for var in ncfile.variables.keys()
                        if file_vars is None or var in file_vars]

        # Fetch requested variable data
        for varName in data['vars']:
            data[varName] = np.array(ncfile.variables[varName])
            # IPE doesn't offer units
            data['units'][varName] = ''

        data['isEnsemble'] = False

    data['times'] = get_date_from_filename(filename)

    return data

# -----------------------------------------------------------------------------
# This reads in a series of vars / files and returns the 3D information
#-----------------------------------------------------------------------------

def read_ipe_all_files(filelist, varlist=[-1], verbose=False):

    filelist = util.any_to_filelist(filelist)

    # Get the prefixes for all entries in filelist;
    prefixes = np.unique([file.split('/')[-1].split('_')[0] \
                          for file in filelist])
    if len(prefixes) > 1: # make sure there is only one output type.
        raise ValueError("Multiple output types cannot be read by this " +
                         "function.\n\tProvided: " + str(prefixes))

    # first read in spatial information:
    gridfile = 'IPE_Grid.nc'
    spatialData = read_ipe_grid_file(gridfile, verbose=False)
    ipeGridShape = spatialData['ipeGridShape']

    nTimes = len(filelist)
    if varlist != [-1]:
       nVars = len(varlist)
    else: # varlist=[-1] means we read in all variables
        header = read_ipe_one_header(filelist[0], verbose=False)
        varlist = header['vars']
        nVars = len(varlist)

    allTimes = []
    # This assumes we have 3D arrays for the coord info.
    lons = spatialData['lon']
    lats = spatialData['lat']
    alts = spatialData['z']
    nLons = len(lons[:, 0, 0])
    nLats = len(lats[0, :, 0])
    nAlts = len(alts[0, 0, :])
    nBlocks = 0
        
    if (nVars == 1):
        allData = np.zeros((nTimes, nLons, nLats, nAlts))
    else:
        allData = np.zeros((nTimes, nVars, nLons, nLats, nAlts))

    for iTime, filename in enumerate(filelist):
        data = read_ipe_one_file(filename, varlist, verbose=verbose)
        allTimes.append(data["times"])
        for iVar, var in enumerate(varlist):
            if (nVars == 1):
                allData[iTime, :, :, :] = reshape_ipe_array(data[var][:, :, :], ipeGridShape)
            else:
                allData[iTime, iVar, :, :, :] = reshape_ipe_array(data[var][:, :, :], ipeGridShape)
                
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
            'nblocks' : nBlocks,
            'nlons' : nLons,
            'nlats': nLats,
            'nalts': nAlts}
    
    return data



