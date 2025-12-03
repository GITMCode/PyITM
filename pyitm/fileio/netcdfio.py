#!/usr/bin/env python3

import os
from datetime import datetime
from struct import unpack
import numpy as np
from pyitm.fileio import util
from pyitm.fileio import variables
from pyitm.general import time_conversion as tc

from netCDF4 import Dataset

class DataArray(np.ndarray):
    def __new__(cls, input_array, attrs={}):
        obj = np.asarray(input_array).view(cls)
        obj.attrs = attrs
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.attrs = getattr(obj, 'attrs', {
            'units': None,
            'long_name': None
        })

        
def read_netcdf_one_file(filename, file_vars = None, verbose = False):
    """Read all data from an Aether netcdf file.

    Parameters
    ----------
    filename : str
        An Aether netCDF filename
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
        raise IOError(f"unknown aether netCDF blocked file: {filename}")

    # NOTE: Includes header information for easy access until
    #       updated package structure is confirmed
    # Initialize data dict with defaults (will remove these defaults later)
    data = {'filename': filename,
            'units': {},
            'long_name': None}
    
    if verbose:
        print('-> Reading netcdf : ', filename, ' --> Vars : ', file_vars)

    with Dataset(filename, 'r') as ncfile:
        # Process header information: nlons, nlats, nalts, nblocks
        data['nlons'] = len(ncfile.dimensions['lon'])
        data['nlats'] = len(ncfile.dimensions['lat'])
        data['nalts'] = len(ncfile.dimensions['z'])
        try:
            data['nblocks'] = len(ncfile.dimensions['block'])
        except:
            data['nblocks'] = 0
            
        # Included for compatibility
        data['vars'] = [var for var in ncfile.variables.keys()
                        if file_vars is None or var in file_vars]

        # Fetch requested variable data
        for key in data['vars']:
            var = ncfile.variables[key]  # key is var name
            data[key] = DataArray(np.array(var), var.__dict__)
            data['units'][key] = var.units if 'units' in var.__dict__ else ''

        if 'since' in ncfile.variables['time'].units:
            t0 = ncfile.variables['time'].units.split('since')[-1].strip()
            t0 = datetime.strptime(t0, '%Y-%m-%d')
            if verbose:
                print('   -> Time conversion using t0 = ', t0)
        else:
            t0 = datetime(1965, 1, 1)

        data['times'] = \
            tc.epoch_to_datetime(np.array(ncfile.variables['time'])[0], t0=t0)

        try:
            data['isEnsemble'] = True if ncfile.isEnsemble == "True" else False
        except:
            data['isEnsemble'] = False
        if data['isEnsemble']:
            data['ensembleNumber'] = int(ncfile.ensembleNumber)
            data['ensembleMembers'] = int(ncfile.ensembleMembers)
            
    return data


def read_netcdf_one_header(filename):
    """Read all keys and such from netcdf file

    Parameters
    ----------
    filename : str
        An Aether netCDF filename
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
        data['nalts'] = len(ncfile.dimensions['z'])
        try:
            data['nblocks'] = len(ncfile.dimensions['block'])
        except:
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
        
        if 'since' in ncfile.variables['time'].units:
            t0 = ncfile.variables['time'].units.split('since')[-1].strip()
            t0 = datetime.strptime(t0, '%Y-%m-%d')
            
        else:
            t0 = datetime(1965, 1, 1)
        data['times'] = \
            tc.epoch_to_datetime(np.array(ncfile.variables['time']), t0=t0)

        try:
            data['isEnsemble'] = True if ncfile.isEnsemble == "True" else False
        except:
            data['isEnsemble'] = False
        if data['isEnsemble']:
            data['ensembleNumber'] = int(ncfile.ensembleNumber)
            data['ensembleMembers'] = int(ncfile.ensembleMembers)

    data['shortname'] = variables.get_short_names(data['vars'])
    
    return data


# -----------------------------------------------------------------------------
# This reads in a series of vars / files and returns the 3D information
#-----------------------------------------------------------------------------

def read_netcdf_all_files(filelist, varlist=[-1], verbose=False):

    filelist = util.any_to_filelist(filelist)

    # Get the prefixes for all entries in filelist;
    prefixes = np.unique([file.split('/')[-1].split('_')[0] \
                          for file in filelist])
    if len(prefixes) > 1: # make sure there is only one output type.
        raise ValueError("Multiple output types cannot be read by this " +
                         "function.\n\tProvided: " + str(prefixes))

    # first read in spatial information:
    vars = ['lon', 'Longitude', 'lat', 'Latitude', 'z', 'Altitude']
    spatialData = read_netcdf_one_file(filelist[0], vars, verbose=False)

    nTimes = len(filelist)
    if varlist != [-1]:
       nVars = len(varlist)
    else: # varlist=[-1] means we read in all variables
        header = read_netcdf_one_header(filelist[0], verbose=False)
        varlist = header['vars']
        nVars = len(varlist)

    allTimes = []
    if (spatialData['nblocks'] == 0):
        # This assumes we have 3D arrays for the coord info.
        # GITM will put 1D arrays into lon/lat/z if it can, which we don't want.   
        lons = spatialData['Longitude' if 'Longitude' in spatialData.keys() else 'lon']
        lats = spatialData['Latitude' if 'Latitude' in spatialData.keys() else 'lat']
        alts = spatialData['Altitude' if 'Altitude' in spatialData.keys() else 'z'] / 1000.0  # Convert from m to km
        nDims = len(np.shape(lons))
        if (nDims == 4):
            nLons = len(lons[0, :, 0, 0])
            nLats = len(lats[0, 0, :, 0])
            nAlts = len(alts[0, 0, 0, :])
        else:
            nLons = len(lons[:, 0, 0])
            nLats = len(lats[0, :, 0])
            nAlts = len(alts[0, 0, :])
        nBlocks = 0
        
        if (nVars == 1):
            allData = np.zeros((nTimes, nLons, nLats, nAlts))
        else:
            allData = np.zeros((nTimes, nVars, nLons, nLats, nAlts))

    else:
            
        # we will now have a block dimension, and the latitude and
        # longitude could be dependent on block, lon, and lat:
        lons = spatialData['lon']
        nLons = len(lons[0, :, 0, 0])
        lats = spatialData['lat']
        nLats = len(lats[0, 0, :, 0])
        alts = spatialData['z'] / 1000.0  # Convert from m to km
        nAlts = len(alts[0, 0, 0, :])
        nBlocks = len(lons[:, 0, 0, 0])
        
        if (nVars == 1):
            allData = np.zeros((nTimes, nBlocks, nLons, nLats, nAlts))
        else:
            allData = np.zeros((nTimes, nVars, nBlocks, nLons, nLats, nAlts))

    for iTime, filename in enumerate(filelist):
        data = read_netcdf_one_file(filename, varlist, verbose=verbose)
        allTimes.append(data["times"])
        for iVar, var in enumerate(varlist):
            if (nBlocks == 0):
                if (nVars == 1):
                    allData[iTime, :, :, :] = data[var][:, :, :]
                else:
                    allData[iTime, iVar, :, :, :] = data[var][:, :, :]
            else:
                if (nVars == 1):
                    allData[iTime, :, :, :, :] = data[var][:, :, :, :]
                else:
                    allData[iTime, iVar, :, :, :, :] = data[var][:, :, :, :]
                
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


