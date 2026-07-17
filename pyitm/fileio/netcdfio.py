#!/usr/bin/env python3

import os
from datetime import datetime
from struct import unpack
import numpy as np
from pyitm.fileio import util
from pyitm.fileio import variables
from pyitm.general import time_conversion as tc

from netCDF4 import Dataset

def check_whether_ipe(filename):

    state = False

    with Dataset(filename, 'r') as ncfile:
        # Process header information: nlons, nlats, nalts, nblocks
        if (('x01' in ncfile.dimensions) and \
            ('x02' in ncfile.dimensions) and \
            ('x03' in ncfile.dimensions)):
            state = True
    ncfile.close

    return state

def check_whether_aether(filename):

    state = False

    with Dataset(filename, 'r') as ncfile:
        # Process header information: nlons, nlats, nalts, nblocks
        if (('lon' in ncfile.dimensions) and \
            ('lat' in ncfile.dimensions) and \
            ('z' in ncfile.dimensions)):
            state = True
    ncfile.close
    return state

class DataArray(np.ndarray):
    def __new__(cls, input_array, attrs={}):
        obj = np.asarray(input_array).view(cls)
        obj.attrs = attrs
        return obj

    def __array__(self, dtype=None, copy=False):
        return np.array(self.view(np.ndarray), dtype=dtype, copy=copy)

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.attrs = getattr(obj, 'attrs', {
            'units': None,
            'long_name': None
        })


def read_nc_times(ncfile, verbose=False):
    """Convert a file's time variable to datetimes, handling
    "[seconds/minutes/hours/days] since [date [time]]" units.
    """
    time_var = ncfile.variables['time']
    units = getattr(time_var, 'units', '')
    t0 = datetime(1965, 1, 1)
    scale = 1.0
    if 'since' in units:
        head, t0_str = [s.strip() for s in units.split('since', 1)]
        scale = {'days': 86400.0, 'hours': 3600.0,
                 'minutes': 60.0, 'seconds': 1.0}.get(head.lower(), 1.0)
        for fmt in ('%Y-%m-%d %H:%M:%S', '%Y-%m-%d'):
            try:
                t0 = datetime.strptime(t0_str, fmt)
                break
            except ValueError:
                continue
        else:
            raise ValueError(f"Unrecognized time units: {units}")
        if verbose:
            print('   -> Time conversion using t0 = ', t0)
    return tc.epoch_to_datetime(time_var[:] * scale, t0=t0)


def read_netcdf_inventory(filelist, verbose = False):
    """Summarize netcdf files without reading any data.

    Returns a dict with:
        filenames - the (normalized) filelist
        times - global sorted list of datetimes across all files
        ifile, ilocal - for each global time, which file and which time
            index within that file it comes from
        ntimes - len(times)
        varinfo - per-variable dict: dims, shape, dtype, units, longname,
            hastime, nbytes (bytes on disk for one time)
        groups - data variables grouped by their dims minus time;
            variables named after a dimension are coordinates and excluded
    """

    filelist = util.any_to_filelist(filelist)

    varinfo = {}
    groups = {}
    entries = []
    for iFile, filename in enumerate(filelist):
        with Dataset(filename, 'r') as ncfile:
            if 'time' not in ncfile.variables:
                raise ValueError(f"no time variable in {filename}")
            for iLocal, t in enumerate(read_nc_times(ncfile)):
                entries.append((t, iFile, iLocal))
            if iFile > 0:
                continue
            for name, var in ncfile.variables.items():
                dims = var.dimensions
                spatial = tuple(d for d in dims if d != 'time')
                nPoints = int(np.prod([len(ncfile.dimensions[d])
                                       for d in spatial], dtype = np.int64))
                varinfo[name] = {'dims': dims,
                                 'shape': var.shape,
                                 'dtype': str(var.dtype),
                                 'units': getattr(var, 'units', ''),
                                 'longname': getattr(var, 'long_name', name),
                                 'hastime': 'time' in dims,
                                 'nbytes': nPoints * var.dtype.itemsize}
                if name not in ncfile.dimensions:
                    groups.setdefault(spatial, []).append(name)

    entries.sort(key = lambda e: (e[0], e[1], e[2]))
    inventory = {'filenames': filelist,
                 'times': [e[0] for e in entries],
                 'ifile': np.array([e[1] for e in entries]),
                 'ilocal': np.array([e[2] for e in entries]),
                 'ntimes': len(entries),
                 'varinfo': varinfo,
                 'groups': groups}

    if verbose:
        print(f' -> Inventory: {len(filelist)} files, {len(entries)} times, '
              f'{len(varinfo)} variables')
        for sig, names in groups.items():
            print('   -> vars on', sig if sig else '(scalar)', ':', names)

    return inventory


def estimate_read_size(inventory, varlist = None, nTimes = None):
    """Estimated bytes on disk to read varlist over nTimes times."""
    if varlist is None or varlist == [-1]:
        varlist = [v for names in inventory['groups'].values() for v in names]
    if nTimes is None:
        nTimes = inventory['ntimes']
    return sum(inventory['varinfo'][v]['nbytes'] for v in varlist) * nTimes


def read_netcdf_one_file(filename, file_vars = None, verbose = False):
    """Read all data from an Aether netcdf file.

    Parameters
    ----------
    filename : str
        A netCDF filename
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

        for dim in ['lon', 'lat', 'z']:
            if dim not in data['vars']:
                data['vars'].append(dim)

        # Fetch requested variable data
        for key in data['vars']:
            if verbose:
                print('   -> Reading variable : ', key)
            try:
                var = ncfile.variables[key]  # key is var name
                data[key] = DataArray(var[:], var.__dict__)
                data['units'][key] = var.units if 'units' in var.__dict__ else ''
            except KeyError:
                if key =='z':
                    data[key] = DataArray(np.array([100]), )
                    data['units'][key] = 'km'
                else:
                    raise 

        data['times'] = read_nc_times(ncfile, verbose=verbose)

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
        data['nvars'] = len(data['vars'])

        data['times'] = read_nc_times(ncfile)

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
    # vars = ['lon', 'Longitude', 'lat', 'Latitude', 'z', 'Altitude']
    # spatialData = read_netcdf_one_file(filelist[0], vars, verbose=False)

    header = read_netcdf_one_header(filelist[0])
    if len(filelist)==1:
        nTimes = len(header['times'])
    else:
        # files may hold one or more times each
        nTimes = sum(len(read_netcdf_one_header(f)['times'])
                     for f in filelist)
    if varlist != [-1]:
       nVars = len(varlist)
    else: # varlist=[-1] means we read in all variables
        varlist = header['vars']
        nVars = len(varlist)

    nBlocks = header['nblocks']
    nLons = header['nlons']
    nLats = header['nlats']
    nAlts = header['nalts']

    allTimes = []

    # Make output holder! its shape is conditional. Order of axis:
    # nTimes, nVars, nBlocks, nLons, nLats, nAlts
    # If nBlocks==1, it's squeezed. Same for nVars (like gitmio)
    out_shape = []
    out_shape.append(nTimes)
    if nVars > 1:
        out_shape.append(nVars)
    if nBlocks > 1:
        out_shape.append(nBlocks)
    out_shape.append(nLons)
    out_shape.append(nLats)
    out_shape.append(nAlts)

    allData = np.zeros(out_shape)

    # We may be reading a file with multiple times...
    # If multiiple times are in one file, we can advance time independent from 
    # the filelist loop
    iAllTimes = 0
    nSpatialDims = len(out_shape) - 2 if nVars > 1 else len(out_shape) - 1
    for filename in filelist:
        data = read_netcdf_one_file(filename, varlist, verbose=verbose)
        for iTime in range(len(data["times"])):
            allTimes.append(data["times"][iTime])
            for iVar, var in enumerate(varlist):
                val = data[var]
                # some files have no time axis on the variables
                if val.ndim > nSpatialDims:
                    val = val[iTime, ...]
                if (nVars == 1):
                    allData[iAllTimes, ...] = val
                else:
                    allData[iAllTimes, iVar, ...] = val
            iAllTimes += 1
    vars = []
    lons = data.pop('Longitude' if 'Longitude' in data.keys() else 'lon')
    lats = data.pop('Latitude' if 'Latitude' in data.keys() else 'lat')
    alts = data.pop('Altitude' if 'Altitude' in data.keys() else 'z')
    # Coordinates may already be multi-dimensional (blocked grids store per-point values).
    # Only meshgrid when all three are 1D coordinate vectors.
    if lons.ndim == 1 and lats.ndim == 1 and alts.ndim == 1:
        lons, lats, alts = np.meshgrid(lons, lats, alts, indexing='ij')
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


