"""
Utilities for reading various files from Madrigal.

NetCDF4 and HDF5 files must be read differently.

"""

from pyitm.general.time_conversion import epoch_to_datetime
from pyitm.modeldata.satellite import calc_wind_dir
import datetime
import numpy as np

# Standardize variable names if possible
var_map = {'gdlat': 'lats',
           'glon': 'lons',
           'gdalt': 'alts',
           'mlong': 'mlon',
           'dtec': 'tec_error',
           'tec': 'tec',
           'hor_ion_v': 'Vi_hor',
           'vert_ion_v': 'Viv',
           'el_m_ener': 'AveE',
           'ion_m_ener': 'AveE_I',
           'el_i_flux': 'eFlux',
           'ion_i_flux': 'eFlux_I',
           'ne': 'e-',
           'ni': 'e-',
           'te': 'Te',
           'ti': 'Ti',
}

def _read_madrigal_one_file(filename, verbose=False):
    """
    Read a single Madrigal file, either netCDF or HDF5. This only reads a single file.
    Use pyitm.fileio.util.read_satfiles to read multiple files.

    Parameters
    ----------
    filename (str): path to the Madrigal file

    verbose (bool): print debugging info? default=False

    Returns
    -------
    data (dict): dictionary containing the data in the file.
    """

    # This can be changed to the same dispatcher in satelliteio,
    # but i don't think we need more than these two readers.
    if filename.endswith('.nc'):
        return _read_madrigal_nc_file(filename, verbose=verbose)
    elif filename.endswith('.hdf5') or filename.endswith('.h5'):
        return _read_madrigal_hdf5_file(filename, verbose=verbose)
    else:
        raise ValueError("File must be a .nc or .hdf5/.h5 file.")

    
def _read_madrigal_nc_file(filepath, verbose=False):
    """
    Read a Madrigal netCDF file. This is called automatically from
      '_read_madrigal_one_file' if the filetype is .nc

    Parameters
    ----------
    filepath (str): path to the netCDF file

    verbose (bool): print debugging info? default=False

    Returns
    -------
    data (dict): dictionary containing the data in the file.
    """

    # in case someone doesn't have netCDF4 installed but wants to use hdf5 reader
    import netCDF4

    if verbose:
        print(f"Reading Madrigal netCDF file: {filepath}")

    with netCDF4.Dataset(filepath, 'r') as ds:
        data = {}
        for varname in ds.variables.keys():
            data[varname] = ds.variables[varname][:]
            if verbose:
                print(f"-> Read variable '{varname}' with shape {data[varname].shape}.")

    # Rename variables, reformatting as needed
    for oldname in data.keys():
        if oldname in var_map:
            data[var_map[oldname]] = data.pop(oldname)
            if verbose:
                print(f"--> Renamed variable '{oldname}' to '{var_map[oldname]}'.")
        
        # Convert time variable to datetime objects if possible
        # Key could be time/times/timestamps
        if 'time' in oldname:
            data['times'] = epoch_to_datetime(data.pop(oldname), datetime.datetime(1970, 1, 1))
            if verbose:
                print(f"--> Found time! Start: {data['times'][0]}, End: {data['times'][-1]}.")

    return data

def extrainfo(v):
    """
    This is for the verbose mode of reading madrigal hdf5 files.
    The whole reader is kinda bonkers. Can't iterate recursively without a function.
    """
    print(' -> ', v)
    return

def _read_madrigal_hdf5_file(filepath, verbose=False):
    """
    Read data from a Madrigal HDF5 file. This is called automatically from
      '_read_madrigal_one_file' if the filetype is .hdf5 or .h5

    Parameters
    ----------
    filepath (str): path to the HDF5 file
    verbose (bool): print debugging info? default=False

    Returns
    -------
    data (dict): dictionary containing the data in the file.

    Notes
    -----
    This is a somewhat rough reader, as the structure of Madrigal HDF5 files can vary
    between data sources. More work may be necessary to make this more robust, but it
    works for TEC files. Run with verbose to see extra info, and reach out with
    questions!


    """

    # in case someone doesn't have h5py installed but wants to use netcdf4 reader
    import h5py

    if verbose:
        print(f"Reading Madrigal HDF5 file: {filepath}")
    
    with h5py.File(filepath, 'r') as f:
        if verbose:
            print(" File structure:")
            f.visit(extrainfo)

        #This is for multi-dimensional data, like TEC or precipitation
        if 'Data/Array Layout' in f:
            if verbose:
                print("Found 'Data/Array Layout' in file, assuming multi-dimensional data.")
            datavars = {} # {datavars:lookup_name}
            for param_dim in range(1, 4):
                lookup = f"Data/Array Layout/{param_dim}D Parameters/Data Parameters"
                if lookup in f:
                    if verbose:
                        print(f"Found {param_dim}D parameters:")

                    for p in f[lookup]:
                        varname = p['mnemonic'].decode()
                        datavars[varname] = '/'.join(lookup.split('/')[:-1])+"/"+varname.lower()
                        if verbose:
                            print(" - ", varname, p['description'].decode())
                else:
                    if verbose:
                        print(f"No {param_dim}D parameters found.")

            data = {}
            for varname, lookup in datavars.items():
                try:
                    newname = var_map[varname.lower()] if varname.lower() in var_map else varname.lower()
                    data[newname] = f[lookup][:]
                    if verbose:
                        print(f"-> Read variable '{varname}'-> {newname} with shape {data[newname].shape}.")
                except Exception as e:
                    if verbose:
                        print(f"-> Could not read variable '{varname}' from '{lookup}': {e}")

            # OK FUN!
            # As asinine as the rest of this was, it's worse when we have to read time.
            # There's no metadata so we can't pull it out with the same methodology.
            # AND the times are stored as 64-bit integers, which can't be used as the
            # seconds component in datetime.timedeltas! :)
            epochtimes = [int(i) for i in f['Data/Array Layout']['timestamps']]
            data['times'] = epoch_to_datetime(epochtimes, datetime.datetime(1970, 1, 1))
        
        # 1D files are just a table, much easier!
        elif 'Data/Table Layout' in f:
            data = {}
            times = {}
            for varname_b in f['Metadata/Data Parameters']:
                varname = varname_b['mnemonic'].decode().lower()
                newname = var_map[varname.lower()] if varname.lower() in var_map else varname.lower()
                # Store time separately for now...
                if varname in ['year', 'month', 'day', 'hour', 'minute', 'min', 'second', 'sec']:
                    try:
                        times[varname] = f['Data/Table Layout'][varname]
                        if verbose:
                            print(f"-> Read time component '{varname}' with shape {times[varname].shape}.")
                    except KeyError:
                        continue
                else:
                    data[newname] = f['Data/Table Layout'][varname]
                    if verbose:
                        print(f"-> Read variable '{varname}'->{newname} with shape {data[newname].shape}.")
            
            # Convert ymdhms columns to times key in data
            for k in ['year', 'month', 'day', 'hour', 'minute', 'second']:
                if k not in times:
                    # If we do not find something (minute), look for first 3 chars (min)
                    first3 = k[:3]
                    if first3 in times:
                        times[k] = times[first3]
                    else:
                        times[k] = np.zeros_like(times[list(times.keys())[0]])
            data['times'] = np.array([datetime.datetime(int(y), int(m), int(d), int(h), int(mi), int(s))
                             for y, m, d, h, mi, s in zip(times['year'], times['month'], times['day'],
                                                  times['hour'], times['minute'], times['second'])])
            if verbose:
                print(f"-> Constructed {len(data['times'])} datetime objects for 'times' key.")

        # Something else entirely! This is basically just a dict of data & is easy-ish
        else:
            if verbose:
                print(f"-> Looks to be a flat file. Iterating through the columns.")
            data = {}
            for varname in f.keys():
                newname = var_map[varname.lower()] if varname.lower() in var_map else varname.lower()
                data[newname] = f[varname][:]

            # Still need to clean up timestamps...
            data['times'] = np.array([datetime.datetime(1970,1,1) 
                                      + datetime.timedelta(seconds=dt) for dt in data.pop('timestamps')])
            if verbose: 
                print(f" --> Found columns: {data.keys()}")

    # More cleaning!! If we found velocity, it's probably horizontal velocity
    # Let's add a direction to it.
    if 'Vi_hor' in data.keys():
        uniteRot, unitnRot = calc_wind_dir(np.array(data["lons"]),
                                           np.array(data["lats"]))
        
        data["Vie"] = uniteRot * data['Vi_hor']
        data["Vin"] = unitnRot * data['Vi_hor']

    return data