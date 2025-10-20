#!/usr/bin/env python3

import datetime
import re, os
from glob import glob

from pyitm.fileio import gitmio, netcdfio, variables, satelliteio, madrigalio
from pyitm.modeldata import utils
import numpy as np
from glob import glob

# ----------------------------------------------------------------------------
# Helpful functions for the io modules
# ----------------------------------------------------------------------------

def notfound(err):
    print(err)
    raise FileNotFoundError

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def determine_filetype(filename):
    """
    Determine the filetype for a given file. 
    There are currently only two supported types:
        - gitm binary files
        - aether netcdf files
    """
    fType = {
        "iGitmBin": 0,
        "iGitmNetcdf": 1,
        "iNetcdf": 2,
        "myfile": -1
    }
    m = re.match(r'(.*)bin', filename)
    if m:
        fType["myfile"] = fType["iGitmBin"]
    else:
        fType["myfile"] = fType["iNetcdf"]
    return fType

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def find_files_in_different_directory(filelist, dir):

    filelist = any_to_filelist(filelist)

    filelistOut = []
    for file in filelist:
        ft = determine_filetype(file)
        if (ft['myfile'] == ft["iGitmBin"]):
            ext = '.bin'
        else:
            ext = '.nc'
        m = re.match(r'(.*)(\d)'+ext, file)
        if m:
            isGood = False
            newfile = glob(dir + '/' + m.group(1) + '?' + ext)
            if (np.isscalar(newfile)):
                if (len(newfile) > 0):
                    filelistOut.append(newfile)
                    isGood = True
            else:
                if (len(newfile[0]) > 0):
                    filelistOut.append(newfile[0])
                    isGood = True
            if (not isGood):
                print('Can not file background file : ',file)

    return filelistOut

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def determine_if_tec(varsToRead):
    isTec = False
    if varsToRead == []:
        return isTec

    if (isinstance(varsToRead, str)):
        if (varsToRead.lower() == 'tec'):
            isTec = True
    else:
        if (isinstance(varsToRead[0], str)):
            if (varsToRead[0].lower() == 'tec'):
                isTec = True
    return isTec

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def determine_if_on2(varsToRead):
    isOn2 = False
    if varsToRead == []:
        return isOn2

    if (isinstance(varsToRead, str)):
        if (varsToRead.lower() == 'on2'):
            isOn2 = True
        if (varsToRead.lower() == 'o/n2'):
            isOn2 = True
    else:
        if (isinstance(varsToRead[0], str)):
            if (varsToRead[0].lower() == 'on2'):
                isOn2 = True
            if (varsToRead[0].lower() == 'o/n2'):
                isOn2 = True
    return isOn2

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def read_all_files(filelist, varsToRead = None, verbose = False):
    
    filelist = any_to_filelist(filelist)
    filetype = determine_filetype(filelist[0])
    header = read_all_headers(filelist[0])
    if (varsToRead == None):
        varsToRead = header['vars']
        isTec = False
        isOn2 = False
    else:
        isTec = determine_if_tec(varsToRead)
        isOn2 = determine_if_on2(varsToRead)

    if (isTec):
        varsToRead = ['e-']
    if (isOn2):
        varsToRead = ['o', 'n2']

    if (filetype["myfile"] == filetype["iGitmBin"]):
        varsToRead = variables.convert_number_to_var(varsToRead, header)
        varsToRead = variables.get_short_names(varsToRead)
        varsToRead = variables.convert_var_to_number(varsToRead, header)
        allData = gitmio.read_gitm_all_files(filelist, varsToRead, verbose=verbose)
    if (filetype["myfile"] == filetype["iNetcdf"]):
        varsToRead = variables.convert_number_to_var(varsToRead, header)
        varsToRead = variables.match_var_name(varsToRead, header)
        if ('NotFound' in varsToRead):
            allData = None
        else:
            allData = netcdfio.read_netcdf_all_files(filelist, varsToRead, verbose=verbose)
    if (isTec):
        tec = utils.calc_tec(allData)
        allData['tec'] = tec
        
    return allData

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def read_all_headers(filelist, verbose = False):

    filelist = any_to_filelist(filelist)
    filetype = determine_filetype(filelist[0])
    isFound = False
    if (filetype["myfile"] == filetype["iGitmBin"]):
        if (verbose):
            print(' -> Reading GITM header(s)')
        header = gitmio.read_gitm_headers(filelist)
        isFound = True
    if (filetype["myfile"] == filetype["iNetcdf"]):
        if (verbose):
            print(' -> Reading NETCDF header --- Can only read one at a time!')
        header = netcdfio.read_netcdf_one_header(filelist[0])        
        isFound = True
        
    if (not isFound):
        print('File type not supported at this point!')
        header = None

    return header


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def list_file_info(filelist):
    
    filelist = any_to_filelist(filelist)
    header = read_all_headers(filelist[0])
    print('Variables in file (num, var, short, long)')
    for iVar, var in enumerate(header['vars']):
        print('%3d' % iVar, '. ', var, ' -> ', \
              header['shortname'][iVar], ' -> ',
              header['longname'][iVar])
    print('Time of file : ', header['times'])
    return

    
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def any_to_filelist(input_data=None):
    """ 
    Take input_data and convert it to a list of file paths
    
    Parameters
    ----------
    input_data (str, list-like, path): raw user input. 
        Can be str, list, str with *'s, etc. This is what we will work on
    
    Returns
    -------
    list of strings: Hopefully the files that the user was trying to input...

    Notes
    -----
    - Will error on unexpected data
    - Will error if no files can be found
    - Can handle:
        - list to one or multiple files
        - str with '*'
        - str to single file
        - str without '*', like "path/to/3DALL"
        - str to directory, we take all '*.bin' files there

    """
    
    from glob import glob
    import os

    if isinstance(input_data, (str, os.PathLike)):
        # Single file that exists:
        if os.path.isfile(input_data) and os.path.exists(input_data):
            return [input_data]
        # directory that exists, return all .bin files inside
        elif os.path.isdir(input_data):
            return sorted(glob(os.path.join(input_data, "*.bin")))
        # User gave us a string to glob:
        elif '*' in input_data:
            outfiles = sorted(glob(input_data))
            # make sure the files actually exist
            if len(outfiles) > 0:
                return outfiles
            else:
                notfound("glob pattern in input_data, but no files found\n"
                         f"Received:\n\t{input_data}")
        # maybe the user wants us to glob for them?
        else:
            # directory is already handled, it must be a path or something
            for fileext in ['*.bin', '*.nc']:
                outfiles = sorted(glob(input_data + fileext))
                if len(outfiles) > 0:
                    return outfiles
            notfound(f"No {fileext} files found within:\n\t{input_data}")

    # we were probably given a list of files already
    else:
        # If it is only one string, send that back in here as string, not list.
        if len(input_data) == 1:
            # This will loop infinitely if we don't first check if that one file exists
            if os.path.exists(input_data[0]):
                return input_data
            else:
                return any_to_filelist(input_data[0])
        else: #list of files is longer that one
            return input_data

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def read_satfiles(filelist=None, satname=None, 
                  satLookup=None, startDate=None, endDate=None,
                  verbose=False):
    """Read multiple satellite files.
    
    Parameters
    ----------
    filelist : list or str  
        List of file paths or a single file path to read.

    satname : str, optional
        Name of the satellite to read data from. If None, will read all satellites.

    satLookup : str, optional
        Path to a CSV file containing a lookup table for satellite names, directories,
        and filenames. This can help in identifying and organizing data from multiple
        satellites.
    
    startDate : datetime, or date-like, optional
        Start date to filter satellite files.
        
    endDate : datetime, or date-like, optional
        End date to filter satellite files. 

    verbose : bool, optional
        If True, print additional information while reading files.
    
    Returns
    -------
    dict
        Dictionary containing combined data from all satellite files. 
        - If only one satellite is read, the dictionary will contain data for that
          satellite. If multiple satellites are read, the dictionary will have keys
          for each satellite name, each containing their respective data.

    """

    if satLookup is not None:
        filelist, satnames = lookup_satfiles(satLookup, date_start=startDate,
                                             date_end=endDate, satname=satname,
                                             verbose=verbose)
    else:
        filelist = any_to_filelist(filelist)
    
    if verbose:
        print(f"Reading {len(filelist)} satellite files...")

    all_data = []
    satnames = []
    read_files = []
    for f in filelist:
        # Do not re-read, some files are monthly.
        if f in read_files:
            continue
        if verbose:
            print(f"-> Reading file: {f}")
        data = satelliteio._read_sat_one_file(f, satname=satname, verbose=verbose)
        satnames.append(data.pop('sat_name'))
        all_data.append(data)
        read_files.append(f)
    
        # Here we separate multiple satellites from a constellation.
        # This has to be manual since naming conventions vary.
        if 'dmsp' in satnames[-1].lower():
            # dms_20150716_18s1.001.* (or) dms_ut_20240515_18.002.hdf5
            # Keep precip/density! cols are different so the reader needs to keep them separate
            satnames[-1] = 'DMSP_F'+f.split('/')[-1].split('_')[-1][:2]+'-'+satnames[-1].split('_')[-1]
        if 'grace' in satnames[-1].lower():
            # ga_dns_...
            satnames[-1] = 'GRACE'+f.split('/')[-1][1]

    combined_data = {}
    unique_names = []
    for name, data in zip(satnames, all_data):
        # Make sure satellite is in combined_data
        if name not in combined_data:
            combined_data[name] = {}
            unique_names.append(name)
        
        # add the satellite data
        for key, value in data.items():
            if key in combined_data[name].keys():
                combined_data[name][key] = np.concatenate((combined_data[name][key], value))
            else:
                combined_data[name][key] = value
    
    if verbose:
        print(f" -> Found sat data from satellites: {combined_data.keys()}")
        for satname in combined_data.keys():
            print(f" --> data['{satname}'] has keys: {combined_data[satname].keys()}")

    return combined_data

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def lookup_satfiles(lookup_file, date_start, date_end=None, satname=None, verbose=False):
    """
    Find satellite files in a directory.

    Parameters
    ----------
    lookup_file : str
        Path to CSV file containing satellite file information. The CSV should have columns
        for satellite name, directory, and filename patterns. Date should be provided to 
        avoid reading in all files.
    date_start : datetime, list-like, or date-like
        Start date to search for files. Or a list of dates!
    date_end : datetime, or date-like, optional
        End date to search for files.
    satname : str, optional
        If given, only return files that contain this string.
    verbose : bool, optional
        If True, print information while searching.

    Returns
    -------
    list of str
        List of paths to satellite files found in the directory.
    """

    # deal with dates
    dates = []
    if (verbose):
        print('Inside lookup_satellite: ')
        print(' -> figuring out dates: ', date_start)
    if isinstance(date_start, (datetime.datetime, datetime.date)):
        dates = [date_start]
        if date_end is not None:
            if (verbose):
                print('   -> end date: ', date_end)
            
            if date_end < date_start:
                raise ValueError("date_end must be after date_start")
            if date_end > date_start:
                tday = date_start + datetime.timedelta(days=1)
                while tday <= date_end:
                    dates.append(tday)
                    tday += datetime.timedelta(days=1)
    else:
        # list-like input?
        try:
            date_start = list(date_start)
        except:
            raise TypeError("date_start must be datetime, or date-like, or list-like")
        

    if (verbose):
        print('  --> dates to search for : ', dates)
        
    # Open the lookup table:
    with open(lookup_file, 'r') as f:
        lines = f.readlines()
    iline = 0
    if '/' in lines[0] or '.' in lines[0]:
        # it's a header.
        iline = 1
    
    satnames = []
    sat_paths = []
    sat_filenames = []

    for iline in range(iline, len(lines)):
        line = lines[iline].strip()
        if line.startswith('#') or line == '':
            continue
        parts = line.split(',')
        if len(parts) < 2:
            continue
        satnames.append(parts[0].strip())
        sat_paths.append(parts[1].strip())
        if len(parts) > 2:
            sat_filenames.append(parts[2].strip())
        else:
            sat_filenames.append('*')
    if verbose:
        print(f"-> Found {len(satnames)} entries in lookup table.")

    if satname is not None:
        if verbose:
            print('  --> Looking at satnames : ')
            for sn in satnames:
                print('    ', sn)
        if satname in satnames:
            isat = satnames.index(satname)
            satnames = satnames[isat]
            sat_paths = sat_paths[isat]
            sat_filenames = sat_filenames[isat]
            if verbose:
                print('  --> Satellite found: ', satname)
                print('      isat : ', isat)
                print('      sat_paths : ', sat_paths)
                print('      sat_filenames : ', sat_filenames)
        else:
            raise KeyError(
                f"From lookup_satfiles '{satname}' not found!"
                f"  Expected one of: {satnames}")

    satfiles_out = []
    satnames_out = []
    sat = satnames
    sdir = sat_paths
    sfname = sat_filenames
    #for sat, sdir, sfname in zip(satnames, sat_paths, sat_filenames):
    if verbose:
        print(f"  --> Searching for files in {sdir}")
    for eachdate in dates:
        p = eachdate.strftime(os.path.join(sdir, sfname))
        if verbose:
            print('  --> search path : ', p)
        days_files = sorted(glob(p))
        if verbose:
            print(f"   --> Found {len(days_files)} files for {eachdate}")
        if len(days_files) == 0:
            continue
        # append to satfiles, making sure it is 1D
        satfiles_out = satfiles_out + days_files
        satnames_out = satnames_out + [sat]*len(days_files)

    return satfiles_out, satnames_out
