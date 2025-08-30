#!/usr/bin/env python3

import re

from pyitm.fileio import gitmio, netcdfio, variables
import numpy as np

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

def read_all_files(filelist, varsToRead = None, verbose = False):

    filelist = any_to_filelist(filelist)
    filetype = determine_filetype(filelist[0])
    header = read_all_headers(filelist[0])
    if (varsToRead == None):
        varsToRead = header['vars']
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
    print('Time of file : ', header['time'])
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
        if os.path.isfile(input_data):
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
            outfiles = sorted(glob(input_data + "*.bin"))
            if len(outfiles) > 0:
                return outfiles
            else:
                notfound(f"No *.bin files found within:\n\t{input_data}")

    # we were probably given a list of files already
    else:
        return input_data
