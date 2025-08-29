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

def read_all_files(filelist, varToPlot):

    filelist = any_to_filelist(filelist)
    filetype = determine_filetype(filelist[0])
    print(filelist[0])
    header = read_header(filelist[0])
    if (filetype["myfile"] == filetype["iGitmBin"]):
        varToPlot = variables.get_short_names(varToPlot)
        varToPlot = variables.convert_var_to_number(varToPlot, header)
        print('vartoplot : ', varToPlot)
        allData = gitmio.read_gitm_all_files(filelist, varToPlot)
    if (filetype["myfile"] == filetype["iNetcdf"]):
        varToPlot = variables.convert_number_to_var(varToPlot, header)
        varToPlot = variables.match_var_name(varToPlot, header)
        if ('NotFound' in varToPlot):
            allData = None
        else:
            allData = netcdfio.read_netcdf_all_files(filelist, varToPlot)
    return allData

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def read_header(filelist):

    filelist = any_to_filelist(filelist)
    filetype = determine_filetype(filelist[0])
    isFound = False
    if (filetype["myfile"] == filetype["iGitmBin"]):
        header = gitmio.read_gitm_headers(filelist[0])
        isFound = True
    if (filetype["myfile"] == filetype["iNetcdf"]):
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
    header = read_header(filelist)
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
