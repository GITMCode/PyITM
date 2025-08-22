#!/usr/bin/env python3

import re

from pyitm.fileio import gitmio
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
        "iAetherNetcdf": 2,
        "myfile": -1
    } 
    m = re.match(r'(.*)bin', filename)
    if m:
        fType["myfile"] = fType["iGitmBin"]
    else:
        fType["myfile"] = fType["iAetherNetcdf"]
    return fType

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def find_string(item, stringList):
    iVal = -1
    if (item in stringList):
        i = 0
        while (i < len(stringList)):
            if (stringList[i] == item):
                iVal = i
                i = len(stringList)
            i += 1
    return iVal

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def convert_var_to_number(varList, header = None):

    if (np.isscalar(varList)):
        if (varList.isnumeric()):
            iVars = [int(varList)]
        else:
            if (header == None):
                print('Non number variables are not supported yet!')
                iVars = [3]
            else:
                iV = find_string(varList, header['shortname'])
                if (iV < 0):
                    iV = find_string(varList, header['vars'])
                if (iV < 0):
                    iV = find_string(varList, header['longname'])
                iVars = [iV]
    else:
        iVars = []
        for var in varList:
            if (var.isnumeric()):
                iVars.append(int(var))
            else:
                if (header == None):
                    print('Non number variables are not supported yet!')
                    iVars = [3]
                else:
                    iV = find_string(varList, header['shortname'])
                    if (iV < 0):
                        iV = find_string(varList, header['vars'])
                    if (iV < 0):
                        iV = find_string(varList, header['longname'])
                    iVars.append(iV)

    return iVars

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def read_all_files(filelist, varToPlot):

    filetype = determine_filetype(filelist[0])
    if (filetype["myfile"] == filetype["iGitmBin"]):
        header = gitmio.read_gitm_headers(filelist[0])
        varToPlot = convert_var_to_number(varToPlot, header)

    print('vartoplot : ', varToPlot)
    allData = gitmio.read_gitm_all_files(filelist, varToPlot)
    return allData

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def read_header(filelist):

    filetype = determine_filetype(filelist[0])
    if (filetype["myfile"] == filetype["iGitmBin"]):
        header = gitmio.read_gitm_headers(filelist[0])
    else:
        print('File type not supported at this point!')

    return header



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