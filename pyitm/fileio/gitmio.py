#!/usr/bin/env python3

import os, glob
from datetime import datetime
from struct import unpack
import numpy as np
from pyitm.fileio import util
from pyitm.fileio import variables

#-----------------------------------------------------------------------------

def read_gitm_single_header(file, verbose=True):
    r""" Grab ancillary information from the GITM file

    Parameters
    ----------
    file - (str/path) name of the file to read and get header from
    verbose - (bool) enable verbose outputs? True by default.

    Returns
    -------
    header: A dictionary containing information about the netcdf file, such
            as nLons, nLons, nAlts, nVars, variable names, time(s)

    Notes
    -----
    - This is a wrapper for read_gitm_headers that will read a single header file

    """
    
    # check if file exists. If yes, just read. if not, attempt to parse
    if os.path.exists(file):
        filename = file
    else:
        filename = util.any_to_filelist(file) # This is a list
        if len(filename) > 1:
            print("Received multiple headers! Only reading the first!!")
        filename = filename[0]

    header = read_gitm_headers(filename, verbose=verbose)

    return header

#-----------------------------------------------------------------------------


def read_gitm_headers(input_files='./3DALL*.bin', verbose=False):
    r""" Grab ancillary information from GITM output files

    Parameters
    ----------
    input_files: (str or list-like) - pattern to glob or list of files to read.
        Default './3DALL*.bin', so only looks in `pwd`
        - Examples: 
            -> "/path/to/GITM/run/data/3DALL*.bin"
            -> "/path/to/GITM/run/data/3DALL" 
            -> ["/path/to/GITM/run/data/3DALL_t110926_000000.bin", 
                "/path/to/GITM/run/data/3DALL_t110926_001000.bin", 
                "/path/to/GITM/run/data/3DALL_t110926_002000.bin"]

    Returns
    -------
    header: A dictionary containing information about the netcdf file, such
            as nLons, nLons, nAlts, nVars, variable names, time(s)

    """

    filelist = util.any_to_filelist(input_files)

    # sanity check to make sure some files exist:
    if verbose:
        print("  -> Found ", len(filelist), "files")
    
    header = {"ntimes": len(filelist), \
              "version": 0.1, #TODO
              "nlons": 0, \
              "nlats": 0, \
              "nalts": 0, \
              "nblocks": 0, \
              "nvars": 0, \
              "vars": [], \
              "shortname": [], \
              "longname": [], \
              "times": [], \
              "filename": [] }

    for file in filelist:

        header["filename"].append(file)

        f=open(file, 'rb')

        # This is all reading header stuff:

        endChar='>'
        rawRecLen=f.read(4)
        recLen=(unpack(endChar+'l',rawRecLen))[0]
        if (recLen>10000)or(recLen<0):
            # Ridiculous record length implies wrong endian.
            endChar='<'
            recLen=(unpack(endChar+'l',rawRecLen))[0]

        # Read version; read fortran footer+header.
        header["version"] = unpack(endChar+'d',f.read(recLen))[0]

        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Read grid size information.
        (header["nlons"],header["nlats"],header["nalts"]) = \
            unpack(endChar+'lll',f.read(recLen))
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Read number of variables.
        header["nvars"]=unpack(endChar+'l',f.read(recLen))[0]
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Collect variable names.
        for i in range(header["nvars"]):
            v = unpack(endChar+'%is'%(recLen),f.read(recLen))[0]
            if (file == filelist[0]):
                header["vars"].append(v.decode('utf-8').replace(" ",""))
            (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Extract time. 
        (yy,mm,dd,hh,mn,ss,ms)=unpack(endChar+'lllllll',f.read(recLen))
        header["times"].append(datetime(yy,mm,dd,hh,mn,ss,ms*1000))

        f.close()

    header["shortname"] = variables.get_short_names(header["vars"])
    header["longname"] = variables.get_long_names(header["vars"])

    return header


#-----------------------------------------------------------------------------


def read_gitm_one_file(file_to_read, varlist=[-1], verbose=True):
    r""" Read list of variables from one GITM file

    Parameters
    ----------
    file_to_read: GITM file to read
    varlist: list of variable NUMBERS to read. Use [-1] (default) for all variables.
    verbosr: bool - print extra info when running? Default = True

    Returns
    -------
    data["times"]: datetime of the file
    data[NUMBER]: data that is read in.
                  NUMBER goes from 0 - number of vars read in (0-3 typical)
    (Also include header information, as described above)
    """

    # Double check input arguments
    if isinstance(varlist, (str, int, float)):
        raise TypeError("read_gitm_one_file must be called with a list of ints!\n"
                        "Maybe try using the default value of [-1] instead.")

    if '~' in file_to_read:
        file_to_read = os.path.expanduser(file_to_read)

    if verbose:
        print('-> Reading file : ' + file_to_read, ' --> Vars : ', varlist)

    data = {"version": 0, \
            "nlons": 0, \
            "nlats": 0, \
            "nalts": 0, \
            "nblocks": 0, \
            "nvars": 0, \
            "times": 0, \
            "vars": [],
            "shortname": [], \
            "longname": [], \
            "data": {},
            }

    print(' -> Reading GITM file: ', file_to_read)
    with open(file_to_read, 'rb') as f:
    
        # This is all reading header stuff:
        endChar='>'
        rawRecLen=f.read(4)
        rec_len = (unpack(endChar + 'l', rawRecLen))[0]
        if rec_len > 10000 or rec_len < 0:
            # Ridiculous record length implies wrong endian.
            endChar='<'
            rec_len = (unpack(endChar + 'l', rawRecLen))[0]
    
        # Read version; read fortran footer+data.
        data["version"] = unpack(endChar + 'd', f.read(rec_len))[0]
    
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))
    
        # Read grid size information.
        (data["nLons"],data["nLats"],data["nAlts"]) = \
            unpack(endChar+'lll',f.read(recLen))
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))
    
        # Read number of variables.
        data["nVars"]=unpack(endChar+'l',f.read(recLen))[0]
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))
    
        if (varlist[0] == -1):
            varlist = np.arange(data['nVars'])
    
        # Collect variable names.
        for i in range(data["nVars"]):
            v = unpack(endChar+'%is'%(recLen),f.read(recLen))[0]
            data["vars"].append(v.decode('utf-8').replace(" ",""))
            (oldLen, recLen)=unpack(endChar+'2l',f.read(8))
    
        # Extract time. 
        (yy,mm,dd,hh,mn,ss,ms)=unpack(endChar+'lllllll',f.read(recLen))
        data["times"] = datetime(yy,mm,dd,hh,mn,ss,ms*1000)
    
        # Header is this length:
        # Version + start/stop byte
        # nLons, nLats, nAlts + start/stop byte
        # nVars + start/stop byte
        # Variable Names + start/stop byte
        # time + start/stop byte
    
        iHeaderLength = 84 + data["nVars"] * 48
    
        nTotal = data["nLons"]*data["nLats"]*data["nAlts"]
        iDataLength = nTotal*8 + 4+4

        for iVar in varlist:
            f.seek(iHeaderLength+iVar*iDataLength)
            s=unpack(endChar+'l',f.read(4))[0]
            data["data"][iVar] = np.array(unpack(endChar+'%id'%(nTotal),f.read(s)))
            data["data"][iVar] = data["data"][iVar].reshape( 
                (data["nLons"],data["nLats"],data["nAlts"]),order="F")

    data["shortname"] = variables.get_short_names(data["vars"])
    data["longname"] = variables.get_long_names(data["vars"])

    return data

# -----------------------------------------------------------------------------
# This reads in a series of vars / files and returns the 3D information
#-----------------------------------------------------------------------------

def read_gitm_all_files(filelist, varlist=[-1], verbose=False):

    filelist = util.any_to_filelist(filelist)

    # Get the prefixes for all entries in filelist;
    prefixes = np.unique([file.split('/')[-1].split('_')[0] for file in filelist])
    if len(prefixes) > 1: # make sure there is only one output type.
        raise ValueError("Multiple output types cannot be read by this function."
                        "\n\tProvided: " + str(prefixes))

    # first read in spatial information:
    vars = [0, 1, 2]
    spatialData = read_gitm_one_file(filelist[0], vars, verbose=False)["data"]

    lons = np.degrees(spatialData[0])  # Convert from rad to deg
    nLons = len(lons[:, 0, 0])
    lats = np.degrees(spatialData[1])  # Convert from rad to deg
    nLats = len(lats[0, :, 0])
    alts = spatialData[2] / 1000.0  # Convert from m to km
    nAlts = len(alts[0, 0, :])

    nTimes = len(filelist)
    if varlist != [-1]:
       nVars = len(varlist)
    else: # varlist=[-1] means we read in all variables
        nVars = read_gitm_headers(filelist[0], verbose=False)['nvars']
        varlist = list(range(nVars))

    if (nVars == 1):
        allData = np.zeros((nTimes, nLons, nLats, nAlts))
    else:
        allData = np.zeros((nTimes, nVars, nLons, nLats, nAlts))

    allTimes = []
    for iTime, filename in enumerate(filelist):
        data = read_gitm_one_file(filename, varlist, verbose=verbose)
        allTimes.append(data["times"])
        for iVar, var in enumerate(varlist):
            if (nVars == 1):
                allData[iTime, :, :, :] = data["data"][var][:, :, :]
            else:
                allData[iTime, iVar, :, :, :] = data["data"][var][:, :, :]

    vars = []
    for var in varlist:
        vars.append(data['vars'][var])

    data = {'times': allTimes,
            'data': allData,
            'vars': vars,
            'shortname': variables.get_short_names(vars),
            'longname': variables.get_long_names(vars),
            'lons': lons,
            'lats': lats,
            'alts': alts,
            'ntimes': nTimes,
            'nvars': nVars,
            'nlons' : nLons,
            'nlats': nLats,
            'nblocks': 0,
            'nalts': nAlts}
    if verbose:
        print('data keys : ', data.keys())
    return data

#-----------------------------------------------------------------------------

def read_logfile(logfile=None, datadir=None, verbose=False):
    """ Read GITM logfile to a dictionary. Either datapath or logfile must be provided.
    
    Parameters
    ----------
    logfile: (str/path) - path to logfile to read. Default None.
    datadir: (str/path) - path to GITM data directory. Default None.
    verbose: (bool) - print extra info when running? Default = False
    
    Returns
    -------
    logdata: (dict) - dictionary containing logfile information

    Notes
    -----
    - If datadir is provided, will look for the most recent logfile in that directory.
    - This should be the directory in which model outputs are stored.
    """

    if datadir: # Grab the most highest numbered logfile
        #TODO: Read & appead all logfiles
        if logfile:
            print("Both logfile and datadir provided. Ignoring logfile and using datadir.")
        logfile = sorted(glob.glob(os.path.join(datadir, 'log*.dat')))[-1]

    if not os.path.exists(logfile):
        raise FileNotFoundError(f"Could not find logfile: {logfile}")
    elif verbose:
        print(f"Reading logfile: {logfile}")

    logdata = {}

    saving = False
    with open(logfile, 'r') as f:
        for n, line in enumerate(f.readlines()):
            if saving:
                # first line (after #START) has the headers
                if len(logdata.keys()) == 0:
                    vars = []
                    for col in line.strip().split():
                        while (col in vars):
                            col = col+'0'
                        vars.append(col)
                        logdata[col] = []
                    if verbose:
                        print(f"  -> Reading columns: {list(logdata.keys())}")
                else:
                    # Read in the data
                    for col, val in zip(logdata.keys(), line.strip().split()):
                        logdata[col].append(val)

            if line.startswith('#START'):
                saving=True
                if verbose:
                    print(f"  -> Found start of data at line {n}")

    if verbose:
        print(f"  -> Read {len(logdata[list(logdata.keys())[0]])} lines of data")

    for col in logdata.keys():
        try:
            logdata[col] = np.array(logdata[col], dtype=float)
        except:
            if verbose:
                print(f"  -> Could not convert column '{col}' to float. Leaving as string.")
            pass

    return logdata
