#!/usr/bin/env python3

import os
from datetime import datetime
from struct import unpack
import numpy as np
from pyitm.fileio import util 

#-----------------------------------------------------------------------------

def read_gitm_single_header(file):
    r""" Grab ancillary information from the GITM file

    Parameters
    ----------
    file - name of the file to read and get header from

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
    
    header = {"nFiles": len(filelist), \
              "version": 0.1, #TODO
              "nLons": 0, \
              "nLats": 0, \
              "nAlts": 0, \
              "nVars": 0, \
              "vars": [], \
              "time": [], \
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
        (header["nLons"],header["nLats"],header["nAlts"]) = \
            unpack(endChar+'lll',f.read(recLen))
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Read number of variables.
        header["nVars"]=unpack(endChar+'l',f.read(recLen))[0]
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Collect variable names.
        for i in range(header["nVars"]):
            v = unpack(endChar+'%is'%(recLen),f.read(recLen))[0]
            if (file == filelist[0]):
                header["vars"].append(v.decode('utf-8').replace(" ",""))
            (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Extract time. 
        (yy,mm,dd,hh,mn,ss,ms)=unpack(endChar+'lllllll',f.read(recLen))
        header["time"].append(datetime(yy,mm,dd,hh,mn,ss,ms*1000))

        f.close()

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
    data["time"]: datetime of the file
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
            "nLons": 0, \
            "nLats": 0, \
            "nAlts": 0, \
            "nVars": 0, \
            "time": 0, \
            "vars": []}

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
        data["time"] = datetime(yy,mm,dd,hh,mn,ss,ms*1000)
    
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
            data[iVar] = np.array(unpack(endChar+'%id'%(nTotal),f.read(s)))
            data[iVar] = data[iVar].reshape( 
                (data["nLons"],data["nLats"],data["nAlts"]),order="F")

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
    spatialData = read_gitm_one_file(filelist[0], vars, verbose=False)

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
        nVars = read_gitm_headers(filelist[0], verbose=False)['nVars']
        varlist = list(range(nVars))

    if (nVars == 1):
        allData = np.zeros((nTimes, nLons, nLats, nAlts))
    else:
        allData = np.zeros((nTimes, nVars, nLons, nLats, nAlts))

    allTimes = []
    for iTime, filename in enumerate(filelist):
        data = read_gitm_one_file(filename, varlist, verbose=verbose)
        allTimes.append(data["time"])
        if (nVars == 1):
            allData[iTime, :, :, :] = data[varlist[0]][:, :, :]
        else:
            for iVar, var in enumerate(varlist):
                allData[iTime, iVar, :, :, :] = data[var][:, :, :]
    vars = []
    for var in varlist:
        vars.append(data['vars'][var])

    data = {'times': allTimes,
            'data': allData,
            'vars': vars,
            'lons': lons,
            'lats': lats,
            'alts': alts,
            'nTimes': nTimes,
            'nVars': nVars,
            'nLons' : nLons,
            'nLats': nLats,
            'nAlts': nAlts}
    
    return data


#-----------------------------------------------------------------------------
# locations should contain :
#   "lats" -> latitudes in degrees
#   "lons" -> longitudes in degrees
#   "alts" -> altitudes in km

def gitm_extract_data(locations, times, variables, dir = './'):

    headers = read_gitm_headers(dir + '3DALL')

    file = headers["filename"][0]

    outdata = {}
    for v in variables:
        outdata[v] = []
    
    # get lon/lat/alt:
    vars = [0,1,2]
    print("Reading first GITM file to get location info")
    data = read_gitm_one_file(file, vars)

    Alts = data[2][0][0]/1000.0
    Lons = np.rad2deg(data[0][:,0,0])
    Lats = np.rad2deg(data[1][0,:,0])

    minLat = np.min(Lats)
    minLon = np.min(Lons)

    dLat = Lats[1]-Lats[0]
    dLon = Lons[1]-Lons[0]

    nFiles = len(headers["time"])
    nPts = len(locations["lats"])

    SatLats = locations["lats"]
    SatLons = locations["lons"]
    SatAlts = locations["alts"]

    nAlts = len(Alts)

    iSatLats = (SatLats-minLat)/dLat
    rSatLats = np.array(iSatLats - np.floor(iSatLats))
    iSatLats = (np.array((np.floor(iSatLats)))).astype(int)

    iSatLons = (SatLons-minLon)/dLon
    rSatLons = iSatLons - np.floor(iSatLons)
    iSatLons = (np.array((np.floor(iSatLons)))).astype(int)

    # We need to figure out Alts, which are not uniform:
    
    SatMinAlt = np.min(SatAlts)
    iSatMinAlt = 0
    while (SatMinAlt > Alts[iSatMinAlt]):
        iSatMinAlt = iSatMinAlt + 1
    if (iSatMinAlt > 0):
        iSatMinAlt = iSatMinAlt-1

    iSatAlts = []
    rSatAlts = []
    for iPt in range(nPts):

        if (SatAlts[iPt] < Alts[0]):
            iSatAlts.append(0)
            rSatAlts.append(0.0)
        else:
            if (SatAlts[iPt] > Alts[nAlts-1]):
                iSatAlts.append(nAlts-2)
                rSatAlts.append(1.0)
            else:
                i = iSatMinAlt
                while (SatAlts[iPt] > Alts[i]):
                    i = i + 1
                i = i - 1
                r = (SatAlts[iPt] - Alts[i])/(Alts[i+1]-Alts[i])
                iSatAlts.append(i)
                rSatAlts.append(r)


    print("Reading files to extract data...")
    for iPt in range(nPts):

        iFileSave = -1
        for iFile in range(nFiles):
            if (times[iPt] >= headers["time"][iFile]):
                iFileSave = iFile
        if (iFileSave == nFiles-1 or iFileSave == -1):
            i = nFiles-2
            r = 1.0
        else:
            i = iFileSave
            r = (times[iPt] - headers["time"][i]).total_seconds() / (headers["time"][i+1] - headers["time"][i]).total_seconds()

        if (iPt == 0):
            iFileLeft = i
            iFileRight = i+1
            fileLeft = headers["filename"][iFileLeft]
            dataLeft = read_gitm_one_file(fileLeft, variables)
            fileRight = headers["filename"][iFileRight]
            dataRight = read_gitm_one_file(fileRight, variables)

        if (iFileLeft != i):
            if (iFileRight == i):
                # Move right to left and grab new data:
                fileLeft = fileRight
                dataLeft = dataRight
                iFileLeft = iFileRight
                iFileRight = i+1
                fileRight = headers["filename"][iFileRight]
                dataRight = read_gitm_one_file(fileRight, variables)
            else:
                # may have skipped a few files or something, so reset:
                iFileLeft = i
                iFileRight = i+1
                fileRight = headers["filename"][iFileRight]
                dataRight = read_gitm_one_file(fileRight, variables)
                fileLeft = headers["filename"][iFileLeft]
                dataLeft = read_gitm_one_file(fileLeft, variables)
                print("Skipping to a new set of files : ", fileLeft, fileRight)

        # At this point, the data should be good, so we do interpolation:
        
        # First get lat, lon, alt for the two data sets:
        
        iLon = iSatLons[iPt]
        rLon = rSatLons[iPt]
        iLat = iSatLats[iPt]
        rLat = rSatLats[iPt]

        iAlt = iSatAlts[iPt]
        rAlt = rSatAlts[iPt]

        data = []

        for iVar in variables:
            dataLeft0 = ( (1.0-rLon) * (1.0-rLat) * (1.0-rAlt) * dataLeft[iVar][iLon  ][iLat  ][iAlt] +
                          (1.0-rLon) * (1.0-rLat) * (    rAlt) * dataLeft[iVar][iLon  ][iLat  ][iAlt+1] +
                          (1.0-rLon) * (    rLat) * (1.0-rAlt) * dataLeft[iVar][iLon  ][iLat+1][iAlt] +
                          (1.0-rLon) * (    rLat) * (    rAlt) * dataLeft[iVar][iLon  ][iLat+1][iAlt+1] +
                          (    rLon) * (1.0-rLat) * (1.0-rAlt) * dataLeft[iVar][iLon+1][iLat  ][iAlt] +
                          (    rLon) * (1.0-rLat) * (    rAlt) * dataLeft[iVar][iLon+1][iLat  ][iAlt+1] +
                          (    rLon) * (    rLat) * (1.0-rAlt) * dataLeft[iVar][iLon+1][iLat+1][iAlt] +
                          (    rLon) * (    rLat) * (    rAlt) * dataLeft[iVar][iLon+1][iLat+1][iAlt+1] )

            dataRight0 = ( (1.0-rLon) * (1.0-rLat) * (1.0-rAlt) * dataRight[iVar][iLon  ][iLat  ][iAlt] +
                           (1.0-rLon) * (1.0-rLat) * (    rAlt) * dataRight[iVar][iLon  ][iLat  ][iAlt+1] +
                           (1.0-rLon) * (    rLat) * (1.0-rAlt) * dataRight[iVar][iLon  ][iLat+1][iAlt] +
                           (1.0-rLon) * (    rLat) * (    rAlt) * dataRight[iVar][iLon  ][iLat+1][iAlt+1] +
                           (    rLon) * (1.0-rLat) * (1.0-rAlt) * dataRight[iVar][iLon+1][iLat  ][iAlt] +
                           (    rLon) * (1.0-rLat) * (    rAlt) * dataRight[iVar][iLon+1][iLat  ][iAlt+1] +
                           (    rLon) * (    rLat) * (1.0-rAlt) * dataRight[iVar][iLon+1][iLat+1][iAlt] +
                           (    rLon) * (    rLat) * (    rAlt) * dataRight[iVar][iLon+1][iLat+1][iAlt+1] )

            outdata[iVar].append((1.0-r) * dataLeft0 + r * dataRight0)
            if (iVar < 2):
                outdata[iVar][-1] = outdata[iVar][-1] * rtod
            if (iVar == 2):
                outdata[iVar][-1] = outdata[iVar][-1] / 1000.0


    return outdata

