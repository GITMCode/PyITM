#!/usr/bin/env python

from datetime import datetime
import numpy as np
import os, glob, re

def read_one_fism_file(fismFile, verbose=False):

    if not os.path.exists(fismFile):
        raise FileNotFoundError(f"Could not find logfile: {fismFile}")
    
    print('> Reading logfile: ', fismFile)

    dataOut = {'times': []}
    nWaves = 59
    with open(fismFile, 'r') as fpin:
        allLines = fpin.readlines()
        nLines = len(allLines)
        # assume first line is #START
        irradiances = np.zeros((nLines-1, nWaves))
        for iLine in range(1, nLines):
            values = allLines[iLine].strip().split()
            t = datetime(int(values[0]),
                         int(values[1]),
                         int(values[2]),
                         int(values[3]),
                         int(values[4]),
                         int(values[5]))
            dataOut['times'].append(t)
            irradiances[iLine - 1, :] = np.array(values[6:]).astype(float)
    dataOut['irradiances'] = irradiances
    return dataOut

def merge_fism_dicts(fismDict1, fismDict2):
    nT1, nWaves = np.shape(fismDict1['irradiances'])
    nT2 = len(fismDict2['times'])
    irradiances = np.zeros((nT1 + nT2, nWaves))
    for iT in range(nT1):
        irradiances[iT, :] = fismDict1['irradiances'][iT,:]
    for iT in range(nT2):
        irradiances[nT1 + iT, :] = fismDict2['irradiances'][iT,:]
    fismDictOut = {'times': np.concatenate((fismDict1['times'], \
                                            fismDict2['times'])), \
                       'irradiances': irradiances}    
    return fismDictOut

def read_all_fism_files(filelist, verbose=False):

    for iFile, file in enumerate(filelist):
        if (iFile == 0):
             fismData = read_one_fism_file(file, verbose = verbose)
        else:
             fismDataNext = read_one_fism_file(file, verbose = verbose)
             print(' -> merging')
             fismData = merge_fism_dicts(fismData, fismDataNext)
    return fismData

def read_f107_file(f107File):
    
    if not os.path.exists(f107File):
        raise FileNotFoundError(f"Could not find logfile: {f107File}")
    
    print('> Reading f107 file: ', f107File)

    with open(f107File, 'r') as fpin:
        allLines = fpin.readlines()
        nLines = len(allLines)

        iOff = 0
        isFound = False
        f107s = []
        times = []
        while (not isFound):
            m = re.match(r'#yyyy-MM-dd', allLines[iOff])
            if m:
                isFound = True
            iOff = iOff + 1
        for iLine in range(iOff, nLines):
            values = allLines[iLine].strip().split()
            ymd = values[0].split('-')
            hm = values[1].split(':')
            t = datetime(int(ymd[0]),
                         int(ymd[1]),
                         int(ymd[2]),
                         int(hm[0]),
                         int(hm[1]),
                         0)
            times.append(t)
            f107s.append(float(values[2]))

    f107 = np.array(f107s)
    nTimes = len(times)
    floatTimes = np.zeros(nTimes)
    for iT in range(nTimes):
        floatTimes[iT] = (times[iT] - times[0]).total_seconds()
    f107a = np.zeros(nTimes)
    halfWindow = 40.5 * 86400.0
    for iT in range(nTimes):
        dt = np.abs(floatTimes - floatTimes[iT])
        f107a[iT] = np.mean(f107[dt < halfWindow])
        
    f107Data = {'f107': f107,
                'f107a': f107a,
                'times': np.array(times),
                'floatTimes': floatTimes}

    return f107Data



    
