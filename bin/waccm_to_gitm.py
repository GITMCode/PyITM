#!/usr/bin/env python3

import numpy as np

import sys
sys.path.insert(0,'/home/ridley/Software/GITM/srcPython/')

from pGITM import *
from pyitm.fileio import util
from pyitm.modeldata import utils

import argparse

# ----------------------------------------------------------------------------
# Get arguments as inputs into the code
#-----------------------------------------------------------------------------

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Convert WACCM netcdf files to GITM output files')
    
    parser.add_argument('-v',  \
                        action='store_true', default = False, \
                        help = 'set verbose')
    
    parser.add_argument('-altmin', default = 70.0, type = float, \
                        help = 'altitude : alt minimum in km') 
    parser.add_argument('-altmax', default = 100.0, type = float, \
                        help = 'altitude : alt maximum in km') 
    parser.add_argument('-dalt', default = 2.5, type = float, \
                        help = 'altitude : alt maximum in km') 
    
    # Get the files to plot:
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to convert (only one accepted!)')

    parser.add_argument('-iStart',  default = 0, type = int, \
                        help = 'start index for time in multi-time file')
    parser.add_argument('-iEnd',  default = 0, type = int, \
                        help = 'end index for time in multi-time file')
    parser.add_argument('-iStep',  default = 1, type = int, \
                        help = 'step size for time in multi-time file')

    args = parser.parse_args()

    return args

# Needed to run main script as the default executable from the command line
if __name__ == '__main__':

    # Get the input arguments
    args = get_args()
    isVerbose = args.v
    filelist = args.filelist
    iStart = args.iStart
    iEnd = args.iEnd
    iStep = args.iStep

    if (iEnd < 0):
        header = util.read_all_headers(filelist, verbose = isVerbose)
        iEnd = len(header['times'])-1
        if (isVerbose):
            print('Manually changing iEnd to ', iEnd)
    
    Ck = 1.380649e-23
    
    # wam variables:
    vars = [ 'U', 'V', 'T']
    iU_ = 0
    iV_ = 1
    iTn_ = 2
    # map to:
    varsOut = [ 'Longitude', 'Latitude', 'Altitude',
                'Tpert', 'Npert', \
                'Pressure', \
                'Ve', 'Vn', 'Tn' ]
    nVarsOffset = len(varsOut) - len(vars)

    minAlt = args.altmin
    maxAlt = args.altmax
    dAlt = args.dalt
    nAlts = (maxAlt - minAlt) / dAlt + 1
    altCuts = np.arange(minAlt, maxAlt + dAlt, dAlt)

    # Read in dummy variable (one time):
    allData = util.read_all_files(filelist[0], ['O'], verbose = isVerbose)

    longitudesOrig = allData['lons']
    latitudesOrig = allData['lats']
    nLons, nLats, nAlts = np.shape(longitudesOrig)
    nAlts = len(altCuts)
    latitudes = np.zeros((nLons, nLats, nAlts))
    longitudes = np.zeros((nLons, nLats, nAlts))
    altitudes = np.zeros((nLons, nLats, nAlts))
    
    for iAlt, altGoal in enumerate(altCuts):
        longitudes[:, :, iAlt] = longitudesOrig[:, :, 0]
        latitudes[:, :, iAlt] = latitudesOrig[:, :, 0] 
        altitudes[:, :, iAlt] = altGoal
    gitmData = {'nLonsTotal': nLons, \
                'nLatsTotal': nLats, \
                'nAltsTotal': nAlts, \
                'version': 999.0, \
                'nVars': len(varsOut), \
                'vars': varsOut, \
                'Longitude': longitudes * np.pi / 180.0, \
                'Latitude': latitudes * np.pi / 180.0, \
                'Altitude': altitudes * 1000.0}
    
    iTimes = list(range(iStart, iEnd+1, iStep))
    nTimes = len(iTimes)
    nVars = len(vars)
    allDataOnUniformGrid = np.zeros((nTimes, nVars, nLons, nLats, nAlts))
    for iVar, var in enumerate(vars):
        print('Processing Variable : ', var)
        allData = util.read_all_files(filelist, var,
                                      iStart = args.iStart, \
                                      iEnd = args.iEnd, \
                                      iStep = args.iStep)
        varOnUniformGrid = np.zeros((nTimes, nLons, nLats, nAlts))
        for iAlt, altGoal in enumerate(altCuts):
            if (isVerbose):
                print('  -> alt: ', altGoal)
            sliceData = utils.data_slice(allData, targetAlt = altGoal)
            allSlices = sliceData['slices']
            allDataOnUniformGrid[:, iVar, :, :, iAlt] = allSlices

    # need to get pressure on the altitude grid:
    print(' --> Mapping pressure to altitudes')
    pressure = allData['pressure']
    allData['vars'] = ['pressure']
    allData['data'] = []
    for iTime in range(nTimes):
        allData['data'].append(pressure)
    allData['data'] = np.array(allData['data'])
    pressureOnUniformGrid = np.zeros((nTimes, nLons, nLats, nAlts))
    for iAlt, altGoal in enumerate(altCuts):
        if (isVerbose):
            print('  -> alt: ', altGoal)
        sliceData = utils.data_slice(allData, targetAlt = altGoal)
        allSlices = sliceData['slices']
        pressureOnUniformGrid[:, :, :, iAlt] = allSlices
        
    weights = np.cos(latitudes * np.pi / 180.0)
    totalWeights = np.sum(weights)
    for iTime in range(nTimes):
        gitmData['time'] = allData['times'][iTime]
        for iVar in range(nVars):
            var = varsOut[iVar + nVarsOffset]
            gitmData[var] = allDataOnUniformGrid[iTime, iVar, :, :, :]
        Tn = allDataOnUniformGrid[iTime, iTn_, :, :, :]
        Pressure = pressureOnUniformGrid[iTime, :, :, :]
        Ntotal = Pressure / (Tn * Ck)
        Npert = Ntotal
        Tpert = Tn
        for iAlt, altGoal in enumerate(altCuts):
          meanN = np.mean(Ntotal[:,:,iAlt])
          Npert[:,:,iAlt] = (Ntotal[:,:,iAlt] - meanN) / meanN
          meanTn = np.mean(Tn[:,:,iAlt])
          Tpert[:,:,iAlt] = Tn[:,:,iAlt] - meanTn
        gitmData['Tpert'] = Tpert
        gitmData['Npert'] = Npert
        gitmData['Pressure'] = Pressure
        file = gitmData['time'].strftime('3DALL_t%Y%m%d_%H%M%S.bin')
        write_gitm_file(file, gitmData, doReverse = True, isVerbose = True)




