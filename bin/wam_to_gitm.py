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
        description = 'Convert WAM netcdf files to GITM output files')
    
    parser.add_argument('-altmin', default = 70.0, type = float, \
                        help = 'altitude : alt minimum in km') 
    parser.add_argument('-altmax', default = 100.0, type = float, \
                        help = 'altitude : alt maximum in km') 
    parser.add_argument('-dalt', default = 2.5, type = float, \
                        help = 'altitude : alt maximum in km') 
    
    # Get the files to plot:
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to convert')

    args = parser.parse_args()

    return args

# Needed to run main script as the default executable from the command line
if __name__ == '__main__':

    # Get the input arguments
    args = get_args()
    filelist = args.filelist

    # wam variables:
    vars = [ 'O', 'O2', 'N2', \
             'u_neutral', 'v_neutral', 'w_neutral', \
             'temp_neutral']
    iO_ = 0
    iO2_ = 1
    iN2_ = 2
    iTn_ = 6
    # map to:
    varsOut = [ 'Longitude', 'Latitude', 'Altitude',
                'Tpert', 'Npert', \
                'O', 'O2', 'N2', \
                'Ve', 'Vn', 'Vv', \
                'Tn' ]
    nVarsOffset = len(varsOut) - len(vars)

    minAlt = args.altmin
    maxAlt = args.altmax
    dAlt = args.dalt
    nAlts = (maxAlt - minAlt) / dAlt + 1
    altCuts = np.arange(minAlt, maxAlt + dAlt, dAlt)
    print(altCuts)

    allData = util.read_all_files(filelist[0], ['O'])
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
    
    nTimes = len(filelist)
    nVars = len(vars)
    allDataOnUniformGrid = np.zeros((nTimes, nVars, nLons, nLats, nAlts))
    for iVar, var in enumerate(vars):
        print('Processing Variable : ', var)
        allData = util.read_all_files(filelist, var)
        varOnUniformGrid = np.zeros((nTimes, nLons, nLats, nAlts))
        for iAlt, altGoal in enumerate(altCuts):
            sliceData = utils.data_slice(allData, targetAlt = altGoal)
            allSlices = sliceData['slices']
            allDataOnUniformGrid[:, iVar, :, :, iAlt] = allSlices

    weights = np.cos(latitudes * np.pi / 180.0)
    totalWeights = np.sum(weights)
    for iTime in range(nTimes):
        gitmData['time'] = allData['times'][iTime]
        for iVar in range(nVars):
            var = varsOut[iVar + nVarsOffset]
            gitmData[var] = allDataOnUniformGrid[iTime, iVar, :, :, :]
        Ntotal = \
            allDataOnUniformGrid[iTime, iO_, :, :, :] + \
            allDataOnUniformGrid[iTime, iO2_, :, :, :] + \
            allDataOnUniformGrid[iTime, iN2_, :, :, :]
        Tn = allDataOnUniformGrid[iTime, iTn_, :, :, :]
        Npert = Ntotal
        Tpert = Tn
        for iAlt, altGoal in enumerate(altCuts):
          meanN = np.mean(Ntotal[:,:,iAlt])
          Npert[:,:,iAlt] = (Ntotal[:,:,iAlt] - meanN) / meanN
          #meanTn = np.sum(Tn * weights) / totalWeights
          meanTn = np.mean(Tn[:,:,iAlt])
          Tpert[:,:,iAlt] = Tn[:,:,iAlt] - meanTn
        gitmData['Tpert'] = Tpert
        gitmData['Npert'] = Npert
        file = gitmData['time'].strftime('3DALL_t%Y%m%d_%H%M%S.bin')
        write_gitm_file(file, gitmData, doReverse = True, isVerbose = True)




