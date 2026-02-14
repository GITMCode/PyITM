#!/usr/bin/env python3

import numpy as np
import sys
import argparse

from pyitm.fileio import gitmio, util
from pyitm.modeldata import utils
from pyitm.general import time_conversion, geometry

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def get_args():

    parser = argparse.ArgumentParser(description = \
                                     'Calculate values and write logfile')

    parser.add_argument('filelist', nargs = '+', \
                        help = 'Files to process')
    parser.add_argument('-ngcs',
                        default = 2, type = int, \
                        help = 'number of ghostcells')
    parser.add_argument('-integral',  \
                        action='store_true', default = False, \
                        help = 'Calculate integral instead of mean')
    
    args = parser.parse_args()

    return args

# Needed to run main script as the default executable from the command line
if __name__ == '__main__':

    # Get the input arguments
    args = get_args()
    filelist = args.filelist
    nGCs = args.ngcs
    doIntegral = args.integral
    
    allData = util.read_all_files(filelist)

    if (not allData):
        util.list_file_info(filelist)
        exit()

    #data = {'times': allTimes,
    #        'data': allData,
    #        'vars': vars,
    #        'shortname': variables.get_short_names(vars),
    #        'longname': variables.get_long_names(vars),
    #        'lons': lons,
    #        'lats': lats,
    #        'alts': alts,
    #        'ntimes': nTimes,
    #        'nvars': nVars,
    #        'nlons' : nLons,
    #        'nlats': nLats,
    #        'nblocks': 0,
    #        'nalts': nAlts}

    nLons = allData['nlons']
    nLats = allData['nlats']
    nAlts = allData['nalts']
    nVars = allData['nvars']
    vars = allData['vars']
    
    area3d = np.zeros((nLons, nLats, nAlts))
    totalarea1d = np.zeros(nAlts)
    lons2d = allData['lons'][:,:,0]
    lats2d = allData['lats'][:,:,0]

    for iAlt in range(nAlts):
        realAlt = allData['alts'][0,0,iAlt]
        area2d = geometry.calc_areas(lons2d, lats2d, realAlt)
        area3d[:,:,iAlt] = area2d
        totalarea1d[iAlt] = np.sum(area2d[nGCs:-nGCs, nGCs:-nGCs])

    dataOut = {
        'version': 2.0,
        'nLonsTotal': 1,
        'nLatsTotal': 1,
        'nAltsTotal': nAlts,
        'nVars': nVars,
        'vars': vars}

    for var in vars:
        dataOut[var] = np.zeros((1,1,nAlts))
    
    for iTime, file in enumerate(filelist):
        fileOut = '1DMEAN_' + file
        dataOut['time'] = allData['times'][iTime]
        
        for iVar, var in enumerate(vars):
            for iAlt in range(nAlts):

                vals2d = allData['data'][iTime, iVar, nGCs:-nGCs, nGCs:-nGCs, iAlt]
                area2d = area3d[nGCs:-nGCs, nGCs:-nGCs, iAlt]
                integral = np.sum(area2d * vals2d)
                if args.integral:
                    dataOut[var][0,0,iAlt] = integral
                else:
                    value = integral / totalarea1d[iAlt]
                    print(var, iAlt, value)
                    dataOut[var][0,0,iAlt] = value

        gitmio.write_gitm_file(fileOut, dataOut, isVerbose = True)
                
        # need to write out:
        #    - version
        #    - nLonsTotal
        #    - nLatsTotal
        #    - nAltsTotal
        #    - nVars
        #    - time
        #    - vars
        #    - all of the variable names (which should be the same as in vars
        #      including lon, lat, alt
        
