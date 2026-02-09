#!/usr/bin/env python3

import numpy as np
import sys
import argparse

from pyitm.fileio import logfile, util
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
    parser.add_argument('-logfile',
                        help = 'output logfile front',
                        default = 'log')
    # select altitude to output:
    parser.add_argument('-alt', metavar = 'alt',
                        default = 400.0, type = float, \
                        help = 'altitude :  alt in km (closest)')
    parser.add_argument('-ngcs',
                        default = 2, type = int, \
                        help = 'number of ghostcells')
    parser.add_argument('-mean',  \
                        action='store_true', default = False, \
                        help = 'Calculate mean instead of integral')
    
    args = parser.parse_args()

    return args

# Needed to run main script as the default executable from the command line
if __name__ == '__main__':

    # Get the input arguments
    args = get_args()
    filelist = args.filelist
    nGCs = args.ngcs
    
    allData = util.read_all_files(filelist)

    if (not allData):
        util.list_file_info(filelist)
        exit()

    altGoal = args.alt
    alts1d = allData['alts'][0, 0, :]
    diff = np.abs(alts1d - altGoal)
    iAlt = np.argmin(diff)
    realAlt = alts1d[iAlt]
        
    lons2d = allData['lons'][:, :, iAlt]
    lats2d = allData['lats'][:, :, iAlt]

    area = geometry.calc_areas(lons2d, lats2d, realAlt)
    area = area[nGCs:-nGCs, nGCs:-nGCs]
    if args.mean:
        totalarea = np.sum(area)
        message = 'Geometric global mean'
    else:
        message = 'Geometric global integral'
    vars = allData['vars']
    logdata = {'times': []}
    for var in vars:
        logdata[var] = []
    for iTime, time in enumerate(allData['times']):
        logdata['times'].append(time)
        for iVar, var in enumerate(vars):
            vals2d = allData['data'][iTime, iVar, nGCs:-nGCs, nGCs:-nGCs, iAlt]
            integral = np.sum(area * vals2d)
            if args.mean:
                value = integral / totalarea
                logdata[var].append(value)
            else:
                logdata[var].append(integral)

    logfile.write_log(logdata, fileHeader = 'log', message = message)
                                    
    
