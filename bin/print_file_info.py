#!/usr/bin/env python3
""" Provide information about the file
"""

import argparse

import sys
sys.path.insert(0,'/home/ridley/Software/PyITM/')

from pyitm.fileio import util

def get_args():

    parser = argparse.ArgumentParser(
        description = 'List file information for Aether / GITM model results')
    
    # Get the files to plot:
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')

    args = parser.parse_args()

    return args

# Needed to run main script as the default executable from the command line
if __name__ == '__main__':

    # Get the input arguments
    args = get_args()
    filelist = args.filelist

    header = util.read_all_headers(filelist)
    print(header)

    print('File information:')
    print(' -> nLats : ', header['nlats'])
    print(' -> nLons : ', header['nlons'])
    print(' -> nAlts : ', header['nalts'])
    print(' -> nVars : ', header['nvars'])
    for i in range(header['nvars']):
        print('  -> %3d : ' % i, header['vars'][i], '->', 
              header['shortname'][i], '->', header['longname'][i])
