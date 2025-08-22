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

    header = util.read_header(filelist)

    print('File information:')
    print(' -> nLats : ', header['nLats'])
    print(' -> nLons : ', header['nLons'])
    print(' -> nAlts : ', header['nAlts'])
    print(' -> nVars : ', header['nVars'])
    for i in range(header['nVars']):
        print('  -> %3d : ' % i, header['vars'][i], '->', 
              header['shortname'][i], '->', header['longname'][i])
