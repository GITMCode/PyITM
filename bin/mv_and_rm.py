#!/usr/bin/env python3

import argparse
import re

from pyitm.fileio import remote
from pyitm.general import system

IsVerbose = False

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args_post():

    parser = argparse.ArgumentParser(
        description = "get files from a remote server",
        formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-remotefile',
                        help = 'File that contains info. about remote system',
                        default = '.remote')

    parser.add_argument('-files',
                        help = 'Files to get ("star" = "*")',
                        default = 'star.png')

    parser.add_argument('-user',
                        help = 'remote user name (default none)',
                        default = 'none')
    
    parser.add_argument('-server',
                        help = 'remote system name (default none)',
                        default = 'none')
    
    parser.add_argument('-dir',
                        help = 'remote directory to use (default none)',
                        default = 'none')
    
    parser.add_argument('-v',
                        help = 'Run with verbose',
                        action = 'store_true')
    
    parser.add_argument('-rm',
                        help = 'remove files after pushing them',
                        action = 'store_true')
    
    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# Main Code
# ----------------------------------------------------------------------

if __name__ == '__main__':  # main code block

    args = parse_args_post()

    # make local variables for arguments:
    files = args.files
    IsVerbose = args.verbose
    DoRemove = args.rm
    
    IsRemote, user, server, dir = load_remote_file(args.remotefile)
    remote = user + '@' + server + ':' + dir + '/'

    push_files(filelist, user, server, dir, DoRemove, IsVerbose = IsVerbose)

