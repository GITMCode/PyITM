#!/usr/bin/env python3

import argparse
import os
import re

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
    
    parser.add_argument('-v', '--verbose',
                        help = 'Run with verbose',
                        action = 'store_true')
    
    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# parse remote file
#   - remote file has the form:
# username
# remote server
# remote directory
# ----------------------------------------------------------------------

def parse_remote_file(file):

    if (IsVerbose):
        print('Reading file ', file)
    
    fpin = open(file, 'r')
    user = fpin.readline()
    server = fpin.readline()
    dir = fpin.readline()
    fpin.close()

    remote = {'user': user.strip(),
              'server': server.strip(),
              'dir': dir.strip()}
    return remote

# ----------------------------------------------------------------------
# do system command
# ----------------------------------------------------------------------

def run_command(command):
    if (IsVerbose):
        print("   -> Running Command : ")
        print("      ", command)
    os.system(command)
    return True

# ----------------------------------------------------------------------
# Check inputs:
# ----------------------------------------------------------------------

def check_inputs(user, server, dir):
    
    IsRemote = True
    if ((len(user) == 0) or (user == 'none')):
        if (IsVerbose):
            print("Can't parse user information")
        IsRemote = False
    if ((len(server) == 0) or (server == 'none')):
        if (IsVerbose):
            print("Can't parse server information")
        IsRemote = False
    if ((len(dir) == 0) or (dir == 'none')):
        if (IsVerbose):
            print("Can't parse dir information")
        IsRemote = False

    return IsRemote

# ----------------------------------------------------------------------
# load remote file
# ----------------------------------------------------------------------

def load_remote_file(args):

    remoteFile = args.remotefile
    
    # figure out remote system information:
    if (os.path.exists(remoteFile)):
        print('Found file: ', remoteFile)
        remote = parse_remote_file(remoteFile)
        user = remote['user']
        server = remote['server']
        dir = remote['dir']
    else:
        user = args.user
        server = args.server
        dir = args.dir

    # Check remote system inputs:
    IsRemote = check_inputs(user, server, dir)
    return IsRemote, user, server, dir
    
# ----------------------------------------------------------------------
# Main Code
# ----------------------------------------------------------------------

if __name__ == '__main__':  # main code block

    args = parse_args_post()

    # make local variables for arguments:
    files = args.files
    IsVerbose = args.verbose

    IsRemote, user, server, dir = load_remote_file(args)
    remote = user + '@' + server + ':' + dir + '/'

    m = re.search('(.*)(star)(.*)', files)
    if m:
        files = '"' + m.group(1) + '*' + m.group(3) + '"'
    remote = remote + files
    command = 'rsync -v ' + remote + ' .'
    run_command(command)
