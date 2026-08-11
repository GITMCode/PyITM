#!/usr/bin/env python3

import os
import re
import time

from pyitm.general import system

# ----------------------------------------------------------------------
# parse remote file
#   - remote file has the form:
# username
# remote server
# remote directory
# ----------------------------------------------------------------------

def parse_remote_file(file, IsVerbose = False):

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
# Check inputs:
# ----------------------------------------------------------------------

def check_inputs(user, server, dir, IsVerbose = False):
    
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

def load_remote_file(remoteFile, IsVerbose = False):

    # figure out remote system information:
    if (os.path.exists(remoteFile)):
        if (IsVerbose):
            print('Found file: ', remoteFile)
        remote = parse_remote_file(remoteFile, IsVerbose = IsVerbose)
        user = remote['user']
        server = remote['server']
        dir = remote['dir']
        # Check remote system inputs:
        IsRemote = check_inputs(user, server, dir, IsVerbose = IsVerbose)
        if (IsVerbose):
            print(' -> Checking if remote file is good : ', IsRemote)
    else:
        IsRemote = False
        
    return IsRemote, user, server, dir

# ----------------------------------------------------------------------
# Checks to see if remote file or directory exists
# ----------------------------------------------------------------------

def test_if_remote_exists(user, server, dir, IsVerbose = False):

    DidWork = True
    
    remote_command = user + "@" + server + " 'ls " + dir + "'"
    # check to see if the remote directory exists:

    if (IsVerbose):
        print(' -> Checking to see if remote file or directory exists')
    command = 'ssh ' + remote_command + ' >& .test_file'
    DidWork = run_command(command)
    DidWork = parse_test_file('.test_file')

    if (DidWork):
        if (IsVerbose):
            print('   --> Remote directory (or file) exists!')
    else:
        print('--> Remote directory (or file) does NOT exist!')
        print('    Need to make this directory!')
        
    return DidWork

# ----------------------------------------------------------------------
# Make a remote directory
# ----------------------------------------------------------------------

def make_remote_dir(user, server, dir):

    DidWork = True
    
    remote_command = user + "@" + server + " 'mkdir " + dir + "'"
    # check to see if the remote directory exists:

    print('Making remote directory : ', dir)
    command = 'ssh ' + remote_command + ' >& .test_file'
    DidWork = run_command(command)
    DidWork = parse_test_file('.mkdir_command')

    return DidWork

# ----------------------------------------------------------------------
# Push files to remote, check if they made it,
#    then delete local (if requested)
# ----------------------------------------------------------------------

def push_files(filelist, user, server, dir, DoRemove, IsVerbose = False):

    DidWork = True
    
    remote = user + '@' + server + ':' + dir

    files = ''
    outfile = '.output_rsync_log'

    if (len(filelist) > 0):
        for file in filelist:
            chmod = 'chmod a+r ' + file
            DidWork = run_command(chmod)
            files = files + ' ' + file
            
        rsync = 'rsync -rav ' + files + ' ' + remote
        if (not IsVerbose):
            rsync = rsync + ' >> ' + outfile + ' 2>&1'
        DidWork = run_command(rsync)
        
    if (DoRemove):
        for file in filelist:
            sep = file.split('/')
            test_file = sep[-1]
            DidTransfer = test_if_remote_exists(user,
                                                server,
                                                dir + '/' + test_file)
            if (DidTransfer):
                if (IsVerbose):
                    print('   --> Remote file (' + test_file + ') exists!' +
                          '  Deleting local!')
                if (DoRm):
                    command = '/bin/rm -f ' + file
                    DidWork = run_command(command)
                    # Systems reject ssh commands if too many happen in
                    # too short of time, so sleep in between the commands. 
                    time.sleep(5)
            else:
                if (IsVerbose):
                    print('Remote file (' + test_file + ') does not exist!')

    return DidWork

# ----------------------------------------------------------------------
# Pull files from remote, check if they made it,
#    then delete remote files (if requested)
# ----------------------------------------------------------------------

def pull_files(files, user, server, dir, DoRemove, IsVerbose = False):

    DidWork = True
    
    remote = user + '@' + server + ':'
    outfile = '.output_rsync_log'

    m = re.search('(.*)(star)(.*)', files)
    if m:
        files = m.group(1) + '*' + m.group(3)
    else:
        didWork = False
        return
    
    remote_command = user + "@" + server + " 'ls " + dir + "/" + files + "'"
    # check to see if the remote directory exists:

    if (IsVerbose):
        print('Getting list of remote files')
    command = 'ssh ' + remote_command + ' > .test_file'
    DidWork = system.run_command(command, verbose = IsVerbose)

    if (DidWork):
        with open('.test_file', 'r') as f:
            allFiles = f.readlines()
            for file in allFiles:
                sep = file.strip().split('/')
                localFile = sep[-1]
                remoteFile = remote + file.strip()
                rsync = 'rsync -av ' + remoteFile + ' .'
                if (not IsVerbose):
                    rsync = rsync + ' >> .rsync_out 2>&1'
                DidWork = system.run_command(rsync, verbose = IsVerbose)
                if (DoRemove):
                    remote_command = \
                        user + "@" + server + " 'rm -f " + file.strip() + "'"
                    command = 'ssh ' + remote_command
                    if (not IsVerbose):
                        command = command + ' >> .rm_out'
                    DidWork = system.run_command(command, verbose = IsVerbose)
                    time.sleep(5)

    return DidWork
