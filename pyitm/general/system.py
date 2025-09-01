#!/usr/bin/env python3

import os

# ----------------------------------------------------------------------
# do system command
# ----------------------------------------------------------------------

def run_command(command, verbose = False):
    if (verbose):
        print("   -> Running Command : ")
        print("      ", command)
    os.system(command)
    return True

