#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from pylab import cm
import sys
import argparse
import datetime as dt

from pyitm.fileio import gitmio
from pyitm.modeldata import utils
from pyitm.plotting import line_plots

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(description = 'Plot AMIE files')
    parser.add_argument('files', metavar = 'file', nargs = '+', \
                        help = 'Files to process')

    parser.add_argument('-plotfile',
                        help = 'output file for plot',
                        default = 'logfile_timeline.png')    

    parser.add_argument('-vars',
                        help = 'var(s) to plot (e.g. -vars 0 or -vars 0 1 2',
                        default = [7], nargs = '+', type = int)

    parser.add_argument('-list',  \
                        action='store_true', default = False, \
                        help = 'List variables in log file')

    args = parser.parse_args()

    return args


def get_logdata(logfile, vars2plot, just_list=False):
    logData = gitmio.read_logfile(logfile=logfile)
    vars = []
    for i, k in enumerate(logData.keys()):
        vars.append(k)
        if (just_list):
            print('%02d. ' % i, vars[-1])

    if (just_list):
        exit()

    iY_ = vars[1]
    iM_ = vars[2]
    iD_ = vars[3]
    iH_ = vars[4]
    iMi_ = vars[5]
    iS_ = vars[6]

    nTimes = len(logData[vars[0]])
    times = []

    for iT in range(nTimes):
        year = int(logData[iY_][iT])
        month = int(logData[iM_][iT])
        day = int(logData[iD_][iT])
        hour = int(logData[iH_][iT])
        minute = int(logData[iMi_][iT])
        second = int(logData[iS_][iT])
        t = dt.datetime(year, month, day, hour, minute, second)
        times.append(t)

    data_to_plot = {'times': times,
                    'note': logfile,
                    'vars': []}
        
    for var in vars2plot:
        print(var, vars[var])
        data_to_plot['vars'].append(vars[var])
        data_to_plot[vars[var]] = logData[vars[var]]
    
    return data_to_plot

# ----------------------------------------------------------------------
# Main Code
# ----------------------------------------------------------------------

args = parse_args()
plt.rcParams.update({'font.size': 14})

plotfile = args.plotfile

logfile = args.files[0]

data_to_plot = get_logdata(logfile, args.vars, just_list=args.list)

if len(args.files) > 1:
    print('multiple')
    # create fig & ax from first logfile
    fig, ax = line_plots.lineplot_data(data_to_plot, outFile=None)
    for iFile in range(1, len(args.files)):
        # add as many lines as we have
        logfile = args.files[iFile]
        data_to_plot = get_logdata(logfile, args.vars)
        if iFile < len(args.files) - 1:
            fig, ax = line_plots.lineplot_data(data_to_plot,
                                               outFile=None, fig=fig, ax=ax)
        else:
            line_plots.lineplot_data(data_to_plot,
                                               outFile=plotfile, fig=fig, ax=ax)

else:
    # just plot single file
    line_plots.lineplot_data(data_to_plot, outFile=None)
