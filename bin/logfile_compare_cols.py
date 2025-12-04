#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from pylab import cm
import sys
import argparse
import datetime as dt

from pyitm.fileio import logfile
from pyitm.modeldata import utils
from pyitm.plotting import line_plots
from pyitm.general import statistics

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(description = 'Plot AMIE files')
    parser.add_argument('files', metavar = 'file', nargs = '+', \
                        help = 'Files to process')

    parser.add_argument('-vars',
                        help = 'var(s) to plot (e.g. -vars 0 or -vars 0 1 2',
                        default = [6,7], nargs = 2, type = int)

    parser.add_argument('-list',  \
                        action='store_true', default = False, \
                        help = 'List variables in log file')

    parser.add_argument('-v',  \
                        action='store_true', default = False, \
                        help = 'Output more information')

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# Main Code
# ----------------------------------------------------------------------

args = parse_args()
plt.rcParams.update({'font.size': 14})

logfilename = args.files[0]
logdata = logfile.read_logfile(logfilename=logfilename, verbose=args.v)

vars = []
for i, var in enumerate(logdata.keys()):
    if (args.list):
        print('%d. ' % i + var)
    vars.append(var)

if (args.list):
    exit()

# Only want to calculate the normalization once, so it doesn't change for each file:
v1 = vars[args.vars[0]]
vals1 = logdata[v1]
norm = statistics.calc_normalization(vals1, doReport = True)

print('Date & $M$ & $M_{mean}$ & $M_{std}$ & RMS & RMS$_{mean}$ & RMS$_{std}$ \\\\ \\hline')
    
for logfilename in args.files:

    logdata = logfile.read_logfile(logfilename=logfilename, verbose=args.v)

    times = logfile.calc_times(logdata)

    v1 = vars[args.vars[0]]
    v2 = vars[args.vars[1]]
    vals1 = logdata[v1]
    vals2 = logdata[v2]

    vals1 = vals1/norm
    vals2 = vals2/norm

    m = statistics.calc_mean_diff(vals1, vals2, doReport = args.v)
    mnm = statistics.calc_normMean_mean_diff(vals1, vals2, doReport = args.v) * 100.0
    mns = statistics.calc_normStd_mean_diff(vals1, vals2, doReport = args.v) * 100.0

    r = statistics.calc_rms_diff(vals1, vals2, doReport = args.v)
    rnm = statistics.calc_normMean_rms_diff(vals1, vals2, doReport = args.v) * 100.0
    rns = statistics.calc_normStd_rms_diff(vals1, vals2, doReport = args.v) * 100.0

    s = times[0].strftime('%Y-%m-%d')
    s = s + ' & %0.3f '% m
    s = s + '& %0.1f\\%% '% mnm
    s = s + '& %0.1f\\%% '% mns
    s = s + '& %0.3f '% r
    s = s + '& %0.1f\\%% '% rnm
    s = s + '& %0.1f\\%% \\\\'% rns

    print(s)
