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
from pyitm.general import time_conversion

# This code will plot log files as a function of time
#    (assuming that there are columns for year, month, day, etc...)
# Simple plots can be made by calling this with a file and a variable
#   and it will make a plot of that variable.
#   call with -list to list the variables in the file.
#
# for example:
# logfile_plotter.py mylogfile.txt -vars 6 7 8 -plotfile=simple.png
#
# Can make dramatically more complicated plots by creating a CSV file
# that contains the file name, the variable, and the plot style of this line
# This allows incredible flexibility in the plots.
# You can use this by doing something like:
# creating a CSV file called plot.csv that contains:
# logfilename1,varToRead1,label,linestyle,linewidth,linecolor
# logfilename2,varToRead2,label,linestyle,linewidth,linecolor
#
# where the filenames can be unique or repeating.
#           vars can be unique or repeating.
# etc.
#
# then call this plotter with:
# logfile_plotter.py -csv=plot.csv none -vars=0 -plotfile=myplot.png
#
# you can also control the start time and end time for the x-axis
# as well as the min and max values for the y-axis


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

    parser.add_argument('-csv',
                        help = 'file to use to describe logfiles and styles',
                        default = 'none')

    parser.add_argument('-xlabel',
                        help = 'label for the x-axis',
                        default = 'none')

    parser.add_argument('-start',
                        help = 'start time of x-axis (yyyy-mm-dd_hh:mm)',
                        default = 'none')
    
    parser.add_argument('-title',
                        help = 'title of plot',
                        default = 'none')

    parser.add_argument('-end',
                        help = 'end time of x-axis (yyyy-mm-dd_hh:mm)',
                        default = 'none')
    
    parser.add_argument('-min', type = float,
                        help = 'Minimum value for y-axis',
                        default = 1e32)
    parser.add_argument('-max', type = float,
                        help = 'Maximum value for y-axis',
                        default = -1e32)

    parser.add_argument('-vars',
                        help = 'var(s) to plot (e.g. -vars 0 or -vars 0 1 2',
                        default = [7], nargs = '+', type = int)

    parser.add_argument('-list',  \
                        action='store_true', default = False, \
                        help = 'List variables in log file')

    args = parser.parse_args()

    return args


def read_logfile_styles(csvFile):

    # format is:
    # logfilename1,varToRead1,label,linestyle,linewidth,linecolor
    # logfilename2,varToRead2,label,linestyle,linewidth,linecolor
    # logfilename3,varToRead3,label,linestyle,linewidth,linecolor

    filenames = []
    vars = []
    labels = []
    widths = []
    styles = []
    colors = []

    print('Reading Controller CSV file :', csvFile)
    
    with open(csvFile, 'r') as f:
        lines = f.readlines()
        f.close()
        if (len(lines) > 5):
            i = 0
            while (i < len(lines)):
                l = lines[i].strip().split(',')
                print('  -> controls : ', l)
                filenames.append(l[0])
                vars.append([int(l[1])])
                labels.append(l[2])
                styles.append(l[3])
                widths.append(float(l[4]))
                colors.append(l[5])
                i = i + 1
                if (i < len(lines)):
                    if (len(lines[i].strip()) < 2):
                        i = len(lines)

    controller = {'files': filenames,
                  'vars': vars,
                  'labels': labels,
                  'widths': widths,
                  'styles': styles,
                  'colors': colors}
    return controller


# ----------------------------------------------------------------------
# Main Code
# ----------------------------------------------------------------------

args = parse_args()
plt.rcParams.update({'font.size': 14})

useController = False
if (args.csv != 'none'):
    controller = read_logfile_styles(args.csv)
    if (len(controller['files']) > 0):
        useController = True

if (useController):
    filenames = controller['files']
    vars = controller['vars']
    color = controller['colors'][0]
    linestyle = controller['styles'][0]
    linewidth = controller['widths'][0]
    labels = controller['labels']
    label = labels[0]
else:
    filenames = args.files
    vars0 = args.vars
    print(vars0)
    vars = []
    color = None
    linestyle = None
    linewidth = None
    if (len(filenames) == 1):
        labels = None
        label = None
    else:
        labels = filenames
        label = labels[0]
    
    for f in filenames:
        vars.append(vars0)

if (args.xlabel != 'none'):
    xlabel = args.xlabel
else:
    xlabel = None
    
plotfile = args.plotfile

logfilename = filenames[0]
print(' -> Reading logfile, with vars: ', logfilename, vars[0])
data_to_plot = logfile.get_logdata(logfilename, vars[0], just_list=args.list)

if (args.title == 'none'):
    title = ''
else:
    title = args.title

if (args.start == 'none'):
    xStart = data_to_plot['times'][0]
else:
    xStart = time_conversion.convert_string_to_datetime(args.start)
if (args.end == 'none'):
    xEnd = data_to_plot['times'][-1]
else:
    xEnd = time_conversion.convert_string_to_datetime(args.end)

if (args.max > -1e32):
    yMax = args.max
else:
    yMax = None
if (args.min < 1e32):
    yMin = args.min
else:
    yMin = None
    
if len(filenames) > 1:
    # create fig & ax from first logfile
    fig, ax = line_plots.lineplot_data(data_to_plot,
                                       outFile=None,
                                       xlabel = xlabel,
                                       label = label,
                                       title = title,
                                       color = color,
                                       linestyle = linestyle,
                                       linewidth = linewidth,
                                       yMin = yMin,
                                       yMax = yMax,
                                       xStart = xStart,
                                       xEnd = xEnd)
    for iFile in range(1, len(filenames)):
        # add as many lines as we have
        logfilename = filenames[iFile]
        print(' -> Reading logfile, with vars: ', logfilename, vars[iFile])
        data_to_plot = logfile.get_logdata(logfilename, vars[iFile])
        if (labels != None):
            label = labels[iFile]            
        if (useController):
            color = controller['colors'][iFile]
            linestyle = controller['styles'][iFile]
            linewidth = controller['widths'][iFile]
            print('  -> style : ', color, linestyle, linewidth)
        if iFile < len(filenames) - 1:
            fig, ax = line_plots.lineplot_data(data_to_plot,
                                               outFile=None, fig=fig, ax=ax,
                                               xlabel = xlabel,
                                               label = label,
                                               title = title,
                                               color = color,
                                               linestyle = linestyle,
                                               linewidth = linewidth,
                                               yMin = yMin,
                                               yMax = yMax,
                                               xStart = xStart,
                                               xEnd = xEnd)
        else:
            if (labels != None):
                nLabels = len(labels)
            else:
                nLabels = 0
            line_plots.lineplot_data(data_to_plot,
                                     outFile=plotfile,
                                     fig=fig,
                                     ax=ax,
                                     xlabel = xlabel,
                                     title = title,
                                     label = label,
                                     nLabels = nLabels,
                                     color = color,
                                     linestyle = linestyle,
                                     linewidth = linewidth,
                                     yMin = yMin,
                                     yMax = yMax,
                                     xStart = xStart,
                                     xEnd = xEnd)
else:
    # just plot single file
    print('plotting...')
    line_plots.lineplot_data(data_to_plot, outFile=args.plotfile,
                             xlabel = xlabel,
                             title = title,
                             label = label,
                             color = color,
                             linestyle = linestyle,
                             linewidth = linewidth,
                             yMin = yMin,
                             yMax = yMax,
                             xStart = xStart,
                             xEnd = xEnd)
