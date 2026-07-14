#!/usr/bin/env python3

import numpy as np

# ----------------------------------------------------------------------------
# this function will read in a controller CSV file that contains
# information about what to plot.
# format is:
# filename1,varToRead1,label,linestyle,linewidth,linecolor
# filename2,varToRead2,label,linestyle,linewidth,linecolor
# filename3,varToRead3,label,linestyle,linewidth,linecolor
# ----------------------------------------------------------------------------

def read_logfile_styles(csvFile):

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
        if (len(lines) > 0):
            i = 0
            while (i < len(lines)):
                l = lines[i].strip().split(',')
                filenames.append(l[0])
                if (l[1].isnumeric()):
                    vars.append([int(l[1])])
                else:
                    vars.append([l[1]])
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

