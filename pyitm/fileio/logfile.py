#!/usr/bin/env python

from datetime import datetime
import numpy as np
import os

#-----------------------------------------------------------------------------
# Write a single ascii line out to an already opened file
#-----------------------------------------------------------------------------

def write_line(fp, sLine):
    sLineN = sLine + '\n'
    fp.write(sLineN.encode())
    return

#-----------------------------------------------------------------------------
# Here is a useful routine for writing out data to a log file
#-----------------------------------------------------------------------------

def write_log(data, fileHeader = 'log', message = ''):

    sTime = data['times'][0].strftime('%Y%m%d')
    fileout = fileHeader + '_' + sTime + '.txt'
    print('--> Writing file : ' + fileout)
    fp = open(fileout, 'wb')

    write_line(fp, '')
    write_line(fp, message)
    
    pwd = os.getcwd()
    write_line(fp, '')
    write_line(fp, '#DIRECTORY')
    write_line(fp, pwd)

    if ('alt' in data):
        if (np.isscalar(data['alt'])):
            alt = data['alt']
        else:
            alt = mean(np.array(data['alt']))
        sLine = '%f' % alt
        write_line(fp, '')
        write_line(fp, '#ALTITUDE')
        write_line(fp, sLine)
    
    write_line(fp, '')
    write_line(fp, '#VARIABLES')
    write_line(fp, 'Year')
    write_line(fp, 'Month')
    write_line(fp, 'Day')
    write_line(fp, 'Hour')
    write_line(fp, 'Minute')
    write_line(fp, 'Second')
    for key in data.keys():
        if ((key != 'times') and (key != 'alt')):
            write_line(fp, key)

    write_line(fp, '')
    write_line(fp, '#START')

    for i, t in enumerate(data['times']):

        sLine = t.strftime(' %Y %m %d %H %M %S')

        for key in data.keys():
            if ((key != 'times') and (key != 'alt')):
                sLine = sLine + ' %e' % data[key][i]
        write_line(fp, sLine)
        
    fp.close()
    
    return


