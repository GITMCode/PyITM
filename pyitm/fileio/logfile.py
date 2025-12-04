#!/usr/bin/env python

from datetime import datetime
import numpy as np
import os, glob

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


def calc_times(logdata):

    times = []
    iYear = -1
    iMonth = -1
    iDay = -1
    iHour = -1
    iMinute = -1
    iSecond = -1
    vars = list(logdata.keys())
    if ('Year' in vars):
        iYear = vars.index('Year')
    if ('Month' in vars):
        iMonth = vars.index('Month')
    if ('Day' in vars):
        iDay = vars.index('Day')
    if ('Hour' in vars):
        iHour = vars.index('Hour')
    if ('Minute' in vars):
        iMinute = vars.index('Minute')
    if ('Second' in vars):
        iSecond = vars.index('Second')
    if ((iYear >= 0) and
        (iMonth >= 0) and
        (iDay >= 0) and
        (iHour >= 0) and
        (iMinute >= 0) and
        (iSecond >= 0)):
        for i, year in enumerate(logdata['Year']):
            t = datetime(int(year),
                         int(logdata['Month'][i]),
                         int(logdata['Day'][i]),
                         int(logdata['Hour'][i]),
                         int(logdata['Minute'][i]),
                         int(logdata['Second'][i]))
            times.append(t)
    return times

#-----------------------------------------------------------------------------

def read_logfile(logfilename=None, datadir=None, verbose=False):
    """ Read GITM logfile to a dictionary. Either datapath or logfile must be provided.
    
    Parameters
    ----------
    logfile: (str/path) - path to logfile to read. Default None.
    datadir: (str/path) - path to GITM data directory. Default None.
    verbose: (bool) - print extra info when running? Default = False
    
    Returns
    -------
    logdata: (dict) - dictionary containing logfile information

    Notes
    -----
    - If datadir is provided, will look for the most recent logfile in that directory.
    - This should be the directory in which model outputs are stored.
    """

    if datadir: # Grab the most highest numbered logfile
        #TODO: Read & appead all logfiles
        if logfilename:
            print("Both logfile and datadir provided. Ignoring logfile and using datadir.")
        logfilename = sorted(glob.glob(os.path.join(datadir, 'log*.dat')))[-1]

    if not os.path.exists(logfilename):
        raise FileNotFoundError(f"Could not find logfile: {logfilename}")
    elif verbose:
        print(f"Reading logfile: {logfilename}")

    logdata = {}

    saving = False
    readVars = True
    with open(logfilename, 'r') as f:
        for n, line in enumerate(f.readlines()):
            if (readVars):
                v = line.strip()
                if (len(v) < 1):
                    readVars = False
                else:
                    logdata[v] = []
                    vars.append(v)

            if line.startswith('#VARIABLES'):
                readVars = True
                vars = []
            if saving:
                # first line (after #START) has the headers
                if len(logdata.keys()) == 0:
                    vars = []
                    for col in line.strip().split():
                        while (col in vars):
                            col = col+'0'
                        vars.append(col)
                        logdata[col] = []
                else:
                    # Read in the data
                    for col, val in zip(logdata.keys(), line.strip().split()):
                        logdata[col].append(val)

            if line.startswith('#START'):
                saving=True
                if verbose:
                    print(f"  -> Found start of data at line {n}")
                    print(f"  -> Reading columns : {list(logdata.keys())}")
    if verbose:
        print(f"  -> Read {len(logdata[list(logdata.keys())[0]])} lines of data")

    for col in logdata.keys():
        try:
            logdata[col] = np.array(logdata[col], dtype=float)
        except:
            if verbose:
                print(f"  -> Could not convert column '{col}' to float. Leaving as string.")
            pass

    return logdata


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

