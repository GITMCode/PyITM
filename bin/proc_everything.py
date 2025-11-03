#!/usr/bin/python3

import os, glob, datetime, subprocess, argparse
import numpy as np
# import time


def debug_mes(message, level):
    if level==0:
        print('\n' + '-'*len(message))
        print(message)
        print('-'*len(message))
    else:
        print(" " + '-'*level + '> ' + str(message))
    return

def dt2str(dt):
    """ time delta to string
    """
    dt = dt.total_seconds()
    hrs = dt//3600
    mins = (dt - hrs*3600)//60
    sec = dt - (hrs * 3600 + mins  * 60)
    return ':'.join([str(int(t)).rjust(2, '0') for t in [hrs, mins, sec]])

def shell_cmd(cmd, dict_key=None, stderr_dict=None, stdout_dict=None):
    """
    run a shell command. Return (appended) stderr and stdout dicts
    """
    
    debug_mes(f"$ {cmd}.", 1)

    out = subprocess.run(cmd, shell=True, capture_output=True)
    if out.returncode != 0:
        debug_mes(f"Error running {cmd}.", 3)
        debug_mes(f"See error dictionary key '{dict_key}'", 3)
        stderr_dict[dict_key][cmd] = out.stderr.decode()

    stdout_dict[dict_key][cmd] = out.stdout.decode()

    return stderr_dict, stdout_dict


def parse_settings(setting_filename):

    with open(setting_filename, 'r') as f:
        all_text = f.readlines()
    
    settings = {}
    mostrecentkey = None
    for nline, text in enumerate(all_text):
        if text.startswith('#'):
            continue
        if text.startswith('-'):
            mostrecentkey = text[1:].strip()
            settings[mostrecentkey] = []
            continue
        if text.strip() == '': # empty! Catches lines with just spaces too 
            continue
        # continuation line
        if text.startswith(' ') and mostrecentkey is not None:
            settings[mostrecentkey][-1] = ' '.join([settings[mostrecentkey][-1], text.strip()])
            continue
        # Otherwise, add it to settings list
        settings[mostrecentkey].append(text.strip())

    return settings

def extract_setting(value, wanttype):
    """
    The settings could be a bool, string, list, etc. but are always read into a list of str's

    Here we take that list of str and convert it to the type it needs to be

    For wanttype, use the python data type
     -> Some things can be bool or str, for those use the string 'boolstr' for wanttype
    
    """

    if wanttype is str:
        if len(value) == 1:
            return value[0]
        else:
            raise TypeError(f"Cannot convert {value} to string")
    elif wanttype is list:
        return value
    elif wanttype is bool:
        # Convert to str
        if type(value) is bool:
            return value
        else:
            value = extract_setting(value, str)
        if value.title() == 'True':
            return True
        elif value.title() == 'False':
            return False
        else:
            raise TypeError(f"Cannot convert {value} to boolean")
    elif wanttype == 'boolstr':
        try:
            return extract_setting(value, bool)
        except TypeError:
            return extract_setting(value, str)
    else:
        raise ValueError(f"Cannot convert to {wanttype} yet. Perhaps you want to write a converter?")



## Get the settings!!

parser = argparse.ArgumentParser()
parser.add_argument("settings_file", help='Path to the settings file to read. Required!')
fname = parser.parse_args().settings_file
print('Reading settings from:', fname)
settings = parse_settings(fname)

# Here's where we actually set the settings, haha
# Use pop to check for unused settings at the end

fmt_str = extract_setting(settings.pop('fmt_str'), str)
dates_to_run = extract_setting(settings.pop('dates_to_run'), list)

each_dir_commands = extract_setting(settings.pop('each_dir_commands'), list)

# These are optional settings. If its present we pop it.
# Otherwise give the default value (False usually)

newname = extract_setting(settings.pop('newname', False), 'boolstr')
oldname = extract_setting(settings.pop('oldname', False), 'boolstr')

move_existing = extract_setting(settings.pop('move_existing', False), 'boolstr')

makedayone = extract_setting(settings.pop('makedayone', False), bool)

dumperr2file = extract_setting(settings.pop('dumperr2file', True), bool)
dumpout2file = extract_setting(settings.pop('dumpout2file', True), bool)



## vars we *need*
ndates = len(dates_to_run)
pwd = os.getcwd()

### The main code!

errs = {}
outs = {}
times = []
t0 = datetime.datetime.now()

for ndir, thisdate in enumerate(dates_to_run):

    tstart = datetime.datetime.now()
    errs[thisdate] = {}
    outs[thisdate] = {}

    thispath = fmt_str.replace('[DATE]', thisdate)
    debug_mes(f"({ndir+1}/{ndates})\tBeginning {thispath} ", 0)

    if not os.path.exists(thispath):
        debug_mes(f"'{thispath}' DOES NOT EXIST!!", 2)
        if newname is not None and oldname is not None:
            if os.path.exists(thispath.replace(newname, oldname)):
                p = thispath[:thispath.rfind('/')]
                cmd = f"mv {p.replace(newname, oldname)}  {p} "
                mv_cmd = subprocess.run(cmd, shell=True, capture_output=True)
                if mv_cmd.returncode == 0:
                    debug_mes(f"{p} moved successfully!!", 1)
                else:
                    debug_mes(f"{p} could not be moved. Check error logs for {thisdate}", 1)
                    errs[thisdate] = mv_cmd.stderr.decode()
                    continue
        else:
            errs[thisdate] = ['Did not exist. No renaming options were specified.']
            continue

    os.chdir(thispath)

    if move_existing:
        # rename run folders (from move_existing)
        f2move = glob.glob('*.txt')
        f2move.extend(glob.glob('*.png'))
        f2move.extend(glob.glob('*.csv'))
        if len(f2move) == 0:
            os.mkdir(move_existing)
            debug_mes(f"Moving {len(f2move)} existing files.", 2)
            for file in f2move:
                errs, outs = shell_cmd(f"mv {file} {os.path.join(move_existing, file)}",
                                    thisdate, errs, outs)


    if (not os.path.exists('dayone')) and makedayone:
        # move outputs from the first day
        # gets the sorted list of files, finds the first date (not time)
        date_str = sorted(glob.glob('3DALL*'))[1].split('_')[1]
        # Glob all files with this date
        fs2move = f'*{date_str}*'
        os.mkdir('dayone')
        debug_mes(f"Moving '{fs2move}' files to 'dayone/' folder", 2)
        errs, outs = shell_cmd(f"mv {fs2move} dayone/", thisdate, errs, outs)


    for plot_script in each_dir_commands:
        # debug_mes(f" $ {plot_script} ", 4)
        errs, outs = shell_cmd(plot_script, thisdate, errs, outs)



    # Timing information & wrapup
    this_dt = datetime.datetime.now() - tstart
    tot_dt = datetime.datetime.now() - t0
    times.append(this_dt)

    debug_mes(f"done with '{thisdate}'. This iteration took: {dt2str(this_dt)}.", 1)
    debug_mes(f'total={dt2str(tot_dt)}. Est remaining: {dt2str(np.mean(times) * (len(dates_to_run) - ndir-1))}',
              3)

    os.chdir(pwd)


# Error information to screen

debug_mes("Error information:", 0)
noerrs = True
for k in errs.keys():
    for cmd in errs[k]:
        print('\t', cmd, k)
        print(errs[k][cmd])
        noerrs = False
if noerrs:
    print(" None :)")

# Error information to file
if dumperr2file:
    with open('stderr_pyitm.txt', 'w') as f:
        for k in errs.keys():
            f.write('='*50 + '\n' + k + '\n')
            for cmd in errs[k]:
                # f.write('-'*50 +'\n'+k+'\n')
                f.write(cmd)
                f.write(errs[k][cmd])
                f.write('\n\n\n')

if dumpout2file:
    with open('stdout_pyitm.txt', 'w') as f:
        for k in outs.keys():
            f.write(k+'\n')
            for cmd in outs[k]:
                f.write(cmd + '\n')
                f.write(outs[k][cmd]+'\n')
