#!/usr/bin/env python3

import numpy as np
import os


# ----------------------------------------------------------------------------
# CSV-based variable mapping loader
# ----------------------------------------------------------------------------

def _load_variable_csv(csv_path=None):
    """Load variables.csv and return structured variable data.

    Parses everything after the ``#START`` marker. Each row becomes a dict
    with keys: shortname, longname, unit, prettyname, modelnames (list).
    Semicolons in model-name columns are replaced with commas so that
    model names like ``Vn(up,O(3P))`` can be stored in a comma-delimited file.

    Args:
        csv_path: Path to the CSV file. Defaults to variables.csv in the
            same directory as this module.

    Returns:
        rows: list of row dicts (one per variable).
        alias_index: dict mapping any known name (lowercase) to its row.
            Keys include shortname, longname, and all model-name aliases.
    """
    if csv_path is None:
        csv_path = os.path.join(os.path.dirname(__file__), 'variables.csv')

    rows = []
    alias_index = {}

    with open(csv_path, 'r') as f:
        lines = f.readlines()

    started = False
    for line in lines:
        line = line.strip()
        if line == '#START':
            started = True
            continue
        if not started or not line or line.startswith('#'):
            continue

        parts = [p.strip() for p in line.split(',')]
        if len(parts) < 2:
            continue

        row = {
            'shortname': parts[0],
            'longname': parts[1] if len(parts) > 1 else parts[0],
            'unit': parts[2] if len(parts) > 2 else '',
            'prettyname': parts[3] if len(parts) > 3 else '',
            'modelnames': []
        }

        # Columns 4+ are model name aliases
        for i in range(4, len(parts)):
            alias = parts[i].strip()
            if alias:
                # Semicolons stand for commas in model names
                alias_real = alias.replace(';', ',')
                row['modelnames'].append(alias_real)

        rows.append(row)

        # Build reverse index: shortname, longname, and all aliases -> row
        alias_index[row['shortname'].lower()] = row
        alias_index[row['longname'].lower()] = row
        for alias in row['modelnames']:
            alias_index[alias.lower()] = row

    return rows, alias_index


# The general idea here is that each code has a bunch of variables that are
# named in different ways.  What we essentially want for each code is that
# you can call a plotter or reader with variable names in different ways:
# - as an actual number (e.g., 15)
# - as a number string  (e.g., '15')
# - as a short variable (e.g., 'Tn')
# - as the actual variable (e.g., 'Temperature_neutral')
# The codes should then be able to figure out what you are asking for.
# These codes then get a bit complicated, since they have to do all sorts
# of interpreting of what you want versus what you provided.

# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------

def convert_filename(filename, convertFile = 'name_convert.csv'):
    """Look up a filename in a CSV mapping file and return the converted name.

    Args:
        filename: The filename to convert.
        convertFile: Path to a two-column CSV mapping old names to new names.

    Returns:
        The converted filename, or the original if no mapping is found.
    """
    filenames = []
    strings = []
    if (os.path.exists(convertFile)):
        with open(convertFile, 'r') as f:
            lines = f.readlines()
            for iLine, line in enumerate(lines):
                fs = line.split(',')
                filenames.append(fs[0])
                strings.append(fs[1])
        if (len(filenames) > 0):
            iFile = find_string(filename, filenames)
            return strings[iFile]
        else:
            return filename
    
# ----------------------------------------------------------------------------
# This code takes something like 'Temperature (K)' and converts it to
# 'Temperature'. This is extremely useful when naming files with variable names
# ----------------------------------------------------------------------------

def strip_varname(varnameIn):
    """Strip a parenthesized unit suffix from a variable name.

    'Temperature (K)' -> 'Temperature '
    'Vn(east)'        -> 'Vn(east)'  (unchanged, no space before paren)
    """
    ind = varnameIn.find('(')
    if (ind > 0):
        varnameOut = varnameIn[0:ind]
    else:
        varnameOut = varnameIn
        
    return varnameOut

# ----------------------------------------------------------------------------
# this returns the index of the array that matches the string, if it is found
# ----------------------------------------------------------------------------

def find_string(item, stringList):
    """Return the index of *item* in *stringList*, or -1 if not found."""
    iVal = -1
    if (item in stringList):
        i = 0
        while (i < len(stringList)):
            if (stringList[i] == item):
                iVal = i
                i = len(stringList)
            i += 1
    return iVal

# ----------------------------------------------------------------------------
# Take variable numbers (as a number or a string of a number) or names and 
# make sure that they are all names
# coming out.  If the user enters a name, it will just use the name.  If the
# user enters a number, this function will take the Nth variable in the
# list of variables returned from the header.
# ----------------------------------------------------------------------------

def convert_number_to_var(varList, header = None):
    """Convert numeric variable indices to variable names using a file header.

    Accepts a single value or list of values. Numbers (int or numeric string)
    are replaced with the corresponding variable name from ``header['vars']``.
    Non-numeric strings are passed through unchanged.
    """
    if (np.isscalar(varList)):
        if (not isinstance(varList, str)):
            varList = '%d' % int(varList)
        if (varList.isnumeric()):
            if (header):
                sVars = [header['vars'][int(varList)]]
        else:
            sVars = [varList]

    else:
        sVars = []
        for var in varList:
            if (not isinstance(var, str)):
                var = '%d' % int(var)
            if (var.isnumeric()):
                if (header):
                    sVars.append(header['vars'][int(var)])
            else:
                sVars.append(var)
            
    return sVars
    
# ----------------------------------------------------------------------------
# This function takes a variable name and tries to figure out what 
# number it is in the file.  In order for this to work, the header
# has to be provided, since it has to look for the variable in the header
# ----------------------------------------------------------------------------

def convert_var_to_number(varList, header = None):
    """Convert variable names to their numeric indices in a file header.

    Accepts a single value or list. Tries matching against header['shortname'],
    header['vars'], and header['longname'] in that order.
    """
    if (np.isscalar(varList)):
        if (varList.isnumeric()):
            iVars = [int(varList)]
        else:
            if (header == None):
                print('Non number variables are not supported yet!')
                iVars = [3]
            else:
                sVar = match_var_name([varList], header)[0]
                print('sVar -> ', sVar)
                iV = find_string(sVar, header['shortname'])
                if (iV < 0):
                    iV = find_string(sVar, header['vars'])
                if (iV < 0):
                    iV = find_string(sVar, header['longname'])
                iVars = [iV]
    else:
        iVars = []
        for var in varList:
            if (var.isnumeric()):
                iVars.append(int(var))
            else:
                if (header == None):
                    print('Non number variables are not supported yet!')
                    iVars = [3]
                else:
                    sVar = match_var_name([var], header)[0]
                    iV = find_string(sVar, header['shortname'])
                    if (iV < 0):
                        iV = find_string(sVar, header['vars'])
                    if (iV < 0):
                        iV = find_string(sVar, header['longname'])
                    iVars.append(iV)

    return iVars

#-----------------------------------------------------------------------------
# take a list of variables, and try to figure out what the user is 
# actually asking for.  First, everything is converted to lower
# case so it can match 'Temperature' with 'temperature'. Then,
# it compares to variables in the header, longnames, and shortnames.
#-----------------------------------------------------------------------------

def match_var_name(varsIn, header):
    """Match user-provided variable names against a file header (case-insensitive).

    Checks header['vars'], header['longname'], and header['shortname']
    in that order. Raises KeyError if a variable cannot be matched.
    """
    varsOut = []

    for varIn in varsIn:
        isFound = False
        for var in header['vars']:
            if (var.lower() == varIn.lower()):
                varsOut.append(var)
                isFound = True
        if (not isFound):
            for iVar, var in enumerate(header['longname']):
                if (var.lower() == varIn.lower()):
                    varsOut.append(header['vars'][iVar])
                    isFound = True
        if (not isFound):
            for iVar, var in enumerate(header['shortname']):
                if (var.lower() == varIn.lower()):
                    varsOut.append(header['vars'][iVar])
                    isFound = True
        if (not isFound):
            varsOut.append('NotFound')
            print('Could not find variable : ', varIn)
            print('  -> Should be able to list variables by putting -list or running with -verbose')
            raise KeyError

    return varsOut

# ----------------------------------------------------------------------------
# CSV-backed variable name mapping functions
# All mappings are now defined in variables.csv
# ----------------------------------------------------------------------------

_ROWS = None
_INDEX = None

def _ensure_loaded():
    global _ROWS, _INDEX
    if _ROWS is None:
        _ROWS, _INDEX = _load_variable_csv()


def remap_variable_names(varsIn):
    """Map model output names -> 'shortname (unit)'.

    Accepts a single string or a list of strings.
    Always returns a list.
    """
    _ensure_loaded()
    if np.isscalar(varsIn):
        varsIn = [varsIn]
    varsOut = []
    for var in varsIn:
        row = _INDEX.get(var.lower())
        if row and row['unit']:
            varsOut.append(f"{row['shortname']} ({row['unit']})")
        elif row:
            varsOut.append(row['shortname'])
        else:
            varsOut.append(var)
    return varsOut


def get_short_name(name):
    """Look up the shortname for any known variable name.

    Accepts shortnames, longnames, or model-output names (case-insensitive).
    Returns the canonical shortname (e.g. 'Tn', 'O+', 'eFlux').
    Falls back to stripping a parenthesized suffix if the name is unknown.
    """
    _ensure_loaded()
    row = _INDEX.get(name.lower())
    return row['shortname'] if row else strip_varname(name)

def get_short_names(varsIn):
    """Look up shortnames for one or more variable names.

    Accepts a single string or a list of strings. Always returns a list.
    """
    if np.isscalar(varsIn):
        varsIn = [varsIn]
    return [get_short_name(v) for v in varsIn]

def get_long_name(name):
    """Look up a descriptive long name with units for a variable.

    Returns e.g. 'Neutral Temperature (K)', 'O+ Density (/m3)'.
    Falls back to the input string if the name is unknown.
    """
    _ensure_loaded()
    row = _INDEX.get(name.lower())
    if row:
        if row['unit']:
            return f"{row['longname']} ({row['unit']})"
        return row['longname']
    return name

def get_long_names(varsIn):
    """Look up long names for one or more variable names.

    Accepts a single string or a list of strings. Always returns a list.
    """
    if np.isscalar(varsIn):
        varsIn = [varsIn]
    return [get_long_name(v) for v in varsIn]

def get_unit(name):
    """Return the unit string for a variable (e.g. 'K', '/m3', 'mW/m2').

    Returns an empty string if the variable is unknown or has no unit.
    """
    _ensure_loaded()
    row = _INDEX.get(name.lower())
    return row['unit'] if row else ''

def get_pretty_name(name):
    """Return matplotlib-ready LaTeX name wrapped in $...$.

    Falls back to shortname (without $ wrapping) if no prettyname is defined.
    """
    _ensure_loaded()
    row = _INDEX.get(name.lower())
    if row and row['prettyname']:
        return f"${row['prettyname']}$"
    return get_short_name(name)

def get_label(name):
    """Return a matplotlib-ready axis label combining prettyname and unit.

    Examples:
        get_label('Tn')        -> '$T_n$ (K)'
        get_label('O+')        -> '$[O^+]$ (/m3)'
        get_label('PedCond')   -> '$\\Sigma_P$ (S)'
        get_label('UnknownVar') -> 'UnknownVar'

    If a prettyname is defined, uses '$prettyname$ (unit)'.
    Otherwise falls back to 'longname (unit)' via get_long_name.
    """
    _ensure_loaded()
    row = _INDEX.get(name.lower())
    if row and row['prettyname']:
        if row['unit']:
            return f"${row['prettyname']}$ ({row['unit']})"
        return f"${row['prettyname']}$"
    return get_long_name(name)

def list_variables():
    """Return a list of all known variable shortnames.

    Useful for discovering what variables are available in the CSV.
    """
    _ensure_loaded()
    return [row['shortname'] for row in _ROWS]

