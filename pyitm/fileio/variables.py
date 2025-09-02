#!/usr/bin/env python3

import numpy as np

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
# This code takes something like 'Temperature (K)' and converts it to
# 'Temperature'. This is extremely useful when naming files with variable names
# ----------------------------------------------------------------------------

def strip_varname(varnameIn):

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
            print('  -> Should be able to list variables by putting -list')

    return varsOut

#-----------------------------------------------------------------------------
# remap gitm variables - turn the output of GITM into human-digestible
# names AND add units to them.
#-----------------------------------------------------------------------------

def remap_variable_names(varsIn):

    mapVars = {
        'Rho' : 'rho (km/m3)',
        '[O(3P)]': '[O] (/m3)',
        '[O2]': '[O2] (/m3)',
        '[N2]': '[N2] (/m3)',
        '[N(4S)]': '[N] (/m3)',
        '[NO]': '[NO] (/m3)',
        '[He]': '[He] (/m3)',
        '[N(2D)]': '[N_2D] (/m3)',
        '[N(2P)]': '[N_2P] (/m3)',
        '[H]': '[H] (/m3)',
        '[CO2]': '[CO2] (/m3)',
        '[O(1D)]': '[O_1D] (/m3)',
        'Temperature': 'Tn (K)',
        'Vn(east)': 'Ve (m/s)',
        'Vn(north)': 'Vn (m/s)',
        'Vn(up)': 'Vv (m/s)',
        'Vn(up,O(3P))': 'Vv_O (m/s)',
        'Vn(up,O2)': 'Vv_O2 (m/s)',
        'Vn(up,N2)': 'Vv_N2 (m/s)',
        'Vn(up,N(4S))': 'Vv_N (m/s)',
        'Vn(up,NO)': 'Vv_NO (m/s)',
        'Vn(up,He)': 'Vv_He (m/s)',
        '[O_4SP_+]': '[O+] (/m3)',
        '[NO+]': '[NO+] (/m3)',
        '[O2+]': '[O2+] (/m3)',
        '[N2+]': '[N2+] (/m3)',
        '[N+]': '[N+] (/m3)',
        '[O(2D)+]': '[O_2D+] (/m3)',
        '[O(2P)+]': '[O_2P+] (/m3)',
        '[H+]': '[H+] (/m3)',
        '[He+]': '[He+] (/m3)',
        '[e-]': '[e-] (/m3)',
        'eTemperature': 'Te (K)',
        'iTemperature': 'Ti (K)',
        'Vi(east)': 'Vie (m/s)',
        'Vi(north)': 'Vin (m/s)',
        'Vi(up)': 'Viv (m/s)'}

    varsOut = []

    for var in varsIn:
        if (var in mapVars):
            varsOut.append(mapVars[var])
        else:
            varsOut.append(var)
    return varsOut

#-----------------------------------------------------------------------------
# remap gitm variables - take the GITM names and match them with 
# names that are very short and terse - that could be used for 
# filenames, for example.
#-----------------------------------------------------------------------------

def get_short_names(varsIn):

    if (np.isscalar(varsIn)):
        varsIn = [varsIn]
    
    mapVars = {
        'Rho' : 'rho',
        '[O(3P)]': 'O',
        '[O2]': 'O2',
        '[N2]': 'N2',
        '[N(4S)]': 'N',
        '[NO]': 'NO',
        '[He]': 'He',
        '[N(2D)]': 'N_2D',
        '[N(2P)]': 'N_2P',
        '[H]': 'H',
        '[CO2]': 'CO2',
        '[O(1D)]': 'O_1D',
        'density_O': 'O',
        'density_O_3P': 'O',
        'density_O2': 'O2',
        'density_N2': 'N2',
        'density_N': 'N',
        'density_N_4S': 'N',
        'density_NO': 'NO',
        'density_He': 'He',
        'density_N_2D': 'N_2D',
        'density_N_2P': 'N_2P',
        'density_H]': 'H',
        'density_CO2': 'CO2',
        'density_O_1D': 'O_1D',
        'Temperature': 'Tn',
        'temperature_neutral': 'Tn',
        'Vn(east)': 'Ve',
        'Vn(north)': 'Vn',
        'Vn(up)': 'Vv',
        'velocity_east_neutral': 'Ve',
        'velocity_north_neutral': 'Vn',
        'velocity_up_neutral': 'Vv',
        'Vn(up,O(3P))': 'Vv_O',
        'Vn(up,O2)': 'Vv_O2',
        'Vn(up,N2)': 'Vv_N2',
        'Vn(up,N(4S))': 'Vv_N',
        'Vn(up,NO)': 'Vv_NO',
        'Vn(up,He)': 'Vv_He',
        '[O_4SP_+]': 'O+',
        '[NO+]': 'NO+',
        '[O2+]': 'O2+',
        '[N2+]': 'N2+',
        '[N+]': 'N+',
        '[O(2D)+]': 'O_2D+',
        '[O(2P)+]': 'O_2P+',
        '[H+]': 'H+',
        '[He+]': 'He+',
        '[e-]': 'e-',
        'density_NO+': 'NO+',
        'density_O+': 'O+',
        'density_O2+': 'O2+',
        'density_N2+': 'N2+',
        'density_N+': 'N+',
        'density_O+_2D': 'O_2D+',
        'density_O+_2P': 'O_2P+',
        'density_H+': 'H+',
        'density_He+': 'He+',
        'density_e-': 'e-',
        'eTemperature': 'Te',
        'iTemperature': 'Ti',
        'temperature_ion': 'Ti',
        'temperature_electron': 'Te',
        'velocity_east_ion': 'Vie',
        'velocity_north_ion': 'Vin',
        'velocity_up_ion': 'Viv',
        'Vi(east)': 'Vie',
        'Vi(north)': 'Vin',
        'Vi(up)': 'Viv',
        'VerticalTEC': 'TEC',
        'Potential': 'pot',
        'PedersenConductance': 'PedCond',
        'HallConductance': 'HalCond',
        'Electron_Average_Energy_Diffuse': 'AveE',
        'Electron_Energy_Flux_Diffuse': 'eFlux',
        'Electron_Average_Energy_Wave': 'AveE_W',
        'Electron_Energy_Flux_Wave': 'eFlux_w',
        'Electron_Average_Energy_Mono': 'AveE_M',
        'Electron_Energy_Flux_Mono': 'eFlux_M',
        'Ion_Average_Energy': 'AveE_I',
        'Ion_Energy_Flux': 'eFlux_I',
        'AltIntJouleHeating(W/m2)': 'JouleHeat',
        'AltIntHeatingTransfer(W/m2)': 'HeatTrans',
        'AltIntEuvHeating(W/m2)': 'EuvHeat',
        'AltIntPhotoElectronHeating(W/m2)': 'PhotoElecHeat',
        'AltIntChamicalHeating(W/m2)': 'ChemHeat',
        'AltIntRadCooling(W/m2)': 'RadCool',
        'AltIntCO2Cooling(W/m2)': 'CO2Cool',
        'AltIntNOCooling(W/m2)': 'NOCool',
        'AltIntOCooling(W/m2)': 'OCool'}


    varsOut = []
    for var in varsIn:
        if (var in mapVars):
            varsOut.append(mapVars[var])
        elif (var.lower() in mapVars):
            varsOut.append(mapVars[var.lower()])
        else:
            varsOut.append(strip_varname(var))
    return varsOut

#-----------------------------------------------------------------------------
# remap gitm variables - Take the GITM variable names and expand them so
# that they could be used for things like publications.
#-----------------------------------------------------------------------------

def get_long_names(varsIn):

    mapVars = {
        'Rho' : 'Mass Density (km/m3)',
        '[O(3P)]': 'Neutral O Density (/m3)',
        '[O2]': 'Neutral O2 Density (/m3)',
        '[N2]': 'Neutral N2 Density (/m3)',
        '[N(4S)]': 'Neutral N Density (/m3)',
        '[NO]': 'Neutral NO Density (/m3)',
        '[He]': 'Neutral He Density (/m3)',
        '[N(2D)]': 'Neutral N(2D) Density(/m3)',
        '[N(2P)]': 'Neutral N(2P) Density (/m3)',
        '[H]': 'Neutral H Density (/m3)',
        '[CO2]': 'Neutral CO2 Density (/m3)',
        '[O(1D)]': 'Neutral O(1D) Density (/m3)',
        'Temperature': 'Neutral Temperature (K)',
        'Vn(east)': 'Neutral Eastward Velocity (m/s)',
        'Vn(north)': 'Neutral Northward Velocity (m/s)',
        'Vn(up)': 'Neutral Vertical Velocity (m/s)',
        'Vn(up,O(3P))': 'Vertical Velocity of O (m/s)',
        'Vn(up,O2)': 'Vertical Velocity of O2 (m/s)',
        'Vn(up,N2)': 'Vertical Velocity of N2 (m/s)',
        'Vn(up,N(4S))': 'Vertical Velocity of N (m/s)',
        'Vn(up,NO)': 'Vertical Velocity of NO (m/s)',
        'Vn(up,He)': 'Vertical Velocity of He (m/s)',
        '[O_4SP_+]': 'O+ Density (/m3)',
        '[NO+]': 'NO+ Density (/m3)',
        '[O2+]': 'O2+ Density (/m3)',
        '[N2+]': 'N2+ Density (/m3)',
        '[N+]': 'N+ Density (/m3)',
        '[O(2D)+]': 'O(2D+) Density (/m3)',
        '[O(2P)+]': 'O(2P)+ Density (/m3)',
        '[H+]': 'H+ Density (/m3)',
        '[He+]': 'He+ Density (/m3)',
        '[e-]': 'Electron Density (/m3)',
        'eTemperature': 'Te (K)',
        'iTemperature': 'Ti (K)',
        'Vi(east)': 'Vie (m/s)',
        'Vi(north)': 'Vin (m/s)',
        'Vi(up)': 'Viv (m/s)'}

    varsOut = []

    for var in varsIn:
        if (var in mapVars):
            varsOut.append(mapVars[var])
        else:
            varsOut.append(var)
    return varsOut

