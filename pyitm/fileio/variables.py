#!/usr/bin/env python3

#-----------------------------------------------------------------------------
# remap gitm variables
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
# remap gitm variables
#-----------------------------------------------------------------------------

def get_short_names(varsIn):

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
        'Temperature': 'Tn',
        'Vn(east)': 'Ve',
        'Vn(north)': 'Vn',
        'Vn(up)': 'Vv',
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
        'eTemperature': 'Te',
        'iTemperature': 'Ti',
        'Vi(east)': 'Vie',
        'Vi(north)': 'Vin',
        'Vi(up)': 'Viv'}

    varsOut = []

    for var in varsIn:
        if (var in mapVars):
            varsOut.append(mapVars[var])
        else:
            varsOut.append(var)
    return varsOut

#-----------------------------------------------------------------------------
# remap gitm variables
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

