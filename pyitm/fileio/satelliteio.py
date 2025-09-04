import numpy as np
from datetime import datetime

def calc_wind_dir(lons, lats):

    dlats = lats[1:] - lats[:-1]
    dlats = np.concatenate((dlats, [dlats[-1]]))

    dlons = lons[1:] - lons[:-1]
    dlons = np.concatenate((dlons, [dlons[-1]]))

    # Longitude can go across the 0 - 360 or 360 - 0 boundary, so
    # we need to correct for this possibility:
    
    dlons[dlons > 180.0] = dlons[dlons > 180.0] - 360.0
    dlons[dlons < -180.0] = dlons[dlons < -180.0] + 360.0

    # Longitudes get closer together near the poles, so we need to
    # correct for that also:
    dlons = dlons * np.cos(lats * np.pi / 180.0)

    # Make a unit vector of the direction of travel:
    mag = np.sqrt(dlats**2 + dlons**2)
    unitn = dlats / mag 
    unite = dlons / mag

    # Rotate the unit vector, so that it points orthogonal to the
    # orbit plane.  This will be the actual wind vector direction:
    uniteRot = unitn
    unitnRot = -unite

    return uniteRot, unitnRot

def _read_goce(file):
    """
    Read a GOCE file

    Inputs
    ------
        file (str) - Path to GOCE file

    Returns
    -------
        (dict) - GOCE data
    
    """

    data = {}
    data["times"] = []
    data["alts"] = []
    data["lats"] = []
    data["lons"] = []
    data["lst"] = []
    data["density"] = []
    data["Ve"] = []
    data["Vn"] = []
    data["Vr"] = []
    data["densityError"] = []
    data["windError"] = []
    data["FlagOver"] = []
    data["FlagEclipse"] = []
    data["FlagAD"] = []
    data["FlagThuster"] = []

    f = open(file, 'r')

    for line in f:

        if (line.find('#') < 0):
            items = line.split()
            ymd = items[0].split('-')
            hms = items[1].split(':')
            s = float(hms[2])
            data["times"].append(datetime(int(ymd[0]),int(ymd[1]),int(ymd[2]),
                                         int(hms[0]),int(hms[1]),int(s)))
            data["alts"].append(float(items[3])/1000.0)
            data["lons"].append(float(items[4]))
            data["lats"].append(float(items[5]))
            data["lst"].append(float(items[6]))
            data["density"].append(float(items[8]))
            data["Ve"].append(float(items[9]))
            data["Vn"].append(float(items[10]))
            data["Vr"].append(float(items[11]))
            data["densityError"].append(float(items[12]))
            data["windError"].append(float(items[13]))
            data["FlagOver"].append(int(items[14]))
            data["FlagEclipse"].append(int(items[15]))
            data["FlagAD"].append(int(items[16]))
            data["FlagThuster"].append(int(items[17]))

    f.close()

    # Here we are calculating the direction of travel of the sat:

    uniteRot, unitnRot = calc_wind_dir(np.array(data["lons"]),
                                       np.array(data["lats"]))

    data["wind_e_dir"] = uniteRot
    data["wind_n_dir"] = unitnRot
    
    return data

#-----------------------------------------------------------------------------
# Read CHAMP data
#-----------------------------------------------------------------------------

def _read_champ(file):
    """
    Read a CHAMP file

    Inputs
    ------
        file (str) - Path to CHAMP file

    Returns
    -------
        (dict) - CHAMP data
    
    """
    
    data = {}
    data["times"] = []
    data["alts"] = []
    data["lats"] = []
    data["lons"] = []
    data["lst"] = []
    data["density"] = []

    f = open(file, 'r')

    for line in f:

        if (line.find('#') < 0):
            items = line.split()
            ymd = items[0].split('-')
            hms = items[1].split(':')
            s = float(hms[2])
            data["times"].append(datetime(int(ymd[0]),int(ymd[1]),int(ymd[2]),
                                         int(hms[0]),int(hms[1]),int(s)))
            data["alts"].append(float(items[3])/1000.0)
            data["lons"].append( (float(items[4])+360.0) % 360.0 )
            data["lats"].append(float(items[5]))
            data["lst"].append(float(items[6]))
            data["density"].append(float(items[8]))

    f.close()

    return data

def _read_grace(file):
    """
    Read a GRACE/GRACE-FO file

    Inputs
    ------
        file (str) - Path to GRACE file

    Returns
    -------
        (dict) - GRACE data
    
    """
    
    data = {}
    data["times"] = []
    data["alts"] = []
    data["lats"] = []
    data["lons"] = []
    data["lst"] = []
    data["arglat"] = []
    data["density_mean"] = []
    data["density_flag"] = []
    data["density_mean_flag"] = []

    f = open(file, 'r')

    for line in f:

        if (line.find('#') < 0):
            items = line.split()
            ymd = items[0].split('-')
            hms = items[1].split(':')
            s = float(hms[2])
            data["times"].append(datetime(int(ymd[0]),int(ymd[1]),int(ymd[2]),
                                         int(hms[0]),int(hms[1]),int(s)))
            data["alts"].append(float(items[4])/1000.0)
            data["lons"].append( (float(items[5])+360.0) % 360.0 )
            data["lats"].append(float(items[6]))
            data["lst"].append(float(items[7]))
            data["arglat"].append(float(items[8]))
            data["density"].append(float(items[9]))
            data["density_mean"].append(float(items[10]))
            data["density_flag"].append(float(items[11]))
            data["density_mean_flag"].append(float(items[12]))

    f.close()

    return data


def _read_grace_winds(file):
    """
    Read a GRACE/GRACE-FO (wind) file

    Inputs
    ------
        file (str) - Path to GRACE file

    Returns
    -------
        (dict) - GRACE data. 
            - crosswind_speed is the magnitude
            - Uv is unit vector in each of the 3 directions
    
    """
    
    data = {}
    data["times"] = []
    data["alts"] = []
    data["lats"] = []
    data["lons"] = []
    data["lst"] = []
    data["arglat"] = []
    data["crosswind_speed"] = []
    data["uVnorth"] = []
    data["uVeast"] = []
    data["uVdown"] = []
    data["validity_flag"] = []

    f = open(file, 'r')

    for line in f:

        if (line.find('#') < 0):
            items = line.split()
            ymd = items[0].split('-')
            hms = items[1].split(':')
            s = float(hms[2])
            data["times"].append(datetime(int(ymd[0]),int(ymd[1]),int(ymd[2]),
                                         int(hms[0]),int(hms[1]),int(s)))
            data["alts"].append(float(items[4])/1000.0)
            data["lons"].append( (float(items[5])+360.0) % 360.0 )
            data["lats"].append(float(items[6]))
            data["lst"].append(float(items[7]))
            data["arglat"].append(float(items[8]))
            data["crosswind_speed"].append(float(items[9]))
            data["uVnorth"].append(float(items[10]))
            data["uVeast"].append(float(items[11]))
            data["uVdown"].append(float(items[12]))
            data["validity_flag"].append(float(items[12]))

    f.close()

    return data




def _read_champ_winds(file):
    """
    Read a CHAMP wind file

    Parameters
    ----------
        file (str) - Path to CHAMP file

    Returns
    -------
        (dict) - CHAMP data
    
    """

    data = {}
    data["times"] = []
    data["alts"] = []
    data["lats"] = []
    data["lons"] = []
    data["lst"] = []
    data["Ve"] = []
    data["Vn"] = []
    data["Vr"] = []
    data["density"] = []

    f = open(file, 'r')

    for line in f:

        if (line.find('#') < 0):
            items = line.split()
            ymd = items[0].split('-')
            hms = items[1].split(':')
            s = float(hms[2])
            data["times"].append(datetime(int(ymd[0]),int(ymd[1]),int(ymd[2]),
                                         int(hms[0]),int(hms[1]),int(s)))
            data["alts"].append(float(items[3])/1000.0)
            data["lons"].append((float(items[4]) + 360.0) % 360.0)
            data["lats"].append(float(items[5]))
            data["lst"].append(float(items[6]))
            data["Ve"].append(float(items[8]))
            data["Vn"].append(float(items[9]))

    f.close()

    # Here we are calculating the direction of travel of the sat:

    uniteRot, unitnRot = calc_wind_dir(np.array(data["lons"]),
                                       np.array(data["lats"]))

    data["wind_e_dir"] = uniteRot
    data["wind_n_dir"] = unitnRot

    return data


def read_sat_file(filename:str, satname=None, verbose=False):
    """
    Generic reader for any satellite file.

    Will infer satellite name and use the correct reader

    Parameters
    ----------
        filename (str): Path to satellite data
        verbose (bool, False): Print extra debugging information?

    Returns
    -------
        (dict): Dictionary containing the satellite data

    Raises
    ------
        ValueError: satellite name wasn't inferred. Either multiple or no matches

    Notes
    -----
    - Satname is optional! If provided, that reader is used. Otherwise it is
        inferred from filename
    
    """

    
    # Define a lookup table for satellite names and patterns.
    # This avoids us having to hardcode if/else for each new satellite type
    # To add a new reader, define it then modify satreaders & satlookup

    # The satellite name & the function used to call it
    satreaders = {'grace_density': _read_grace,
                  'grace_wind': _read_grace_winds,
                  'goce': _read_goce,
                  'champ': _read_champ}
    # satellite name & patterns that should be checked against filename
    satlookup = {'goce': ['go'],
                 'champ': ['ch'],
                 'grace_density': ['gr_dns', 'ga_dns', 'gb_dns', 'gc_dns'],
                 'grace_wind': ['gr_wnd', 'ga_wnd', 'gb_wnd', 'gc_wnd'],
                 }

    if satname is None:
        # Infer satellite name
        # Make sure we only look at the filename
        if '/' in filename:
            sat_filename = filename.split('/')[-1].lower()
        else:
            sat_filename = filename.lower()

        # See if it matches any of the patterns we expect.
        # Use list to check if there are multiple matches
        satnames = []
        for name in satlookup.keys():
            for pattern in satlookup[name]:
                if pattern in sat_filename:
                    satnames.append(name)
                    if verbose:
                        print(f"Match for {filename} found as {pattern}: {name}")

        if len(satnames) > 1:
            raise ValueError(
                f"Satellite file {filename} matches {satnames}. "
                "Maybe pass satellite name manually")
        elif len(satnames) < 1:
            raise ValueError(f"No reader found for satellite file: '{filename}'\n"
                             f"\t>> Supported satnames are: {list(satreaders.keys())}")
        else:
            satName = satnames[0]
    else:
        # We were handed a name. Try using the reader
        satName = satname.lower()
        if satName not in satreaders.keys():
            raise NotImplementedError(
                f"Reader for {satName} is not yet created! "
                f"Supported satname's are: {satreaders.keys()}"
            )
        if verbose:
            print(f"Attempting to read {filename} with reader for {satName}")

    # Dispatch the reader based on the inferred name:
    satData = satreaders[satName](filename)
    satData['sat_name'] = satName

    return satData

