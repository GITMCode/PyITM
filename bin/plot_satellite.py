import argparse, os, scipy
import numpy as np
from pyitm.fileio import variables
from pyitm.fileio.util import any_to_filelist, read_all_headers, read_all_files, read_satfiles
from pyitm.modeldata import satellite
from pyitm.general import time_conversion
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(
        description="Script to interpolate model outputs to a satellite trajectory."
        "By default, both a plot and file with the interpolated data are created.")

    # Make separate group for required arguments for more useful help message
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-model', '--model_data', required=True, nargs='*',
                          default='./',
                          help='Path to model data, or glob pattern. (default = ./)')
    
    required.add_argument('-sat', '--sat_files', nargs='+',
                        help='Path to satellite files. As many as you want!\n'
                        'If a path is given that is not a file, we will recursively '
                        ' look for satellite files.')
    
    required.add_argument('-satmap', '--sat_lookup', default=None,
                        help='Path to csv file with lookup table on satellite:'
                        ' name, directory, and filename. Can use datetime strftime'
                        ' and glob patterns!')
    
    # optional arguments:
    parser.add_argument('-vars', '--variables', nargs='*', default=None,
                        help='Which variable(s) to plot? Default is 3 (rho), '
                        'just neutral density. If not set, all variables are '
                        'interpolated & written to file, but only rho is plotted.')
    
    parser.add_argument('-out', '--save_path', default='./',
                        help="Location where plot and interpolated data should be "
                        "saved. Default='./' (pwd)")
    
    parser.add_argument('-v', '--verbose', help='Print extra debug info?',
                        action='store_true')
    
    parser.add_argument('-np', '--noplot', action='store_true', 
                        help="Do not save plot, just data file?")
    
    parser.add_argument('-nd', '--nodata', action='store_true',
                        help="Do not save data, just make & save plot?\n\n")
    
    parser.add_argument("--sat_names", nargs='*',
                        help="If the script cannot infer the satellite names correctly,"
                             " specify them here")

    parser.add_argument('-winds', action='store_true', 
                        help='Plot the winds? Default is to only plot density. '
                        'Specifying this will turn off plotting all other columns.')
    
    args = parser.parse_args()

    return args


#-----------------------------------------------------------------------------  
def gridData4map(satDataDict, var='sat_rho',
                 lat_lim=None, dlat=2.5,
                 verbose=False):
    """
    Interpolate satellite data variable to a time/latitude grid for mapping.

    Separates satellite data into ascending and descending orbital passes at
    mid-latitudes, then interpolates the specified variable onto a regular
    time-latitude grid.

    Parameters
    ----------
    satDataDict : dict
        Dictionary containing satellite data with keys:
        - 'times' : array-like
            Datetime objects representing observation times
        - 'lats' : array-like
            Latitude values in degrees
        - **var** : array-like
            Variable data to interpolate (key name specified by var parameter)
        - 'orbital_period' : float
            Satellite orbital period in seconds
    var : str, optional
        Variable name to interpolate from satDataDict (default: 'sat_rho')
    lat_lim : float, optional
        Latitude limits for output grid in degrees. Set to the min/max observed by
        satellite if None (default: None)
    dlat : float, optional
        Latitude step for output grid in degrees (default: 2.5)
    verbose : bool, optional
        If True, print debug information (default: False)

    Returns
    -------
    dict
        Dictionary containing gridded data with keys:
        - 'desc' : numpy.ndarray
            Gridded data for descending orbital passes
        - 'asc' : numpy.ndarray
            Gridded data for ascending orbital passes
        - 'times' : numpy.ndarray
            Time coordinate array (hours from start)
        - 'lats' : numpy.ndarray
            Latitude coordinate array (degrees)
        - 'asc_n' : float
            Local time of ascending node in hours (0 if not found, orbit too short)
        - 'desc_n' : float
            Local time of descending node in hours (0 if not found)

    """

    # Make sure we have the orbital period
    if 'orbital_period' not in satDataDict.keys():
        if verbose:
            print("-> Griddata4map: No orbital period found. Calling calc_period()")
        satDataDict = satellite.calc_period(satDataDict, verbose=verbose)
        if 'orbital_period' not in satDataDict.keys():
            raise ValueError("No orbital period found. Cannot grid data for map.")

    # find where sat is descending & at midlatitudes

    if lat_lim is None:
        low_lat_lim = np.min(satDataDict['lats'])
        high_lat_lim = np.max(satDataDict['lats'])
        if verbose:
            print(f"-> Griddata4map: No lat_lim given. Setting to ({high_lat_lim:.1f},{low_lat_lim:.1f}) deg")
    else:
        low_lat_lim = -lat_lim
        high_lat_lim = lat_lim
        if verbose:
            print(f"-> Griddata4map: Setting lat_lim to +/- {lat_lim:.1f} deg")
            
    desc_indices = np.where((np.diff(satDataDict['lats']) < 0) & (
                            (satDataDict['lats'] < high_lat_lim)[1:]) & (
                            (satDataDict['lats'] > low_lat_lim)[1:]))[0]

    asc_indices = np.where((np.diff(satDataDict['lats']) > 0) & (
                           (satDataDict['lats'] < high_lat_lim)[1:]) & (
                           (satDataDict['lats'] > low_lat_lim)[1:]))[0]

    t0, t1 = min(satDataDict['times']), max(satDataDict['times'])

    Xs = np.arange(0, (t1 - t0).total_seconds()/3600.0, 
                   satDataDict['orbital_period'].total_seconds()/3600.0)
    Ys = np.arange(low_lat_lim, high_lat_lim+dlat, dlat)
    
    gridX, gridY = np.meshgrid(Xs, Ys, indexing='ij')

    ts = []
    for dt in satDataDict['times'] - t0:
        ts.append(dt.total_seconds()/3600.0)
    ts = np.array(ts)

    desc = scipy.interpolate.griddata(
        ((ts[desc_indices]),
         satDataDict['lats'][desc_indices]),
        satDataDict[var][desc_indices],
        (gridX, gridY),
        method='linear')

    asc = scipy.interpolate.griddata(
        ((ts[asc_indices]),
         satDataDict['lats'][asc_indices]),
        satDataDict[var][asc_indices],
        (gridX, gridY),
        method='linear')
    
    # get local time
    if 'lst' in satDataDict.keys():
        lts = satDataDict['lst']
    elif 'lsts' in satDataDict.keys():
        lts = satDataDict['lsts']
    else:
        lts = time_conversion.ut_to_lt(satDataDict['times'],
                                       satDataDict['lons'])
        
    # find where lat crosses 0 deg north going north (ascending node)
    an_idx = np.where((np.diff(satDataDict['lats']) > 0) 
                      & ((satDataDict['lats'][:-1] < 0)))[0]
    if len(an_idx) > 0:
        an = np.median(lts[an_idx])
    else:
        an = 0
        if verbose:
            print("-> Griddata4map: No ascending node found.")
    
    dn_idx = np.where((np.diff(satDataDict['lats']) < 0) 
                      & ((satDataDict['lats'][:-1] > 0)))[0]
    if len(dn_idx) > 0:
        dn = np.median(lts[dn_idx])
    else:
        dn = 0
        if verbose:
            print("-> Griddata4map: No descending node found.")

    return dict(desc =desc,
                asc = asc,
                times = Xs,
                lats = Ys,
                asc_n = an,
                desc_n = dn)

#-----------------------------------------------------------------------------
def plot_asc_desc(satDict, satName, varname='rho',
                  savepath='./', verbose=False):
    """
    Makes a lat vs time contour plot of satellite data & model data. Saves.
    
    """

    satvarname = 'sat_' + varname
    modelvarname = 'model_' + varname

    satGridded = gridData4map(satDict, var=satvarname, verbose=verbose)
    modelGridded = gridData4map(satDict, var=modelvarname, verbose=verbose)

    for direction in ['desc', 'asc']:
        # for the title - direction name
        dirname = direction.title() + 'ending'

        fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

        vmin = min(np.nanmin(satGridded[direction]), np.nanmin(modelGridded[direction]))
        vmax = max(np.nanmax(satGridded[direction]), np.nanmax(modelGridded[direction]))
        if vmin > 0:
            cmap='plasma'
        else:
            cmap='bwr'
            # make sure limits are symmetric about 0
            minlim = min(np.abs([vmin, vmax]))
            vmin = -minlim
            vmax = minlim

        pcm = axs[0].pcolormesh(satGridded['times'], satGridded['lats'], 
                                satGridded[direction].T,
                                cmap=cmap, vmin=vmin, vmax=vmax)
        
        axs[0].set_ylabel('Latitude (deg)')
        if f"{direction}_n" in satGridded.keys():
            axs[0].set_title(satName.upper() + " " + dirname +
                             " Node (LT  : %4.1f hours)" % satGridded[f"{direction}_n"])
        else:
            axs[0].set_title(satName.upper() + " " + dirname)
        fig.colorbar(pcm, ax=axs[0], label=varname + ' (kg/m3)')

        pcm = axs[1].pcolormesh(satGridded['times'], satGridded['lats'], 
                                modelGridded[direction].T,
                                cmap=cmap, vmin=vmin, vmax=vmax)
        
        axs[1].set_ylabel('Latitude (deg)')
        axs[1].set_title("Modeled " + dirname)
        axs[1].set_xlabel('Time (hours from ' + 
                        satDict['times'][0].strftime('%Y-%m-%d %H:%M') + ')')
        fig.colorbar(pcm, ax=axs[1], label=varname + ' (kg/m3)')

        outdate = satDict['times'][0].strftime('%Y%m%d')
        savename = os.path.join(savepath, 
                                f"{satName}_{varname}_{dirname}_{outdate}.png")

        print('--> Saving plot as: ' + savename)
        plt.savefig(savename)
        plt.close()


    return



def makesatplot(satDataDict, satName, varname, savepath, verbose=False):
    """
    Makes line plots of raw & % difference of satellite density vs model density. Saves.
    
    Parameters
    ----------
    satDataDict : dict
        Dictionary containing satellite data with keys:
        - 'times' : array-like
            Datetime objects representing observation times
        - 'sat_[varname]' : array-like
            Satellite variable data (e.g., 'sat_rho' for density)
        - 'model_[varname]' : array-like
            Corresponding model variable data (e.g., 'model_rho' for density)
        - Optional: 'smoothed_sat_[varname]' and 'smoothed_model_[varname]'
            Smoothed versions of the satellite and model variable data
    satName : str
        Name of the satellite (e.g., 'CHAMP')
    varname : str
        Variable name to plot (e.g., 'rho' for density)
    savepath : str
        Directory path where the plot will be saved
    verbose : bool, optional
        If True, print debug information (default: False)

    Returns
    -------
    None

    """

    svar = 'sat_' + varname
    mvar = 'model_' + varname
    if svar not in satDataDict.keys() or mvar not in satDataDict.keys():
        print(satDataDict.keys())
        raise ValueError(f"satDataDict must contain keys '{svar}' and '{mvar}'")
    svars = 'smoothed_' + svar
    mvars = 'smoothed_' + mvar

    fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    didSmooth = svars in satDataDict.keys() and mvars in satDataDict.keys()
    
    if didSmooth:
        axs[0].plot(satDataDict['times'], satDataDict[mvar], 'b',
                    alpha=0.1)
        axs[0].plot(satDataDict['times'], satDataDict[svar], 'r',
                    alpha=0.1)
        axs[0].plot(satDataDict['times'], satDataDict[mvars], 'b',
                    label="Modeled Density")
        axs[0].plot(satDataDict['times'], satDataDict[svars], 'r',
                    label=satName.upper() + " Density")
    else:
        axs[0].plot(satDataDict['times'], satDataDict[mvar], 'b',
                    label=f"Modeled {mvar}")
        axs[0].plot(satDataDict['times'], satDataDict[svar], 'r',
                    label=satName.upper() + f" {svar}")
    axs[0].legend()

    
    if didSmooth:
        axs[1].plot(satDataDict['times'],
                    100*(satDataDict[mvar] - satDataDict[svar])
                    /satDataDict[svar],
                    alpha=0.1, color='k')
        axs[1].plot(satDataDict['times'],
                    100*(satDataDict[mvars]-satDataDict[svars])
                    /satDataDict[svars],
                    color='k')
    else:
        axs[1].plot(satDataDict['times'],
                    100*(satDataDict[mvar]- satDataDict[svar])
                    /satDataDict[svar],
                    color='k')
    
    axs[1].hlines(0, min(satDataDict['times']), max(satDataDict['times']),
                  linestyle='--', color='k')

    axs[1].set_xlabel(satDataDict['times'][0].strftime('%b %d, %Y %H:%M UT') + ' - ' + \
                      satDataDict['times'][-1].strftime('%b %d, %Y %H:%M UT'))

    ylabel = variables.get_long_names([varname])[0]
    if verbose and ylabel == varname:
            print(f" variables.get_long_names() returned no match for {varname}")
    axs[0].set_ylabel(ylabel)
    axs[1].set_ylabel(' Diff (%)')
    fig.suptitle(satName.upper())

    outdate = satDataDict['times'][0].strftime('%Y%m%d')
    savename = os.path.join(savepath, f"{satName}_{varname}_line_{outdate}.png")
        
    print('--> Saving plot as: ' + savename)
    plt.tight_layout()
    plt.savefig(savename)

    return

def write_out_data(satDataDict, savepath, varname, satName, verbose=False):
    """
    Writes SatDataDict (from plot_satellite / extract1d to a csv file)
    """

    outdate = satDataDict['times'][0].strftime('%Y%m%d')
    savename = os.path.join(savepath, f"{satName}_{varname}_compare_{outdate}.csv")

    if os.path.exists(savename) and verbose:
        print(f"-> File {savename} exists. Overwriting!")
    else:
        print('--> Saving Data file as: ' + savename)
    
    with open(savename, 'w') as f:
        # Header
        f.write(f"# Satellite name: {satName}\n")
        f.write(f"# Orbital Period: {satDataDict.pop('orbital_period', None)}\n")
        f.write("\n")

        # First line
        f.write(','.join([k for k in satDataDict.keys()]) +'\n')
        for iline in range(len(satDataDict['times'])):
            f.write(','.join([str(satDataDict[k][iline]) for k in satDataDict.keys()]) +'\n')
    
    return
            


def main(satfiles, modeldatapath, vars2plot=3, satNames=None, sat_lookup=None,
         isWind=False,
         saveData=True, savePlot=True, savePath='./',
         verbose=False):
    
    if not saveData and not savePlot: # Save compute cycles
        print("Not saving data & not saving plot. Nothing to do here.")
        return
    
    # Check if savepath exists. Make it if not. 
    # Good to leave early & save time if there was a typo
    if not os.path.isdir(savePath):
        # This will error if the parent directories aren't found
        os.mkdir(savePath)
    
    # We should make sure the model files exist.
    model_filelist = any_to_filelist(modeldatapath)
    # Needs to be done this way in case we have netcdf files (can't read multiple headers)
    header0 = read_all_headers(model_filelist[0], verbose=verbose)
    header1 = read_all_headers(model_filelist[-1], verbose=verbose)
    # And we need the dates to read in satellite files!
    date0, date1 = header0['times'][0], header1['times'][-1]
    if verbose:
        print(f"Model data ranges from {date0} to {date1}")

    # Read in satellite data
    print("Reading in satellite data...")

    satData = read_satfiles(satfiles, 
                            satname=satNames,
                            satLookup=sat_lookup,
                            startDate=date0, 
                            endDate=date1,
                            verbose=verbose)
    
    print(f"Found data from satellites: {list(satData.keys())}.")

    # If we're plotting winds, only read ve/vn from model data:
    if isWind:
        vars2plot = []
        # 2 versions of this - one for us to iterate over later, one to tell the reader what to read
        winds2plot = []
        # Check if we have neutral winds...
        neutralwind = False
        ionwind = False #weird name
        for satname in satData.keys():
            if 'Ve' in satData[satname].keys() and 'Ve' in satData[satname].keys():
                neutralwind = True
            if 'Vie' in satData[satname].keys() and 'Vie' in satData[satname].keys():
                ionwind = True
        if neutralwind:
            winds2plot.append(['Ve', 'Vn'])
            vars2plot.extend(['Ve', 'Vn'])
        if ionwind:
            winds2plot.append(['Vie', 'Vin'])
            vars2plot.extend(['Vie', 'Vin'])

    # Read in the first model data file, which allows the model reader to handle the
    # variable name remapping instead of us. Allows for more flexibility with TEC or
    # multple models.
    modelt0 = read_all_files(model_filelist[0], varsToRead=vars2plot, verbose=verbose)
    varnames = modelt0['shortname']
    # varnames.extend(modelt0['vars'])

    vars_found = []

    for varname in varnames:
        if verbose:
            print(f"-> looking for Variable {varname}")
        for sat in satData.keys():
            if varname in satData[sat].keys():
                if verbose:
                    print(f"-> Variable {varname} found in {sat}.")
                if varname not in vars_found:
                    vars_found.append(varname)
    if verbose:
        print(f"Variables to plot: {vars_found}")

    ############################
    # Now read in all the model data & do the plotting
    ############################

    if isWind:
        print(f"-> Reading in model data for {vars2plot}. This may take a while...")
            # Read in model data. We pass handling of the variable names to read_all_files
        if len(modeldatapath) == 1:
            modeldatapath = modeldatapath[0]
        modelData = read_all_files(modeldatapath, varsToRead=vars2plot, verbose=verbose)
        for name, data in satData.items():
            # will be some combo/order of [ion velo, neutral velo]
            for windvars in winds2plot:
                # Make sure satellite has the wind
                if all([vname in data for vname in windvars]):

                    # Get indices of wind in model data
                    interp_ind = [varnames.index(vname) for vname in windvars]

                    interp_data = satellite.extract_1d(data, modelData, interpVar=interp_ind, 
                                                       verbose=verbose)
                    windData = satellite.calc_zon_merid_wind(interp_data)

                    #Call plotting functions. We have already determined if ion/neutral are present.
                    for varname in windvars:
                        if savePlot:
                            plot_asc_desc(windData, satName=name, varname=varname,
                                            savepath=savePath, verbose=verbose)
                        if saveData:
                            write_out_data(windData, savepath=savePath, varname=varname,
                                        satName=name, verbose=verbose)


    else:
        for varname in vars_found:

            print(f"-> Reading in model data for variable: {varname}. This may take a while...")
            # Read in model data. We pass handling of the variable names to read_all_files
            if len(modeldatapath) == 1:
                modeldatapath = modeldatapath[0]
            modelData = read_all_files(modeldatapath, varsToRead=varname, verbose=verbose)

            print("Model data has been read in! Continuing")

            for name, data in satData.items():

                if varname not in data.keys():
                    if verbose:
                        print(f"-> Variable {varname} not found in satellite data for {name}. Skipping...")
                    continue
                if verbose:
                    print(f"Processing satellite: {name} for variable: {varname}")

                sat_data_tmp = satellite.extract_1d(data, modelData, verbose=verbose,
                                                    interpVar=None)
                outData = satellite.orbit_average(sat_data_tmp,
                                                varlist=['sat_'+varname, 'model_'+varname],
                                                verbose=verbose)

                if savePlot:
                    makesatplot(outData, satName=name, varname=varname,
                                savepath=savePath, verbose=verbose)
                    plot_asc_desc(outData, satName=name, varname=varname,
                                savepath=savePath, verbose=verbose)
                if saveData:
                    write_out_data(outData, savepath=savePath, varname=varname,
                                satName=name, verbose=verbose)



if __name__ == "__main__":
    args = parse_args()

    main(args.sat_files, args.model_data, args.variables, args.sat_names,
         sat_lookup=args.sat_lookup,
         isWind=args.winds,
         saveData=(not args.nodata), 
         savePlot=(not args.noplot),
         savePath=args.save_path,
         verbose=args.verbose)