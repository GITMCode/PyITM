import argparse
import os
import numpy as np
from pyitm.fileio import satelliteio
from pyitm.fileio.util import read_all_files
from pyitm.modeldata import satellite
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(
        description="Script to interpolate model outputs to a satellite trajectory."
        "By default, both a plot and file with the interpolated data are created.")

    # Make separate group for required arguments for more useful help message
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-model', '--model_data', required=True, nargs='+',
                        help='Path to model data, or glob pattern')
    
    required.add_argument('-sat', '--sat_files', required=True, nargs='+',
                        help='Path to satellite files. As many as you want!')
    
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
    
    args = parser.parse_args()

    return args


def makesatplot(satDataDict, savepath, verbose=False, saveName=None):

    fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    didSmooth = 'smoothed_sat_rho' in satDataDict.keys() 
    
    if didSmooth:
        axs[0].plot(satDataDict['times'], satDataDict['model_rho'], 'b', 
                    alpha = 0.1)
        axs[0].plot(satDataDict['times'], satDataDict['sat_rho'], 'r', 
                    alpha = 0.1)
        axs[0].plot(satDataDict['times'], satDataDict['smoothed_model_rho'], 'b', 
                    label="Modeled Density")
        axs[0].plot(satDataDict['times'], satDataDict['smoothed_sat_rho'], 'r', 
                    label=satDataDict['sat_name'].upper() + " Density")
    else:
        axs[0].plot(satDataDict['times'], satDataDict['model_rho'], 'b', 
                    label="Modeled Density")
        axs[0].plot(satDataDict['times'], satDataDict['sat_rho'], 'r', 
                    label=satDataDict['sat_name'].upper() + " Density")
    axs[0].legend()

    
    if didSmooth:
        axs[1].plot(satDataDict['times'], 
                    100*(satDataDict['sat_rho']- satDataDict['model_rho'])
                    /satDataDict['sat_rho'],
                    alpha = 0.1, color='k')
        axs[1].plot(satDataDict['times'], 
                    100*(satDataDict['smoothed_sat_rho']-satDataDict['smoothed_model_rho'])
                    /satDataDict['smoothed_sat_rho'],
                    color='k')
    else:
        axs[1].plot(satDataDict['times'], 
                    100*(satDataDict['sat_rho']- satDataDict['model_rho'])
                    /satDataDict['sat_rho'],
                    color='k')
    
    axs[1].hlines(0, min(satDataDict['times']), max(satDataDict['times']),
                  linestyle='--', color='k')

    axs[1].set_ylabel(' Density (kg/m3)')
    axs[1].set_ylabel(' Diff (%)')
    fig.suptitle(satDataDict['sat_name'].upper())

    i = 0
    outdate = satDataDict['times'][0].strftime('%Y%m%d')
    savename = os.path.join(savepath, f"{satDataDict['sat_name']}_density_{outdate}_000.png")
    
    old_end = str(i).rjust(3, '0') + '.png'
    while os.path.exists(savename):
        # start adding zeros!
        i += 1
        new_end = str(i).rjust(3, '0') + '.png'
        savename = savename.replace(old_end, new_end)
        old_end = new_end
    
    print('--> Saving plot as: ' + savename)
    plt.savefig(savename)

    return

def write_out_data(satDataDict, savepath, verbose=False):
    """
    Writes SatDataDict (from plot_satellite / extract1d to a csv file)
    """


    i = 0
    outdate = satDataDict['times'][0].strftime('%Y%m%d')
    savename = os.path.join(savepath, f"{satDataDict['sat_name']}_{outdate}_000.csv")
    
    old_end = str(i).rjust(3, '0') + '.csv'
    while os.path.exists(savename):
        # start adding zeros!
        i += 1
        new_end = str(i).rjust(3, '0') + '.csv'
        savename = savename.replace(old_end, new_end)
        old_end = new_end


    print('--> Saving Data file as: ' + savename)
    with open(savename, 'w') as f:
        # Header
        f.write(f"# Satellite name: {satDataDict.pop('sat_name', None)}\n")
        f.write(f"# Orbital Period: {satDataDict.pop('orbital_period', None)}\n")
        f.write("\n")

        # First line
        f.write(','.join([k for k in satDataDict.keys()]) +'\n')
        for iline in range(len(satDataDict['times'])):
            # l = ''
            # for k in satDataDict.keys():
            #     l += str(satDataDict[k][iline]) + ','
            # l += '\n'
            # print(l)
            # f.write(l)
            f.write(','.join([str(satDataDict[k][iline]) for k in satDataDict.keys()]) +'\n')
    
    return
            


def main(satfiles, modeldatapath, vars2plot=3, satNames=None, 
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

    if satNames is not None:
        # Ensure the length of satnames and satfiles are the same
        if len(satNames) == 1:
            satNames *= len(satfiles)
        if len(satNames) != len(satfiles):
            raise ValueError("Length of satellite names & files do not match!")
    else:
        # Let reader infer
        satNames = [None] * len(satfiles)

    # Read in satellite data
    satData = []
    for satFileName, satName in zip(satfiles, satNames):
        if verbose:
            print(f"Attempting to read satellite file: '{satFileName}'")
            print(f" with name '{satName}'" if satName is not None else "")

        satData.append(satelliteio.read_sat_file(satFileName, satname=satName,
                                                 verbose=verbose))

    # Read in model data
    print("Reading in model data. This may take a while...")
    modelData = read_all_files(modeldatapath, varsToRead=vars2plot, verbose=verbose)

    print("Model data has been read in! Continuing")

    outData = []
    for oneSat in satData:
        sat_data_tmp = satellite.extract_1d(oneSat, modelData, verbose=verbose,
                                            interpVar=vars2plot)
        outData.append(satellite.orbit_average(sat_data_tmp,
                                                verbose=verbose))

    for oneSat in outData:
        if savePlot:
            makesatplot(oneSat, savepath=savePath, verbose=True)
        if saveData:
            write_out_data(oneSat, savepath=savePath, verbose=True)



if __name__ == "__main__":
    args = parse_args()

    main(args.sat_files, args.model_data, args.variables, args.sat_names,
         saveData=(not args.nodata), 
         savePlot=(not args.noplot),
         savePath=args.save_path,
         verbose=args.verbose)