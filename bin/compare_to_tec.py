#!/usr/bin/env python3
""" This is to compare output files to Madrigal TEC files
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse

from pyitm.fileio import util, madrigalio, gitmio, logfile
from pyitm.modeldata import satellite
from pyitm.general import time_conversion
from pyitm.plotting import plotutils
from pathlib import Path

# ----------------------------------------------------------------------------
# Get arguments as inputs into the code
#-----------------------------------------------------------------------------

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Compare Aether / GITM model results to Madrigal TEC')

    home_directory = str(Path.home())
    
    # Get the file that points to where to find the data:
    parser.add_argument('-lookup', default = home_directory + '/.pyitm_data', \
                        help = 'file that points to data locations')
    
    # Get the files to plot:
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')

    parser.add_argument('-v',  \
                        action='store_true', default = False, \
                        help = 'set verbose to true')
    
    parser.add_argument('-nomaps',  \
                        action='store_true', default = False, \
                        help = 'dont plot all of the maps')
    
    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------------
# Need a weighted average
# ----------------------------------------------------------------------------

def weighted_mean(var, weights):
    top = np.sum(var * weights)
    bot = np.sum(weights)
    ave = top / bot
    return ave

# ----------------------------------------------------------------------------
# Needed to run main script as the default executable from the command line
# ----------------------------------------------------------------------------

if __name__ == '__main__':
    
    # Get the input arguments
    args = get_args()
    filelist = args.filelist
    filelist = util.any_to_filelist(input_data = filelist)
    lookup = args.lookup
    verbose = args.v
    
    f1 = util.read_all_headers(filelist[0])
    f2 = util.read_all_headers(filelist[-1])

    # Read in 'gps' data. If we have vista, this may need to be changed.
    # This hands back a dict with keys for each data type, so only take 'gps'
    tecData = util.read_satfiles(satLookup=lookup, satname='gps', 
                                 startDate=f1['times'][0],
                                 endDate=f2['times'][0],
                                 verbose=verbose)['gps']

    # make some assumptions here about the TEC grid:
    nLats, nLons, nTimes = np.shape(tecData['tec'])
    dLat = 180.0/nLats
    dLon = 360.0/nLons
    lats1d = np.arange(-90.0 + dLat/2, 90.0, dLat)
    lons1d = (np.arange(0.0 + dLon/2, 360.0, dLon) + 180.0) % 360.0
    lons2d, lats2d = np.meshgrid(lons1d, lats1d)

    gitmData = util.read_all_files(filelist, varsToRead = 'tec', verbose = verbose)
    # need to reformulate the gitm data to fit into the canned functions:
    nTimesGitm, nLonsGitm, nLatsGitm = np.shape(np.squeeze(gitmData['tec']))

    # code needs a 1 for number of variables and number of altitudes:
    gitmData['data'] = gitmData['tec'].reshape(nTimesGitm, 1, nLonsGitm, nLatsGitm, 1)
    gitmData['vars'] = ['tec']

    # just say we have one alt at 100 km:
    gitmData['alts'] = 100.0

    iTimes = time_conversion.find_closest_times(tecData['times'], gitmData['times'])

    # need to re-arrange the data to use canned functions:
    for iT in iTimes:
        tec = tecData['tec'][:, :, iT]
        tecSave = tec[tec > 0]
        # get all valid points
        lons = lons2d[tec > 0]
        lats = lats2d[tec > 0]
        nPts = len(lons)
        # make a time array for all points
        t = np.array(tecData['times'][iT])
        ts = np.repeat(t, nPts)
        alts = np.zeros(nPts) + 100.0
        if (iT == iTimes[0]):
            tecInfo = {'times': ts,
                       'tec': tecSave,
                       'lons': lons,
                       'lats': lats,
                       'alts': alts}
            tecLocs = {'iStart': [0],
                       'iEnd': [nPts]}
        else:
            tecInfo['times'] = np.concatenate((tecInfo['times'], ts))
            tecInfo['tec'] = np.concatenate((tecInfo['tec'], tecSave))
            tecInfo['lons'] = np.concatenate((tecInfo['lons'], lons))
            tecInfo['lats'] = np.concatenate((tecInfo['lats'], lats))
            tecInfo['alts'] = np.concatenate((tecInfo['alts'], alts))
            tecLocs['iStart'].append(tecLocs['iEnd'][-1])
            tecLocs['iEnd'].append(len(tecInfo['times']))
            
    
    gitmTEC = satellite.extract_1d(tecInfo, gitmData, extrapolate=False, verbose=verbose, interpVar=[0])

    dataMinMax = plotutils.get_min_max_data(tecInfo['tec'], None, \
                                            color = 'red', stdevFac = 7.0)

    dpi = 150
    americanCutData = np.zeros((nTimesGitm, nLats))
    americanCutModel = np.zeros((nTimesGitm, nLats))
    globalDiff = np.zeros(nTimesGitm)
    globalTecMean = np.zeros(nTimesGitm)
    globalModelMean = np.zeros(nTimesGitm)
    globalTecMedian = np.zeros(nTimesGitm)
    globalModelMedian = np.zeros(nTimesGitm)
    globalRms = np.zeros(nTimesGitm)
    globalnRms = np.zeros(nTimesGitm)
    americanDiff = np.zeros(nTimesGitm)
    americanTecMean = np.zeros(nTimesGitm)
    americanModelMean = np.zeros(nTimesGitm)
    americanRms = np.zeros(nTimesGitm)
    americannRms = np.zeros(nTimesGitm)
    for iT in range(len(iTimes)):

        # Get the data:
        iS = tecLocs['iStart'][iT]
        iE = tecLocs['iEnd'][iT]
        lons = tecInfo['lons'][iS:iE]
        lats = tecInfo['lats'][iS:iE]
        tecData = tecInfo['tec'][iS:iE]
        tecGitm = gitmTEC['model_tec'][iS:iE]
        # Calculate some quantities to save:
        weights = np.cos(lats * np.pi/180.0)
        ave = weighted_mean((tecData - tecGitm)**2, weights)
        globalRms[iT] = np.sqrt(ave)
        globalDiff[iT] = weighted_mean(tecData - tecGitm, weights)
        globalTecMean[iT] = weighted_mean(tecData, weights)
        globalModelMean[iT] = weighted_mean(tecGitm, weights)
        globalTecMedian[iT] = np.median(tecData)
        globalModelMedian[iT] = np.median(tecGitm)
        globalnRms[iT] = globalRms[iT] / globalTecMean[iT]

        if (not args.nomaps):
        
            # Make some plots:
            fig = plt.figure(figsize=(10, 10), dpi = dpi)
            # 1 is bottom, 2 is top
            ax1 = fig.add_axes([0.08, 0.05, 0.95, 0.45])
            ax2 = fig.add_axes([0.08, 0.55, 0.95, 0.45])
        
            con1 = ax1.scatter(lons, lats, c = tecGitm, \
                               cmap = dataMinMax['cmap'], \
                               vmin = dataMinMax['mini'], \
                               vmax = dataMinMax['maxi'])
            con2 = ax2.scatter(lons, lats, c = tecData, \
                               cmap = dataMinMax['cmap'], \
                               vmin = dataMinMax['mini'], \
                               vmax = dataMinMax['maxi'])
            sTime= tecInfo['times'][iS].strftime("%d %b %Y %H:%M:%S UT")
            ax1.set_aspect(1.0)
            ax1.set_xlabel('Longitude (deg)')
            ax1.set_ylabel('Latitude (deg)')
            ax1.set_title('Model TEC Results at ' + sTime)
            ax2.set_aspect(1.0)
            ax2.set_xlabel('Longitude (deg)')
            ax2.set_ylabel('Latitude (deg)')
            ax2.set_title('TEC Measurements at ' + sTime)
                          
            cbar1 = fig.colorbar(con1, ax = ax1, shrink = 0.5, pad = 0.02)
            cbar1.set_label('TEC Model Results', rotation=90)
            cbar2 = fig.colorbar(con2, ax = ax2, shrink = 0.5, pad = 0.02)
            cbar2.set_label('TEC Measurements', rotation=90)

            filenamePrefix = 'gitm_compare_tec_'
            sTimeOut = tecInfo['times'][iS].strftime('%Y%m%d_%H%M%S')
            outFile = filenamePrefix + sTimeOut + '.png'
            print(" ==> Writing file : ", outFile)
            fig.savefig(outFile, dpi = dpi)
            plt.close(fig)

        # Focus on the American sector and make some data there:
        centerLon = 285
        dLon = 15.0
        americanWeights = np.zeros(nLats)
        for iL, lat in enumerate(lats1d):
            # hack, in two steps:
            dlats = np.abs(lats - lat)
            latsSmall = lats[dlats < 0.1]
            if (len(latsSmall) > 1):
                lonsSmall = lons[(lats == lat)]
                dataSmall = tecData[(lats == lat)]
                modelSmall = tecGitm[(lats == lat)]
                dlons = np.abs(lonsSmall - centerLon)
                if (np.min(dlons) < dLon):
                    americanCutData[iT, iL] = \
                        np.mean(dataSmall[dlons <= dLon])
                    americanCutModel[iT, iL] = \
                        np.mean(modelSmall[dlons <= dLon])
                    americanWeights[iL] = len(dataSmall[dlons <= dLon]) * np.cos(lat * np.pi / 180.0)
        ave = weighted_mean((americanCutData[iT, :] - americanCutModel[iT, :])**2, americanWeights)
        americanRms[iT] = np.sqrt(ave)
        americanDiff[iT] = \
            weighted_mean(americanCutData[iT, :] - americanCutModel[iT, :], americanWeights)
        americanTecMean[iT] = \
            weighted_mean(americanCutData[iT, :], americanWeights)
        americanModelMean[iT] = \
            weighted_mean(americanCutModel[iT, :], americanWeights)
        americannRms[iT] = globalRms[iT] / americanTecMean[iT]
    
    # Plot stuff in the american sector:
    l2d, t2d = np.meshgrid(lats1d, gitmData['times'])
    fig = plt.figure(figsize=(10, 10), dpi = dpi)
    # 1 is bottom, 2 is top
    ax1 = fig.add_axes([0.1, 0.05, 0.90, 0.4])
    ax2 = fig.add_axes([0.1, 0.55, 0.90, 0.4])
    con1 = ax1.pcolor(t2d, l2d, americanCutModel, \
                       cmap = dataMinMax['cmap'], \
                       vmin = dataMinMax['mini'], \
                       vmax = dataMinMax['maxi'])
    con2 = ax2.pcolor(t2d, l2d, americanCutData, \
                       cmap = dataMinMax['cmap'], \
                       vmin = dataMinMax['mini'], \
                       vmax = dataMinMax['maxi'])
    ax1.set_ylabel('Latitude (deg)')
    ax1.set_title('Model TEC Results in American Sector')
    ax2.set_ylabel('Latitude (deg)')
    ax2.set_title('TEC Measurements in American Sector')

    sTime = tecInfo['times'][0].strftime("%d %b %Y %H:%M:%S UT")
    eTime = tecInfo['times'][-1].strftime("%d %b %Y %H:%M:%S UT")
    ax1.set_xlabel('Time from ' + sTime + ' to ' + eTime)
    
    cbar1 = fig.colorbar(con1, ax = ax1, shrink = 0.5, pad = 0.02)
    cbar1.set_label('TEC Model Results', rotation=90)
    cbar2 = fig.colorbar(con2, ax = ax2, shrink = 0.5, pad = 0.02)
    cbar2.set_label('TEC Measurement Data', rotation=90)
    filenamePrefix = 'gitm_compare_tec_american_'
    sTimeOut = tecInfo['times'][0].strftime('%Y%m%d')
    outFile = filenamePrefix + sTimeOut + '.png'
    print(" ==> Writing file : ", outFile)
    fig.savefig(outFile, dpi = dpi)
    plt.close(fig)

    # Make summary plots of global means:
    
    fig = plt.figure(figsize=(10, 10), dpi = dpi)
    # 1 is bottom, 2 is top
    ax1 = fig.add_axes([0.1, 0.05, 0.85, 0.4])
    ax2 = fig.add_axes([0.1, 0.55, 0.85, 0.4])
    # top plot:
    ax2.plot(gitmData['times'], \
             globalTecMean, color = 'b', label = 'TEC Measurement')
    ax2.plot(gitmData['times'], \
             globalModelMean, color = 'r', label = 'TEC Model Results')
    ax2.plot(gitmData['times'], \
             globalTecMedian, color = 'b', \
             linestyle = ':', label = 'TEC Measurement (Median)')
    ax2.plot(gitmData['times'], \
             globalModelMedian, color = 'r', \
             linestyle = ':', label = 'TEC Model Results (Median)')
    ax1.plot(gitmData['times'], \
               globalRms, color = 'c', label = 'RMS Difference')
    ax1.plot(gitmData['times'], \
               globalDiff, color = 'k', label = 'Difference')
    ax1.axhline(0.0, linestyle = ':', color = 'grey')
                   
    ax1.set_ylabel('Global TEC Differences (TECU)')
    ax1.set_title('Mean and RMS Differences Across Globe')
    ax2.set_ylabel('TEC Means (TECU)')
    ax2.set_title('TEC Means/Medians Across Globe')
    ax1.legend()
    ax2.legend()
    
    sTime = tecInfo['times'][0].strftime("%d %b %Y %H:%M:%S UT")
    eTime = tecInfo['times'][-1].strftime("%d %b %Y %H:%M:%S UT")
    ax1.set_xlabel('Time from ' + sTime + ' to ' + eTime)

    # Write out a log file so other plots can be made at a later time:
    filenamePrefix = 'gitm_compare_tec_'
    sTimeOut = tecInfo['times'][0].strftime('%Y%m%d')
    outFile = filenamePrefix + sTimeOut + '.png'
    print(" ==> Writing file : ", outFile)
    fig.savefig(outFile, dpi = dpi)
    plt.close(fig)

    logData = {'times' : gitmData['times'],
               'TEC_Mean_Data': globalTecMean,
               'TEC_Mean_Model': globalModelMean,
               'TEC_Median_Data': globalTecMedian,
               'TEC_Median_Model': globalModelMedian,
               'RMS': globalRms,
               'nRMS': globalnRms,
               'Diff': globalDiff,
               'American_Mean_Data': americanTecMean,
               'American_Mean_Model': americanModelMean,
               'americanRMS': americanRms,
               'AmericannRMS': americannRms,
               'AmericanDiff': americanDiff}

    logfile.write_log(logData, fileHeader = 'gitm_tec')
    
