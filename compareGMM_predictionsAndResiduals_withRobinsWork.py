# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 16:28:08 2018

@author: cde62

Computes response spectral ratios and smoothed fourier ratios for strong motion station in station.txt relative to reference station POTS

"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import matplotlib.lines as mlines
from obspy.signal import konnoohmachismoothing
import os.path
import pandas as pd
from glob import glob
from scipy.stats.mstats import gmean
from math import ceil, log
import eqsig

# define the run ID to differentiate runs
# runID = 'v1_wCPLB_genericZvals'
# runID = 'v1_wCPLB_siteSpecificZVals'
runID = 'v1_replicateRobins_siteSpecificZVals'

#create a flag for processing only NGAwest2 and Bradley 2013 models
NGAonly = 0

#define the model you want to compare
#modelToCompare = 'S_22_Crustal'
modelToCompare = 'A_22_Crustal'

#define sites to compare
sites2Plot = ['PVCS', 'PGMS', 'NBSS', 'TFSS', 'PIPS', 'TEPS', 'POTS']

#load a dataframe with geomorphology categories for Wellington
geomorphFilePath = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\residualAnalysis\geomorphology\WellingtonGeomorphology_fromAyushi_Final_ForResidualsPaper_v2.csv'
df_geom = pd.read_csv(geomorphFilePath)

#define directory for saving figures
figDir = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\nonlinearBasinEffects\compareEmpiricalPreds_withRobins'

#create a pdf to save all figures
pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join(figDir, 'compareEmpiricalPredictions_wRobins_%s.pdf' %modelToCompare))

#change plotting settings
plt.style.use('classic')
plt.rcParams["font.family"] = "Times New Roman"

#define directory where GMPE results are saved
#resultsDir = r'C:\Users\cde84\PycharmProjects\EmpiricalGMMs_Files\output\residuals'
resultsDirRoot = r'C:\Users\cde84\PycharmProjects\EmpiricalGMMs_Files\output'
resultsDir = os.path.join(resultsDirRoot, runID, 'models')

#define directory to robins results
resultsDirRobin = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\residualAnalysis\results\GMPE_Residuals_Compiled\Master_September2022_ALL_CompiledOnly_Copy'

#create a dataframe of model IDs and weights used in the NSHM
# gmmWeights_Path = '/home/cde62/NSHM_WelliSiteResponse/residualAnalysis/GMM_NSHM_weights.txt'
gmmWeights_Path = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\residualAnalysis\metadata\GMM_NSHM_weights.txt'
df_weights = pd.read_csv(gmmWeights_Path, delim_whitespace=True)

#define path to stations list file
stationsFilePath = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\nonlinearBasinEffects\sites\stations_ForPaper_wRefStations.txt'
stationsFile = open(stationsFilePath, 'r')

#create a dataframe of the site database based on Jesse's DB v3.2
#welliDB_Path = '/home/cde62/NSHM_WelliSiteResponse/metadata/Liams_GeonetSiteMetadataSummary_v1.4_working.csv'
welliDB_Path = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\gmDatabases\GMDB_3.2\Tables\site_table.csv'
df_siteDB = pd.read_csv(welliDB_Path)

#create list with file headers
#get SA headers from file
#df_forHeaders = pd.read_csv(os.path.join(resultsDir, 'Br_13_rrup_fmin_mw_cmt_HNBN_newrrup_1p0_crustal_mw3p5_rrup300_0p12/Residuals/PJSreStationBiased_sim.csv'))
# df_forHeaders = pd.read_csv(os.path.join(resultsDir, 'A_22_Crustal.csv'), index_col=0)
df_forHeaders = pd.read_csv(os.path.join(resultsDir, 'CY_14_Crustal.csv'), index_col=0)

#SAHeaders = list(df_forHeaders['imName'].values)[1:]
SAHeaders = list(df_forHeaders.columns.values[4:])

#get a list of periods used in the regression from dataframe
periods = [ float(t[4:]) for t in SAHeaders]

print(periods)

#get a list of regions from the geomorph cats file
regionsList = df_geom['Region'].unique()

#create a list of "full GMM ID" (including GMM ID and tec type) from files in resultsDir
filePaths = glob(os.path.join(resultsDir, '*.csv'))
GMM_IDs_List = [path.split('\\')[-1][:-4] for path in filePaths]

#loop through all stations in station file
for i, station in enumerate(stationsFile.read().splitlines()):

    if station not in sites2Plot:
        print('Dont want to plot site: %s' %station)
        continue

    #get the region from geom cat df
    regionForFiles = str(df_geom['Region'][df_geom['stat_id'] == station].values[0]
                         ).replace(' ', '')

    print('station = %s' %station)

    # skip stations that aren't in the site database
    if not station in df_siteDB['sta'].values:
        print('station: %s not in site database' %station)
        continue

    #create list to save all events for a station to compute mean and sigma
    IMratios_AllGMPEs = []

    # get Vs30 and Tsite from site database
    Vs30 = df_siteDB['Vs30'][df_siteDB['sta'] == station].values[0]
    Tsite = df_siteDB['Tsite'][df_siteDB['sta'] == station].values[0]

    # get geomorphic catergory from dataframe
    geomCat = df_geom['Updated Geomorphology'][df_geom['stat_id'] == station].values[0]

    #Loop through different GMPE IDs (full ID including Tec Type)
    for iGMPE, GMPE in enumerate(GMM_IDs_List):

        if GMPE != modelToCompare:
            continue

        # get the GMPE ID without tec type
        gmpe_ID = '_'.join(GMPE.split('_')[:-1])
        gmpe_Tectype = GMPE.split('_')[-1]

        if NGAonly and gmpe_ID not in ['Br_13', 'BSSA_14', 'ASK_14', 'CY_14', 'CB_14']:
            continue

        print('gmpe_ID = %s' %gmpe_ID)
        print('GMPE = %s' % GMPE)

        #load a dataframe with GMPE predictions
        # df_dS2S = pd.read_csv(os.path.join(resultsDir, GMPE, 'Residuals', 'PJSreStationBiased_sim.csv'))
        # df_biasComp = pd.read_csv(os.path.join(resultsDir, GMPE, 'Residuals', 'PJSvarCompsBiased_sim.csv'))
        df_preds = pd.read_csv(os.path.join(resultsDir, '%s.csv' % GMPE), index_col=0)

        #load dataframe of robins predictions
        predDirRobin = os.path.join(resultsDirRobin, GMPE)
        df_predsRobin = pd.read_csv(os.path.join(predDirRobin, 'im_sim.csv'))

        #load station and event files for Robins work to map to residual IDs
        df_stations_Robin = pd.read_csv(os.path.join(predDirRobin, 'stations.csv'))
        df_events_Robin = pd.read_csv(os.path.join(predDirRobin, 'events.csv'))

        # get station ID used in Robins regression
        statID = df_stations_Robin['stat_id'][df_stations_Robin['sta'] == station].values[0]
        refStatID = df_stations_Robin['stat_id'][df_stations_Robin['sta'] == 'POTS'].values[0]

        # skip stations that aren't in the dS2S database
        if not station in df_preds['sta'].values.tolist():
            print('station: %s not in site database' % station)
            continue

        #extract the predictions from all gms for station
        df_preds_Station_allGMs = df_preds[df_preds['sta'] == station]

        IMratios_AllEvents = []

        # loop through events at station
        for evid in df_preds_Station_allGMs['evid'].values:

            print('evid = %s' %evid)

            # get event ID used for Robins regression
            eventID = df_events_Robin['event_id'][df_events_Robin['evid'] == evid].values[0]

            # check if event recorded at POTS
            if evid not in df_preds[df_preds['sta'] == 'POTS']['evid'].values.tolist():
            #if not list(df_preds[np.logical_and(df_preds['sta'] == 'POTS', df_preds['evid'] == evid)].values):
                print('evid %s not recorded at POTS' %evid)
                continue

            # create a figure for plotting
            #fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True)
            fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True)

            # get prediction at POTS
            pred_POTS_evid = df_preds[np.logical_and(df_preds['sta'] == 'POTS', df_preds['evid'] == evid)].values[0][4:].astype(float)
            pred_Station_evid = df_preds[np.logical_and(df_preds['sta'] == station, df_preds['evid'] == evid)].values[0][4:].astype(float)

            IMratio = pred_Station_evid / pred_POTS_evid

            IMratios_AllEvents.append(IMratio)

            # get prediction at POTS
            pred_POTS_evid_Robin = df_predsRobin[np.logical_and(df_predsRobin['stat_id'] == refStatID, df_predsRobin['event_id'] == eventID)].values[0][4:].astype(float)
            pred_Station_evid_Robin = df_predsRobin[np.logical_and(df_predsRobin['stat_id'] == statID, df_predsRobin['event_id'] == eventID)].values[0][4:].astype(float)

            IMratioRobin = pred_Station_evid_Robin / pred_POTS_evid_Robin

            ax.loglog(periods, pred_Station_evid_Robin, color='k', label='Robin')
            ax.loglog(periods, pred_Station_evid, color='r', label='Chris')

            '''
            axs[1].loglog(periods, pred_POTS_evid_Robin, color='k', label='Robin')
            axs[1].loglog(periods, pred_POTS_evid, color='r', label='Chris')

            axs[2].semilogx(periods, IMratio, color='r', label='Chris')
            axs[2].semilogx(periods, IMratioRobin, color='k', label='Robin')
            '''

            ax.set_xlabel('Period (s)')
            ax.set_ylabel('GMM predicted SA (g)')
            ax.set_title('GMM: %s, site: %s,  event: %s' %(GMPE, station, evid))
            ax.legend()
            fig.tight_layout()

            pdf.savefig(fig)

        if not IMratios_AllEvents:
            print('no SA rations for station: %s and GMPE: %s' %(station, GMPE))
            continue

        '''
        #compute the geomean of predictions at station, across all gms
        IMratios_Station_MedianOfEvents = np.exp(np.nanmean(np.log(IMratios_AllEvents), axis=0))

        #append the results to list of all GMPE results for each site
        IMratios_AllGMPEs.append(IMratios_Station_MedianOfEvents)
        '''


    '''
    #if lists are empty, then no events were recorded so skip station
    if not IMratios_AllGMPEs:
        print('No results for station: %s' %station)
        continue

    #compute median and sigmas across all events
    #print('ratios_AllEvents_SA = %s' %ratios_AllEvents_SA)
    if IMratios_AllGMPEs:
        IMratios_AllGMPEs = np.vstack(IMratios_AllGMPEs).astype(float)
        IMratio_Median_AllGMPEs = np.exp(np.nanmean(np.log(IMratios_AllGMPEs), axis=0))
        IMratioSigma_AllGMPEs = np.nanstd(np.log(IMratios_AllGMPEs), axis=0)
    else:
        IMratioMedian_AllGMPEs = np.ones(len(SAHeaders) + 1) * np.nan
        IMratioSigma_AllGMPEs = np.ones(len(SAHeaders) + 1) * np.nan

    #break
    '''

pdf.close()