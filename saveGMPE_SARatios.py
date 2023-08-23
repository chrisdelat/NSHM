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
runID = 'v1_wCPLB_POTS_760mps_siteSpecificZVals'
# runID = 'v1_replicateRobins_siteSpecificZVals'

#create a flag for processing only NGAwest2 and Bradley 2013 models
NGAonly = 1

#load a dataframe with geomorphology categories for Wellington
geomorphFilePath = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\residualAnalysis\geomorphology\WellingtonGeomorphology_fromAyushi_Final_ForResidualsPaper_v2.csv'
df_geom = pd.read_csv(geomorphFilePath)

#define directory for saving files
saveDir = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\nonlinearBasinEffects\IMratios_GMPEs_wrtPOTS'

#define directory where GMPE results are saved
#resultsDir = r'C:\Users\cde84\PycharmProjects\EmpiricalGMMs_Files\output\residuals'
resultsDirRoot = r'C:\Users\cde84\PycharmProjects\EmpiricalGMMs_Files\output'
resultsDir = os.path.join(resultsDirRoot, runID, 'models')

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

#resample period norm period array so that all sites have the same values (for averaging)
Headers = ['stat_id'] + ['gmmID'] + ['Region'] + ['GeoMorph'] + ['TectonicType'] + ['PGA'] + SAHeaders

HeadersMedian = ['stat_id'] + ['Region'] + ['GeoMorph'] + ['TectonicType']+ ['PGA'] + SAHeaders

if NGAonly:
    #create a file to save residuals from all GMPE and mean for each site
    IMratioFile = open(os.path.join(saveDir, 'PGAandSAratios_wrtPOTS_GMPEs_AllSites_AllModels_NGAonly_%s.txt' %runID), 'w')

    #create a file to save median and standard deviation of residuals
    IMratioMedianFile = open(os.path.join(saveDir, 'PGAandSAratios_wrtPOTS_GMPEs_Medians_AllSites_NGAonly_%s.txt' %runID), 'w')
    IMratioSigmaFile = open(os.path.join(saveDir, 'PGAandSAratios_wrtPOTS_GMPEs_Sigmas_AllSites_NGAonly_%s.txt' %runID), 'w')
else:
    # create a file to save residuals from all GMPE and mean for each site
    IMratioFile = open(os.path.join(saveDir, 'PGAandSAratios_wrtPOTS_GMPEs_AllSites_AllModels_%s.txt' %runID), 'w')

    # create a file to save median and standard deviation of residuals
    IMratioMedianFile = open(os.path.join(saveDir, 'PGAandSAratios_wrtPOTS_GMPEs_Medians_AllSites_%s.txt' %runID), 'w')
    IMratioSigmaFile = open(os.path.join(saveDir, 'PGAandSAratios_wrtPOTS_GMPEs_Sigmas_AllSites_%s.txt' %runID), 'w')

#write headers on files
for header in Headers[0:-1]:
    IMratioFile.write('%s ' %header)
IMratioFile.write('%s\n' %Headers[-1])

for header in HeadersMedian[0:-1]:
    IMratioMedianFile.write('%s ' %header)
    IMratioSigmaFile.write('%s ' %header)

IMratioMedianFile.write('%s\n' %HeadersMedian[-1])
IMratioSigmaFile.write('%s\n' %HeadersMedian[-1])

#get a list of regions from the geomorph cats file
regionsList = df_geom['Region'].unique()

#create a list of "full GMM ID" (including GMM ID and tec type) from files in resultsDir
filePaths = glob(os.path.join(resultsDir, '*.csv'))
GMM_IDs_List = [path.split('\\')[-1][:-4] for path in filePaths]

#loop through all stations in station file
for i, station in enumerate(stationsFile.read().splitlines()):

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

            # check if event recorded at POTS
            if evid not in df_preds[df_preds['sta'] == 'POTS']['evid'].values.tolist():
            #if not list(df_preds[np.logical_and(df_preds['sta'] == 'POTS', df_preds['evid'] == evid)].values):
                print('evid %s not recorded at POTS' %evid)
                continue

            # get prediction at POTS
            pred_POTS_evid = df_preds[np.logical_and(df_preds['sta'] == 'POTS', df_preds['evid'] == evid)].values[0][3:].astype(float)
            pred_Station_evid = df_preds[np.logical_and(df_preds['sta'] == station, df_preds['evid'] == evid)].values[0][3:].astype(float)

            IMratio = pred_Station_evid / pred_POTS_evid

            IMratios_AllEvents.append(IMratio)

        if not IMratios_AllEvents:
            print('no SA rations for station: %s and GMPE: %s' %(station, GMPE))
            continue

        #compute the geomean of predictions at station, across all gms
        IMratios_Station_MedianOfEvents = np.exp(np.nanmean(np.log(IMratios_AllEvents), axis=0))

        #append the results to list of all GMPE results for each site
        IMratios_AllGMPEs.append(IMratios_Station_MedianOfEvents)

        #write event IM and residuals to files
        IMratioFile.write('%s %s %s %s %s ' %(station, gmpe_ID, regionForFiles, geomCat, gmpe_Tectype) )

        #loop through SA lists to write to file
        for IMratio in IMratios_Station_MedianOfEvents[:-1]:
            IMratioFile.write('%s ' %IMratio)
        IMratioFile.write('%s\n' %IMratios_Station_MedianOfEvents[-1])

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

    '''
    #append median and sigmas of each model to list
    ratios_AllModelMedians_SA.append(Ratio_Median_SA)
    ratios_AllModelSigmas_SA.append(Ratio_eventSigma_SA)
    ratios_AllModelMedians_PGA.append(Ratio_Median_PGA)
    ratios_AllModelSigmas_PGA.append(Ratio_eventSigma_PGA)
    '''

    #write median and sigma to file for All GMPEs
    IMratioMedianFile.write('%s %s %s All ' %(station, regionForFiles, geomCat) )
    IMratioSigmaFile.write('%s %s %s All ' %(station, regionForFiles, geomCat) )

    #loop through SA lists to write to file
    for (IMratioMed, IMratioSig) in zip(IMratio_Median_AllGMPEs[:-1], IMratioSigma_AllGMPEs[:-1]):
        IMratioMedianFile.write('%s ' %IMratioMed)
        IMratioSigmaFile.write('%s ' %IMratioSig)

    IMratioMedianFile.write('%s\n' %IMratio_Median_AllGMPEs[-1])
    IMratioSigmaFile.write('%s\n' %IMratioSigma_AllGMPEs[-1])



    ##Region mean
    #skip stations that don't exist in Database
    if station == 'WSTS':
    #if station == 'WSTS' or station == 'NBSS':
        continue
    #skip stations that are likely not basin sites
    if station == 'WTES' or station == 'WCFS':
        continue

    #break

IMratioFile.close()
IMratioMedianFile.close()
IMratioSigmaFile.close()
