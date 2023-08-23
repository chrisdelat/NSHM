# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 16:28:08 2018

@author: cde62

Computes response spectral ratios and smoothed fourier ratios for strong motion station in station.txt relative to reference station POTS

"""
import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import matplotlib.lines as mlines
from obspy.signal import konnoohmachismoothing
import os.path
import pandas as pd
from scipy.stats.mstats import gmean
from math import ceil, log
import eqsig



#Create a list with paths to all different station ID files
stationsDir = '/home/cde62/NSHM_WelliSiteResponse/stations'

# create a dataframe for geomorphic categories of all sites (to get a site list)
geomorphFilePath = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\Site Databases\WellingtonGeomorphology_fromAyushi_Final_ForResidualsPaper.csv'
df_geom = pd.read_csv(geomorphFilePath)

# create a list of all Wellington sites
stationList = list(df_geom['stat_id'].values)
print('stationList = %s' %stationList)

# define the reference station for SSR
refStation = 'POTS'

#create a counting variable to only compute smoothing matrix once
count = 0

#define the factor appled to Tmax
TmaxFact = 1.0

#define smoothing flag and b value
smoothFlag = 0
bVal = 80

#define periods at which FAS was computed
#periods = np.power(10,np.linspace(-2, 1, 100))

#define directory for saving files
saveDir = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\SSRh\FAS_EQ_SSR_ChuanbinsDB'
os.makedirs(saveDir, exist_ok=True)

metadataDir = '/home/cde62/NSHM_WelliSiteResponse/metadata'

# create a dataframe with FAS values of full database
#path_FASdatabase = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\gmDatabases\Chaunbins_FAS_01.03.2023\0_EAS_0p5_25p0_crustal_b40.xlsx'
path_FASdatabase = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\gmDatabases\Chaunbins_FAS_01.03.2023\0_EAS_0p1_25p0_crustal_b40.xlsx'
df_FAS_database = pd.read_excel(path_FASdatabase)

# load array of frequencies for EAS
freqFilePath = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\gmDatabases\Chaunbins_FAS_01.03.2023\freq_0p1_25p0.txt'
FAS_freqs = np.loadtxt(freqFilePath)
FAS_headers = [str(freq) + '_Hz' for freq in FAS_freqs]
                
#create tags for mean TF files to identify smoothing based on inputs
if smoothFlag:
    analID = 'smoothFAS_b%s' %bVal
else:
    analID = 'unsmooth'

#smooth FAS before computing TF if smooth=1
if smoothFlag:
    print('computing smoothing matrix')
    smooth_matrix = konnoohmachismoothing.calculate_smoothing_matrix(FAS_freqs, bandwidth=bVal, normalize=True)
    #smooth_matrix = np.ones(len(ft_freq))
    print('smoothing matrix complete')

#create a file to save AF ratios for all events
IMRatiosFile = open(os.path.join(saveDir, 'Obs_FAS_SSR_wrt%s_%s_AllSites_AllEvents.csv' %(refStation, analID)), 'w')

#create a file to save median and standard deviation of AF ratios for all sites
IMRatiosMedianFile = open(os.path.join(saveDir, 'Obs_FAS_SSR_wrt%s_%s_Median_AllSites.csv' %(refStation, analID)), 'w')

#create a file to save median and standard deviation of AF ratios for all sites
IMRatiosSigmaFile = open(os.path.join(saveDir, 'Obs_FAS_SSR_wrt%s_%s_Sigma_AllSites.csv' %(refStation, analID)), 'w')

#create list with file headers
Headers = ['stat_id'] + ['Region'] + ['ref_stat_id'] + ['event_id'] + ['PGARock_g'] + ['PGASoil_g'] + FAS_headers

HeadersMedian = ['stat_id'] + ['Region'] + ['ref_stat_id'] + ['event_id'] + ['numEvents'] + FAS_headers

#write headers on files
for header in Headers[0:-1]:
    IMRatiosFile.write('%s,' %header)
IMRatiosFile.write('%s\n' %Headers[-1])

for header in HeadersMedian[0:-1]:
    IMRatiosMedianFile.write('%s,' %header)
    IMRatiosSigmaFile.write('%s,' %header)
IMRatiosMedianFile.write('%s\n' %HeadersMedian[-1])
IMRatiosSigmaFile.write('%s\n' %HeadersMedian[-1])

#loop through different station files (i.e., regions)
for i, station in enumerate(stationList):
    stationRegion = df_geom['Region'][df_geom['stat_id'] == station].values[0]
        
    #create a list of all events for station, to check for repeats between databases
    EQs_4Station = []

    #create list to save all events for a station to compute mean and sigma
    ratios_AllEvents_FAS = []
    ratios_AllEvents_PGA = []
    ratios_AllEvents_PGV = []

    print(station)

    #create a sub df with only data for station
    sub_df_FAS_database = df_FAS_database[df_FAS_database['sta'].isin([station]) & df_FAS_database['chan'].isin(['HH', 'HN', 'BN'])]

    # Loop through all events recorded for station
    for df_index_EQ, df_row_EQ in  sub_df_FAS_database.iterrows():

        EQ = df_row_EQ['evid']

        print(EQ)


        #print("list(df_FAS_database[np.logical_and(df_FAS_database['Station']==station, df_FAS_database['Event'] == EQ)].values) = %s"\
        #      %list(df_FAS_database[np.logical_and(df_FAS_database['Station']==station, df_FAS_database['Event'] == EQ)].values))
        #print("list(df_FAS_database[np.logical_and(df_FAS_database['Station']==refStation, df_FAS_database['Event'] == EQ)].values) = %s"\
        #      %list(df_FAS_database[np.logical_and(df_FAS_database['Station']==refStation, df_FAS_database['Event'] == EQ)].values))

        # check if EQ was recorded at reference station
        if not list(df_FAS_database[np.logical_and(df_FAS_database['sta']==refStation, df_FAS_database['evid'] == EQ)].values):
            print('event: %s not recorded at ref station: %s' %(EQ, refStation))
            continue

        print('event: %s was recorded at ref station: %s and station: %s' % (EQ, refStation, station))

        #row_site = df_FAS_database[np.logical_and(df_FAS_database['sta']==station, df_FAS_database['evid'] == EQ)].values
        #row_ref = df_FAS_database[np.logical_and(df_FAS_database['sta']==refStation, df_FAS_database['evid'] == EQ)].values

        #len_row_site = len(row_site)
        #len_row_ref = len(row_ref)

        FAS_site = df_FAS_database[np.logical_and(df_FAS_database['sta']==station, df_FAS_database['evid'] == EQ)].values[0][33:]
        FAS_ref = df_FAS_database[np.logical_and(df_FAS_database['sta']==refStation, df_FAS_database['evid'] == EQ)].values[0][33:]

        ratio_FAS = np.array(FAS_site / FAS_ref).astype(float)

        #ratio_FAS = np.array(df_FAS_database[np.logical_and(df_FAS_database['sta']==station, df_FAS_database['evid'] == EQ)].values[0][33:]\
        #           / df_FAS_database[np.logical_and(df_FAS_database['sta']==refStation, df_FAS_database['evid'] == EQ)].values[0][33:]).astype(float)
        #ratio_PGA = np.array(df_FAS_database[np.logical_and(df_FAS_database['stat_id']==station, df_FAS_database['event_id'] == EQ)].values[0][13]\
        #           / df_FAS_database[np.logical_and(df_FAS_database['stat_id']==refStation, df_FAS_database['event_id'] == EQ)].values[0][13])

        #print('ratio_FAS = %s' %ratio_FAS)

        # smooth the FAS SSR
        if smoothFlag:
            ratio_FAS = np.exp(np.dot(ratio_FAS, smooth_matrix))

        PGARock = df_FAS_database[np.logical_and(df_FAS_database['sta']==refStation, df_FAS_database['evid'] == EQ)].values[0][6]
        PGASoil = df_FAS_database[np.logical_and(df_FAS_database['sta']==station, df_FAS_database['evid'] == EQ)].values[0][6]

        #append ratios for event to list for all events
        ratios_AllEvents_FAS.append(ratio_FAS)

        print(EQ)

        #write event IM and residuals to files
        IMRatiosFile.write('%s,%s,%s,%s,%s,%s,' %(station, stationRegion, refStation, EQ, PGARock, PGASoil) )

        #loop through FAS lists to write to file
        for iFAS_Ratio, FAS_Ratio in enumerate(ratio_FAS[:-1]):
            IMRatiosFile.write('%s,' %FAS_Ratio)

        IMRatiosFile.write('%s\n' %ratio_FAS[-1])


    #if arrays are empty, then no events were recorded so skip station
    #print(ratios_AllEvents_FAS)

    if not ratios_AllEvents_FAS:
        continue

    # get the number of events for each site
    numEvents = len(ratios_AllEvents_FAS)
    #compute median and sigmas across all events
    #print('ratios_AllEvents_FAS = %s' %ratios_AllEvents_FAS)
    ratios_AllEvents_FAS = np.vstack(ratios_AllEvents_FAS).astype(np.float)

    Ratio_Median_FAS = np.exp(np.nanmean(np.log(ratios_AllEvents_FAS), axis=0))
    #Ratio_Median_FAS = gmean(ratios_AllEvents_FAS, axis=0)

    Ratio_eventSigma_FAS = np.nanstd(np.log(ratios_AllEvents_FAS), axis=0)
    #Ratio_eventSigma_FAS = np.std(np.log(ratios_AllEvents_FAS), axis=0)

    IMRatiosMedianFile.write('%s,%s,%s,%s,%s,' %(station, stationRegion, refStation, numEvents, 'Median') )
    IMRatiosSigmaFile.write('%s,%s,%s,%s,%s,' %(station, stationRegion, refStation, numEvents, 'btwnEvent_Sigma') )

    #loop through FAS lists to write to file
    for iFAS_MedianRatio, FAS_MedianRatio in enumerate(Ratio_Median_FAS[:-1]):
        IMRatiosMedianFile.write('%s,' %FAS_MedianRatio)
        IMRatiosSigmaFile.write('%s,' %Ratio_eventSigma_FAS[iFAS_MedianRatio])

    IMRatiosMedianFile.write('%s\n' %Ratio_Median_FAS[-1])
    IMRatiosSigmaFile.write('%s\n' %Ratio_eventSigma_FAS[-1])

    # break
IMRatiosFile.close()
IMRatiosMedianFile.close()
IMRatiosSigmaFile.close()


