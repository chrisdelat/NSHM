# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 16:28:08 2018

@author: cde62

Computes the median and std dev of eHVSR from Chuanbin's database

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
saveDir = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\gmDatabases\Chaunbins_FAS_01.03.2023\EQ_HVSR\MedianAndStd'
os.makedirs(saveDir, exist_ok=True)

metadataDir = '/home/cde62/NSHM_WelliSiteResponse/metadata'

# create a dataframe with FAS values of full database
#path_FASdatabase = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\gmDatabases\Chaunbins_FAS_01.03.2023\0_EAS_0p5_25p0_crustal_b40.xlsx'
path_eHVSRdatabase = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\gmDatabases\Chaunbins_FAS_01.03.2023\EQ_HVSR\0_HVSR_0p1_25p0_crustal_b40.xlsx'
df_eHVSR_database = pd.read_excel(path_eHVSRdatabase)

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

#create a file to save median and standard deviation of AF ratios for all sites
IMRatiosMedianFile = open(os.path.join(saveDir, 'Obs_eHVSR_FAS_Median_AllSites_%s.csv' %(analID)), 'w')

#create a file to save median and standard deviation of AF ratios for all sites
IMRatiosSigmaFile = open(os.path.join(saveDir, 'Obs_eHVSR_FAS_Sigma_AllSites_%s.csv' %(analID)), 'w')

#create list with file headers
#HeadersMedian = ['stat_id'] + ['Region'] + ['numEvents'] + FAS_headers
HeadersMedian = ['stat_id'] + ['numEvents'] + FAS_headers

#write headers on files
for header in HeadersMedian[0:-1]:
    IMRatiosMedianFile.write('%s,' %header)
    IMRatiosSigmaFile.write('%s,' %header)
IMRatiosMedianFile.write('%s\n' %HeadersMedian[-1])
IMRatiosSigmaFile.write('%s\n' %HeadersMedian[-1])

# create a list of all sites
stationList = list(df_eHVSR_database['sta'].unique())
print('stationList = %s' %stationList)


#loop through different station files (i.e., regions)
for i, station in enumerate(stationList):

    #stationRegion = df_geom['Region'][df_geom['stat_id'] == station].values[0]
        
    #create a list of all events for station, to check for repeats between databases
    EQs_4Station = []

    #create list to save all events for a station to compute mean and sigma
    ratios_AllEvents_FAS = []
    ratios_AllEvents_PGA = []
    ratios_AllEvents_PGV = []

    print(station)

    #create a sub df with only data for station
    sub_df_eHVSR_database = df_eHVSR_database[df_eHVSR_database['sta'].isin([station]) & df_eHVSR_database['chan'].isin(['HH', 'HN', 'BN'])]



    # turn the sub df into an array
    eHVSR_site_allEvents = np.array(sub_df_eHVSR_database.iloc[:, 33:])

    if len(eHVSR_site_allEvents) == 0:
        continue

    #print(eHVSR_site_allEvents)
    #print(len(eHVSR_site_allEvents))

    # get the number of events for each site
    numEvents = len(eHVSR_site_allEvents)
    # compute median and sigmas across all events
    # print('ratios_AllEvents_FAS = %s' %ratios_AllEvents_FAS)
    #ratios_AllEvents_FAS = np.vstack(ratios_AllEvents_FAS).astype(np.float)

    Ratio_Median_FAS = np.exp(np.nanmean(np.log(eHVSR_site_allEvents), axis=0))
    # Ratio_Median_FAS = gmean(ratios_AllEvents_FAS, axis=0)

    Ratio_eventSigma_FAS = np.nanstd(np.log(eHVSR_site_allEvents), axis=0)
    # Ratio_eventSigma_FAS = np.std(np.log(ratios_AllEvents_FAS), axis=0)

    #IMRatiosMedianFile.write('%s,%s,%s,' % (station, stationRegion, numEvents))
    #IMRatiosSigmaFile.write('%s,%s,%s,' % (station, stationRegion, numEvents))

    IMRatiosMedianFile.write('%s,%s,' % (station, numEvents))
    IMRatiosSigmaFile.write('%s,%s,' % (station, numEvents))

    # loop through FAS lists to write to file
    for iFAS_MedianRatio, FAS_MedianRatio in enumerate(Ratio_Median_FAS[:-1]):
        IMRatiosMedianFile.write('%s,' % FAS_MedianRatio)
        IMRatiosSigmaFile.write('%s,' % Ratio_eventSigma_FAS[iFAS_MedianRatio])

    IMRatiosMedianFile.write('%s\n' % Ratio_Median_FAS[-1])
    IMRatiosSigmaFile.write('%s\n' % Ratio_eventSigma_FAS[-1])

    # break
IMRatiosMedianFile.close()
IMRatiosSigmaFile.close()
