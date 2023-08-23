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
from scipy.stats.mstats import gmean
from math import ceil, log
import eqsig

from vel2acc import acc2vel
from readASCII import read_ascii
from computeFourier import computeFourier


#Create a list with paths to all different station ID files
stationsDir = '/home/cde62/NSHM_WelliSiteResponse/stations'

#create a counting variable to only compute smoothing matrix once
count = 0

#define periods at which SA was computed
periods = np.power(10,np.linspace(-2, 1, 100))

#define directory for saving files
saveDir = '/home/cde62/NSHM_WelliSiteResponse/results/GMPE_Ratios'
metadataDir = '/home/cde62/NSHM_WelliSiteResponse/metadata'

#define directory where GMPE results are saved
resultsDir = '/nesi/nobackup/uc02594/OpenSees/NSHM_WelliBasin/GMPE_Predictions_Compiled'

#create a dataframe with observed SA values of full database
path_SAdatabase = '/scale_wlg_nobackup/filesets/nobackup/nesi00213/RunFolder/Validation/WelliDatabase_2022-01-24/SA_obs.csv'
df_SA_database = pd.read_csv(path_SAdatabase)

#get a list of all unique events in the database
EQList_WelliDatabase = df_SA_database['event_id'].unique()


#create list with file headers
#get SA headers from file
df_forHeaders = pd.read_csv(os.path.join(resultsDir, 'ASK14/im_sim.csv'))
SAHeaders = list(df_forHeaders.keys()[4:])

Headers = ['stat_id'] + ['Region'] + ['SubRegion'] + ['ref_stat_id'] + ['event_id'] + ['PGARock_g'] + ['PGASoil_g'] + ['PGA'] + SAHeaders

HeadersMedian = ['stat_id'] + ['Region'] + ['SubRegion'] + ['ref_stat_id'] + ['event_id'] + ['PGA'] + SAHeaders

'''
Headers = ['stat_id'] + ['Region'] + ['SubRegion'] + ['ref_stat_id'] + ['event_id'] + ['PGARock_g'] + ['PGASoil_g'] + ['PGA'] +\
          ['SA_' + str(x).replace('.', 'p') for x in np.round(periods, 3)]

HeadersMedian = ['stat_id'] + ['Region'] + ['SubRegion'] + ['ref_stat_id'] + ['event_id'] + ['PGA'] +\
          ['SA_' + str(x).replace('.', 'p') for x in np.round(periods, 3)]
'''

#Loop through different GMPE folders
for GMPE in os.listdir(resultsDir):
    
    #load dataframes of event and station IDs
    df_stations = pd.read_csv(os.path.join(resultsDir, GMPE, 'stations.csv'))
    df_events = pd.read_csv(os.path.join(resultsDir, GMPE, 'events.csv'))
    
    #load a dataframe with predicted IMs
    df_IMs = pd.read_csv(os.path.join(resultsDir, GMPE, 'im_sim.csv'))

    #create a file to save AF ratios for all events
    IMRatiosFile = open(os.path.join(saveDir, 'RatiosGMPE_SA_AFs_AllSites_AllEvents_%s.txt' %GMPE ), 'w')

    #create a file to save median and standard deviation of AF ratios for all sites
    IMRatiosMedianFile = open(os.path.join(saveDir, 'RatiosGMPE_Median_SA_AFs_AllSites_%s.txt' % GMPE), 'w')

    #create a file to save median and standard deviation of AF ratios for all sites
    IMRatiosSigmaFile = open(os.path.join(saveDir, 'RatiosGMPE_Sigma_SA_AFs_AllSites_%s.txt' %GMPE), 'w')
    
    #write headers on files
    for header in Headers[0:-1]:
        IMRatiosFile.write('%s ' %header)
    IMRatiosFile.write('%s\n' %Headers[-1])

    for header in HeadersMedian[0:-1]:
        IMRatiosMedianFile.write('%s ' %header)
        IMRatiosSigmaFile.write('%s ' %header)
    IMRatiosMedianFile.write('%s\n' %HeadersMedian[-1])
    IMRatiosSigmaFile.write('%s\n' %HeadersMedian[-1])

    #loop through different station files (i.e., regions)
    for stationFileName in os.listdir(stationsDir):
        
        stationFileSplit = stationFileName.split('_')

        stationRegion = stationFileSplit[1]
        
        print(stationRegion)

        stationSubRegion = stationFileSplit[-1].split('.')[0]
        print('sub region is %s' %stationSubRegion)
        
        stationFilePath = os.path.join(stationsDir, stationFileName)
        print('sub region  station file path is: %s' %stationFilePath)
        stationsFile = open(stationFilePath, 'r')

        if stationRegion == 'Welli':
            refStation = 'POTS'
        elif stationRegion == 'Hutt':
            refStation = 'POTS'
        else:
            if stationSubRegion == 'ReferenceRock':
                refStation = 'POTS'
            elif  stationSubRegion == 'Porirua':
                refStation = 'PWES'
            elif stationSubRegion == 'Wainuiomata':
                refStation = 'POTS'
            elif stationSubRegion == 'LakeWairarapaBasin':
                refStation = 'KIRS'
        
        #get station ID used in regression for reference station
        refStatID = df_stations['stat_id'][df_stations['sta'] == refStation].values[0]
        print('refStatID = %s' %refStatID)

        #loop through all stations in station file
        for i, station in enumerate(stationsFile.read().splitlines()):
            
            print(station)
            
            #check if station exists in GMPE database
            if not df_stations['stat_id'][df_stations['sta'] == station].values:
                continue

            #get station ID used in regression
            statID = df_stations['stat_id'][df_stations['sta'] == station].values[0]
            print('statID = %s' %statID)

            #create a list of all events for station, to check for repeats between databases
            EQs_4Station = []
            
            #create list to save all events for a station to compute mean and sigma
            ratios_AllEvents_SA = []
            ratios_AllEvents_PGA = []
            ratios_AllEvents_PGV = []
            

            #Loop through all event folder in given directory
            #for eqDatabaseDir in eqDatabasePaths:
            #    print('starting database: %s' %eqDatabaseDir)
            #for EQ in os.listdir(eqDatabaseDir):
            for EQ in EQList_WelliDatabase:
                
                print(EQ)
                

                #print("list(df_SA_database[np.logical_and(df_SA_database['Station']==station, df_SA_database['Event'] == EQ)].values) = %s"\
                #      %list(df_SA_database[np.logical_and(df_SA_database['Station']==station, df_SA_database['Event'] == EQ)].values))
                #print("list(df_SA_database[np.logical_and(df_SA_database['Station']==refStation, df_SA_database['Event'] == EQ)].values) = %s"\
                #      %list(df_SA_database[np.logical_and(df_SA_database['Station']==refStation, df_SA_database['Event'] == EQ)].values))
                
                if not list(df_SA_database[np.logical_and(df_SA_database['stat_id']==station, df_SA_database['event_id'] == EQ)].values)\
                or not list(df_SA_database[np.logical_and(df_SA_database['stat_id']==refStation, df_SA_database['event_id'] == EQ)].values):
                    print('event: %s not recorded at ref station: %s or station: %s' %(EQ, refStation, station))
                    continue
                
                print('event: %s was recorded at ref station: %s and station: %s' %(EQ, refStation, station))
                
                #check if event has already been processed for this station (i.e., if there are repeats between databases)
                if EQ in EQs_4Station:
                    continue
                else:
                    EQs_4Station.append(EQ)
                
                #check if event exists in GMPE database (which only includes crustal events)
                if not df_events['event_id'][df_events['evid'] == EQ].values:
                    continue
                
                #get event ID used for regression
                eventID = df_events['event_id'][df_events['evid'] == EQ].values[0]
                print('eventID = %s' %eventID)
                
                
                #print("df_SA_database[np.logical_and(df_SA_database['Station']==station, df_SA_database['Event'] == EQ)].values[0] = %s"\
                #       %(df_SA_database[np.logical_and(df_SA_database['Station']==station, df_SA_database['Event'] == EQ)].values[0]))
                
                #print("df_SA_database[np.logical_and(df_SA_database['Station']==station, df_SA_database['Event'] == EQ)].values[0][4:] = %s"\
                #       %(df_SA_database[np.logical_and(df_SA_database['Station']==station, df_SA_database['Event'] == EQ)].values[0][4:]))

                ratio_SA = np.array(df_IMs[np.logical_and(df_IMs['stat_id']==statID, df_IMs['event_id'] == eventID)].values[0][4:]\
                           / df_IMs[np.logical_and(df_IMs['stat_id']==refStatID, df_IMs['event_id'] == eventID)].values[0][4:])
                ratio_PGA = np.array(df_IMs[np.logical_and(df_IMs['stat_id']==statID, df_IMs['event_id'] == eventID)].values[0][3]\
                           / df_IMs[np.logical_and(df_IMs['stat_id']==refStatID, df_IMs['event_id'] == eventID)].values[0][3])
                
                PGARock = df_IMs[np.logical_and(df_IMs['stat_id']==refStatID, df_IMs['event_id'] == eventID)].values[0][3]
                PGASoil = df_IMs[np.logical_and(df_IMs['stat_id']==statID, df_IMs['event_id'] == eventID)].values[0][3]
                
                '''
                #extract fmins from database
                fmin_Y = df_SA_database[np.logical_and(df_SA_database['stat_id']==station, df_SA_database['event_id'] == EQ)].values[0][10]
                fmin_X = df_SA_database[np.logical_and(df_SA_database['stat_id']==station, df_SA_database['event_id'] == EQ)].values[0][11]
                fmin_Y_ref = df_SA_database[np.logical_and(df_SA_database['stat_id']==refStation, df_SA_database['event_id'] == EQ)].values[0][10]
                fmin_X_ref = df_SA_database[np.logical_and(df_SA_database['stat_id']==refStation, df_SA_database['event_id'] == EQ)].values[0][11]
                
                print('f_mins site = %s and %s' %(fmin_Y, fmin_X))
                print('f_mins ref = %s and %s' %(fmin_Y_ref, fmin_X_ref))
                
                #identify the highest fmin of both components of site and reference motions
                f_min_4Ratio = np.max([fmin_Y, fmin_X, fmin_Y_ref, fmin_X_ref])
                
                #calculate Tmax from fmin
                Tmax = 1.0 / f_min_4Ratio

                #truncate ratio arrays for periods longer than Tmax
                ratio_SA[np.where(periods > Tmax)] = np.nan
                '''

                #append ratios for event to list for all events
                ratios_AllEvents_SA.append(ratio_SA)
                ratios_AllEvents_PGA.append(ratio_PGA)

                print(EQ)
                
                #write event IM and residuals to files
                IMRatiosFile.write('%s %s %s %s %s %s %s %s ' %(station, stationRegion, stationSubRegion, refStation, EQ, PGARock, PGASoil, ratio_PGA) )
                
                #loop through SA lists to write to file
                for iSA_Ratio, SA_Ratio in enumerate(ratio_SA[:-1]):
                    IMRatiosFile.write('%s ' %SA_Ratio)
                
                IMRatiosFile.write('%s\n' %ratio_SA[-1])
            
            #if arrays are empty, then no events were recorded so skip station
            if not ratios_AllEvents_SA:
                continue

            #compute median and sigmas across all events
            #print('ratios_AllEvents_SA = %s' %ratios_AllEvents_SA)
            ratios_AllEvents_SA = np.vstack(ratios_AllEvents_SA).astype(np.float)
            ratios_AllEvents_PGA = np.vstack(ratios_AllEvents_PGA).astype(np.float)
            
            Ratio_Median_SA = np.exp(np.nanmean(np.log(ratios_AllEvents_SA), axis=0))
            #Ratio_Median_SA = gmean(ratios_AllEvents_SA, axis=0)
            Ratio_Median_PGA = gmean(ratios_AllEvents_PGA)[0]

            Ratio_eventSigma_SA = np.nanstd(np.log(ratios_AllEvents_SA), axis=0)
            #Ratio_eventSigma_SA = np.std(np.log(ratios_AllEvents_SA), axis=0)
            Ratio_eventSigma_PGA = np.std(np.log(ratios_AllEvents_PGA))
            
            '''
            #append median and sigmas of each model to list
            ratios_AllModelMedians_SA.append(Ratio_Median_SA)
            ratios_AllModelSigmas_SA.append(Ratio_eventSigma_SA)
            ratios_AllModelMedians_PGA.append(Ratio_Median_PGA)
            ratios_AllModelSigmas_PGA.append(Ratio_eventSigma_PGA)
            '''

            IMRatiosMedianFile.write('%s %s %s %s %s %s ' %(station, stationRegion, stationSubRegion, refStation, 'Median', Ratio_Median_PGA) )
            IMRatiosSigmaFile.write('%s %s %s %s %s %s ' %(station, stationRegion, stationSubRegion, refStation, 'btwnEvent_Sigma', Ratio_eventSigma_PGA) )
            
            #loop through SA lists to write to file
            for iSA_MedianRatio, SA_MedianRatio in enumerate(Ratio_Median_SA[:-1]):
                IMRatiosMedianFile.write('%s ' %SA_MedianRatio)
                IMRatiosSigmaFile.write('%s ' %Ratio_eventSigma_SA[iSA_MedianRatio])
            
            IMRatiosMedianFile.write('%s\n' %Ratio_Median_SA[-1])
            IMRatiosSigmaFile.write('%s\n' %Ratio_eventSigma_SA[-1])

        #break
    IMRatiosFile.close()
    IMRatiosMedianFile.close()
    IMRatiosSigmaFile.close()


