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

from vel2acc import acc2vel
from readASCII import read_ascii
from computeFourier import computeFourier


#Create a list with paths to all different station ID files
stationsDir = '/home/cde62/NSHM_WelliSiteResponse/stations'

#define periods at which SA was computed
periods = np.power(10,np.linspace(-2, 1, 100))

#define directory for saving files
saveDir = '/home/cde62/NSHM_WelliSiteResponse/results/GMPE_Ratios'

#define directory where GMPE results are saved
ratiosDir = '/home/cde62/NSHM_WelliSiteResponse/results/GMPE_Ratios'
resultsDir = '/nesi/nobackup/uc02594/OpenSees/NSHM_WelliBasin/GMPE_Predictions_Compiled'


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
    
#create a files to save median and standard deviation of AF ratios across all models
IMRatiosMedianFile_allModels = open(os.path.join(saveDir, 'RatiosGMPE_MedianOfAllModels_SA_AFs.txt'), 'w')
#IMRatiosSigmaFile_allModels = open(os.path.join(saveDir, 'RatiosGMPE_Sigma_SA_AFs_AllSites_%s.txt' %GMPE), 'w')
    

for header in HeadersMedian[0:-1]:
    IMRatiosMedianFile_allModels.write('%s ' %header)
IMRatiosMedianFile_allModels.write('%s\n' %HeadersMedian[-1])
    

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

    #loop through all stations in station file
    for i, station in enumerate(stationsFile.read().splitlines()):
        #create lists to save all medians and sigmas from all GMPEs
        ratios_AllModelMedians_SA = []
        ratios_AllModelMedians_PGA = []
        #ratios_AllModelSigmas_SA = []
        #ratios_AllModelSigmas_PGA = []
        
        print(station)
        
        #Loop through all GMPE median files
        for GMPEfilePath in glob(os.path.join(ratiosDir, '*_Median_*')):
            print('GMPEfilePath = %s' %GMPEfilePath)
            
            #load dataframe with median ratios for all sites (one GMPE)
            df_Medians_1GMPE = pd.read_csv(GMPEfilePath, delimiter=' ')
            
            #check if station exists in ratio database
            if not list(df_Medians_1GMPE['stat_id'][df_Medians_1GMPE['stat_id'] == station].values):
                continue
            
            MedianRatio_SA = df_Medians_1GMPE[df_Medians_1GMPE['stat_id'] == station].values[0][6:]
            MedianRatio_PGA = df_Medians_1GMPE[df_Medians_1GMPE['stat_id'] == station].values[0][5]
            
            #append ratios for event to list for all events
            ratios_AllModelMedians_SA.append(MedianRatio_SA)
            ratios_AllModelMedians_PGA.append(MedianRatio_PGA)

        #skip station if there was no values found for it
        if not ratios_AllModelMedians_SA:
            continue

        #compute median and sigmas across all events
        print('ratios_AllModelMedians_SA = %s' %ratios_AllModelMedians_SA)
        print('ratios_AllModelMedians_PGA = %s' %ratios_AllModelMedians_PGA)
        ratios_AllModelMedians_SA = np.array(ratios_AllModelMedians_SA).astype(np.float)
        ratios_AllModelMedians_PGA = np.array(ratios_AllModelMedians_PGA)
        
        Ratio_Median_SA = np.exp(np.nanmean(np.log(ratios_AllModelMedians_SA), axis=0))
        Ratio_Median_PGA = gmean(ratios_AllModelMedians_PGA)

        IMRatiosMedianFile_allModels.write('%s %s %s %s %s %s ' %(station, stationRegion, stationSubRegion, refStation, 'Median', Ratio_Median_PGA) )
        
        #loop through SA lists to write to file
        for iSA_MedianRatio, SA_MedianRatio in enumerate(Ratio_Median_SA[:-1]):
            IMRatiosMedianFile_allModels.write('%s ' %SA_MedianRatio)
        
        IMRatiosMedianFile_allModels.write('%s\n' %Ratio_Median_SA[-1])

    #break
IMRatiosMedianFile_allModels.close()


