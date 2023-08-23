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

#create a flag for processing only NGAwest2 and Bradley 2013 models
NGAonly = 0

#load a dataframe with geomorphology categories for Wellington
geomorphFilePath = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\residualAnalysis\geomorphology\WellingtonGeomorphology_fromAyushi_Final_ForResidualsPaper_v2.csv'
df_geom = pd.read_csv(geomorphFilePath)

#specify directory for saving figs
figDir = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\JournalPaper_WelliSiteResp\nonlinearBasinEffects\R1_Submission2\figs'

#define directory where GMPE results are saved
#resultsDir = r'C:\Users\cde84\PycharmProjects\EmpiricalGMMs_Files\output\residuals'
# resultsDirRoot = r'C:\Users\cde84\PycharmProjects\EmpiricalGMMs_Files\output'
# resultsDir = os.path.join(resultsDirRoot, runID, 'residuals')
resultsDir = r'C:\Users\cde84\PycharmProjects\EmpiricalGMMs_Files\output\v1_wCPLB_siteSpecificZVals\residuals'
resultsDir_760 = r'C:\Users\cde84\PycharmProjects\EmpiricalGMMs_Files\output\v1_wCPLB_POTS_760mps_siteSpecificZVals\residuals'


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
df_forHeaders = pd.read_csv(os.path.join(resultsDir, 'CY_14_Crustal_site.csv'), index_col=0)

#SAHeaders = list(df_forHeaders['imName'].values)[1:]
SAHeaders = list(df_forHeaders.columns.values[1:])

#get a list of periods used in the regression from dataframe
periods = [ float(t[4:]) for t in SAHeaders]

print(periods)

#create a list of "full GMM ID" (including GMM ID and tec type) from files in resultsDir
filePaths = glob(os.path.join(resultsDir, '*_site.csv'))
GMM_IDs_List = [path.split('\\')[-1][:-9] for path in filePaths]

station = 'POTS'

#create list to save all events for a station to compute mean and sigma
dS2S_AllGMPEs = []
dS2S_AllGMPEs_760 = []


# get Vs30 and Tsite from site database
Vs30 = df_siteDB['Vs30'][df_siteDB['sta'] == station].values[0]
Tsite = df_siteDB['Tsite'][df_siteDB['sta'] == station].values[0]

# get geomorphic catergory from dataframe
geomCat = df_geom['Updated Geomorphology'][df_geom['stat_id'] == station].values[0]

#change plotting settings
plt.style.use('classic')
plt.rcParams["font.family"] = "Times New Roman"

#create figure
fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True)

#Loop through different GMPE IDs (full ID including Tec Type)
for iGMPE, GMPE in enumerate(GMM_IDs_List):

    # get the GMPE ID without tec type
    gmpe_ID = '_'.join(GMPE.split('_')[:-1])
    gmpe_Tectype = GMPE.split('_')[-1]

    if NGAonly and gmpe_ID not in ['Br_13', 'BSSA_14', 'ASK_14', 'CY_14', 'CB_14']:
        continue

    print('gmpe_ID = %s' %gmpe_ID)
    print('GMPE = %s' % GMPE)

    #load a dataframe with GMPE site terms
    df_dS2S = pd.read_csv(os.path.join(resultsDir, '%s_site.csv' % GMPE), index_col=0)

    df_dS2S_760 = pd.read_csv(os.path.join(resultsDir_760, '%s_site.csv' % GMPE), index_col=0)

    # extract site residual for station from df
    dS2S = df_dS2S.loc[station].values[1:]

    dS2S_760 = df_dS2S_760.loc[station].values[1:]

    #plot dS2S for GMPE
    if iGMPE == 0:
        pLabel = 'Individual GMM'
    else:
        pLabel = None
    axs[0].semilogx(periods, dS2S, color='0.5', label=pLabel)
    axs[1].semilogx(periods, dS2S_760, color='0.5', label=pLabel)

    #append the results to list of all GMPE results for each site
    dS2S_AllGMPEs.append(dS2S)
    dS2S_AllGMPEs_760.append(dS2S_760)

#compute median and sigmas across all events
#print('ratios_AllEvents_SA = %s' %ratios_AllEvents_SA)
dS2S_AllGMPEs = np.vstack(dS2S_AllGMPEs).astype(float)
dS2S_Median_AllGMPEs = np.nanmean(dS2S_AllGMPEs, axis=0)
dS2S_Sigma_AllGMPEs = np.nanstd(dS2S_AllGMPEs, axis=0)

dS2S_AllGMPEs_760 = np.vstack(dS2S_AllGMPEs_760).astype(float)
dS2S_Median_AllGMPEs_760 = np.nanmean(dS2S_AllGMPEs_760, axis=0)
dS2S_Sigma_AllGMPEs_760 = np.nanstd(dS2S_AllGMPEs_760, axis=0)

axs[0].semilogx(periods, dS2S_Median_AllGMPEs, color='k', zorder=50, lw=3, label='Geomean')
axs[1].semilogx(periods, dS2S_Median_AllGMPEs_760, color='k', zorder=50, lw=3, label='Geomean')


axs[1].fill_between(periods, np.array(dS2S_Median_AllGMPEs_760 + dS2S_Sigma_AllGMPEs_760),
                                               np.array(dS2S_Median_AllGMPEs_760 - dS2S_Sigma_AllGMPEs_760), color='0.5', alpha=0.4, zorder=3,
                                               label='$\pm$ Between-model $\sigma _{B-m,s}$')

axs[0].fill_between(periods, np.array(dS2S_Median_AllGMPEs + dS2S_Sigma_AllGMPEs),
                                               np.array(dS2S_Median_AllGMPEs - dS2S_Sigma_AllGMPEs), color='0.5', alpha=0.4, zorder=3,
                                               label='$\pm$ Between-model $\sigma$')



for ax in axs:
    ax.set_ylim(-1, 1)
    ax.set_xlim(0.01, 10)

    ax.yaxis.grid(True)
    #ax.xaxis.grid(True)
    ax.tick_params(which='major', direction='inout', length=10)
    ax.tick_params(which='minor', direction='inout', length=5)
    ax.tick_params(axis='both', labelsize=16)
    ax.axhline(0, color='k')

fig.text(0.013, 0.55, 'Site-to-site residual, $\delta S2S$', ha='left', va='center', rotation='vertical', fontsize=16)
axs[0].text(0.5, 0.96, 'POTS with $V_{S30} = 453$ m/s', transform=axs[0].transAxes, fontsize=14, va='top', ha='center')
axs[1].text(0.5, 0.96, 'POTS with $V_{S30} = 760$ m/s', transform=axs[1].transAxes, fontsize=14, va='top', ha='center')
axs[1].set_xlabel('Vibration Period (s)', fontsize=16)

axs[0].legend(ncol=1, fancybox=True, loc='best', prop={'size':10}, handlelength=2)
axs[1].text(0.97, 0.03, 'Overprediction', transform=axs[1].transAxes, fontsize=14, va='bottom', ha='right')
axs[1].text(0.97, 0.97, 'Underprediction', transform=axs[1].transAxes, fontsize=14, va='top', ha='right')

# ax.set_ylabel('Site-to-site residual', fontsize=20)

# fig.text(0.54, 0.01, 'Shear Wave Velocity, $V_S$ (m/s)', ha='center', va='bottom', fontsize=20)

fig.tight_layout()
fig.subplots_adjust(right=0.95, left=0.107, bottom=0.15, top=0.97, wspace=0.06, hspace=0.1)

fig.savefig(os.path.join(figDir, 'Fig_dS2S_POTS.pdf'))
fig.savefig(os.path.join(figDir, 'Fig_dS2S_POTS.png'), dpi=500)
