import numpy as np
from glob import glob
import os.path
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as colors

def get_depth_versus_Vs_for_step_plot(layerThicks, layerVs, startDepth=0):

    #get the depth to bottom of layers
    depth_BoL = np.cumsum(layerThicks)

    depth4plot = np.concatenate([[x, x] for x in depth_BoL])
    depth4plot = np.insert(depth4plot, 0, 0)
    depth4plot = np.delete(depth4plot, -1)

    depth4plot += startDepth

    # make array of Vs for plotting stepwise Vs

    Vs4plot = np.concatenate([[x, x] for x in layerVs])

    return Vs4plot, depth4plot


# specify path to Vs profiles
VsDir = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\WellingtonSiteResponse\VsProfiles_forPaperPlots'

#specify directory for saving figs
figDir = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\NSHM_WellingtonBasin\JournalPaper_WelliSiteResp\nonlinearBasinEffects\R1_Submission2\figs'

#create a dataframe of the site database based on Jesse's DB v3.2
welliDB_Path = r'C:\Users\cde84\PycharmProjects\EmpiricalGMMs_Files\gm_DB_input\site_table_v3.2.csv'
df_siteDB = pd.read_csv(welliDB_Path)

#change plotting settings
plt.style.use('classic')
plt.rcParams["font.family"] = "Times New Roman"

#create plotting colors for basins
#pColorsList = ['blue', 'maroon', 'yellow', 'Orange']
pColorsList_Hex = ['#4363d8', '#800000', '#ffe119', '#f58231']
#pColorsList = [pygame.Color(pColor) for pColor in pColorsList_Hex]
pColorsList = [colors.hex2color(pColor) for pColor in pColorsList_Hex]

# create a figure
fig, axs = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=True, figsize=(9.5,11))
axs = np.concatenate(axs)

# create a list of siteIDs to plot
siteIDs = ['WNKS', 'FKPS', 'VUWS', 'WEMS', 'TEPS', 'PIPS', 'WNAS', 'MISS', 'CPLB']

# loop through sites
for iSite, siteID in enumerate(siteIDs):
    # Glob paths for Vs files for each sites
    VsPaths = glob(os.path.join(VsDir, '%s*' %siteID))

    print(VsPaths)

    # get Vs30 and site period
    Vs30 = df_siteDB['Vs30'][df_siteDB['sta'] == siteID].values[0]
    Tsite = df_siteDB['Tsite'][df_siteDB['sta'] == siteID].values[0]

    # loop through Vs profile paths
    for iVs, VsPath in enumerate(VsPaths):
        VsArray = np.loadtxt(VsPath, skiprows=1).T
        layerThick = VsArray[0]
        Vs = VsArray[1]

        plotVs, plotDepths = get_depth_versus_Vs_for_step_plot(layerThick, Vs)
        print(plotVs)
        print(plotDepths)
        axs[iSite].plot(plotVs, plotDepths, color=pColorsList[iVs], lw=2)

    axs[iSite].invert_yaxis()
    #axs[iSite].text(0.03, 0.03, '%s' % (siteID),
    #                    transform=axs[iSite].transAxes, fontsize=16, va='bottom', ha='left')
    axs[iSite].text(0.03, 0.03, '%s\n$T_{site}=$ %ss\n$V_{S30}=$ %s m/s' % (siteID, Tsite, Vs30),
                        transform=axs[iSite].transAxes, fontsize=18, va='bottom', ha='left')
    axs[iSite].set_xlim(0, 550)

for ax in  axs:
    ax.yaxis.grid(True)
    ax.xaxis.grid(True)
    ax.tick_params(which='major', direction='inout', length=10)
    ax.tick_params(which='minor', direction='inout', length=5)
    ax.tick_params(axis='both', labelsize=20)

axs[3].set_ylabel('Depth (m)', fontsize=20)
fig.text(0.54, 0.01, 'Shear Wave Velocity, $V_S$ (m/s)', ha='center', va='bottom', fontsize=20)

fig.tight_layout()
fig.subplots_adjust(right=0.98, left=0.1, bottom=0.08, top=0.99, wspace=0.06, hspace=0.07)

fig.savefig(os.path.join(figDir, 'Fig_VsProfiles.pdf'))
fig.savefig(os.path.join(figDir, 'Fig_VsProfiles.png'), dpi=500)

plt.close(fig)

#create a figure for POTS
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5.5,8))

siteID = 'POTS'
VsPath = glob(os.path.join(VsDir, '%s*' %siteID))[0]

# get Vs30 and site period
Vs30 = df_siteDB['Vs30'][df_siteDB['sta'] == siteID].values[0]
Tsite = df_siteDB['Tsite'][df_siteDB['sta'] == siteID].values[0]

# load Vs profile path
VsArray = np.loadtxt(VsPath, skiprows=1).T
layerThick = VsArray[0]
Vs = VsArray[1]

plotVs, plotDepths = get_depth_versus_Vs_for_step_plot(layerThick, Vs, 1)
print(plotVs)
print(plotDepths)
ax.plot(plotVs, plotDepths, color='k', lw=3)

ax.invert_yaxis()
#axs[iSite].text(0.03, 0.03, '%s' % (siteID),
#                    transform=axs[iSite].transAxes, fontsize=16, va='bottom', ha='left')
ax.text(0.03, 0.03, '%s\n$T_{site}=$ %ss\n$V_{S30}=$ %s m/s' % (siteID, Tsite, Vs30),
                    transform=ax.transAxes, fontsize=18, va='bottom', ha='left')
ax.set_xlim(0, 1200)

ax.yaxis.grid(True)
ax.xaxis.grid(True)
ax.tick_params(which='major', direction='inout', length=10)
ax.tick_params(which='minor', direction='inout', length=5)
ax.tick_params(axis='both', labelsize=16)

ax.set_ylabel('Depth (m)', fontsize=16)
ax.set_xlabel('Shear Wave Velocity, $V_S$ (m/s)', fontsize=16)
# fig.text(0.54, 0.01, 'Shear Wave Velocity, $V_S$ (m/s)', ha='center', va='bottom', fontsize=20)

fig.tight_layout()
fig.subplots_adjust(right=0.95, left=0.18, bottom=0.12, top=0.97, wspace=0.06, hspace=0.07)

fig.savefig(os.path.join(figDir, 'Fig_VsProfile_POTS.pdf'))
fig.savefig(os.path.join(figDir, 'Fig_VsProfiles_POTS.png'), dpi=500)