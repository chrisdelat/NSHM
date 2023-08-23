import pandas as pd
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

root_path = "C:/Users/cde84/Dropbox (Personal)/PostDocWork/WellingtonSiteChar/LowerHutt_Summer2023/HVSR_dataShared/HVSR_data/All_Cut_Waveforms"

# define directory to save figures
figDir = 'C:/Users/cde84/Dropbox (Personal)/PostDocWork/NSHM_WellingtonBasin/SSRh/figures'
os.makedirs(figDir, exist_ok=True)

metadataPath = os.path.join(root_path, "MetadataHVSR_refStations_AllTests.xlsx")
df = pd.read_excel(metadataPath)

#create a pdf to save all figures
pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join(figDir, 'HVSR_meanCurvesAtSMS.pdf'))

# create a df of mHVSR mean curves for all tests
mHVSR_meanCurvesFilePath = r'C:\Users\cde84\Dropbox (Personal)\PostDocWork\WellingtonSiteChar\LowerHutt_Summer2023\HVSR_dataShared\HVSR_data\All_Cut_Waveforms\HVSRpy\120secWindows_fmin0.05_fmax25.00_peakfmax0.30_peakfmin20.00\HVSR_meanCurve_AllTests_GridAndMAM.txt'
df_mHVSR_meanCurve = pd.read_csv(mHVSR_meanCurvesFilePath, delim_whitespace=True)

# get the headers and frequencies for earthquake SSR from file
mHVSR_meanCurve_headers = list(df_mHVSR_meanCurve.columns.values)[1:]
mHVSR_meanCurve_freqs = np.array([float(freq.split('_')[0]) for freq in mHVSR_meanCurve_headers]).astype(float)

# create a list of all SMS
SMSList = ['NBSS', 'TAIS', 'PVCS', 'SEVS', 'LRSS', 'PGMS', 'SOCS', 'LNBS', 'LHES']

for SMS_ID in SMSList:

    # get all the testIDs corresponding to SMS
    testIDs_SMS = df['Test_ID'][df['Test_ID'].str.contains(SMS_ID)].values

    # create a figure to plot all mean curves
    fig, ax = plt.subplots(nrows=1, ncols=1)

    # print(testIDs_SMS)

    for testID in testIDs_SMS:

        # extract the quality flag from df and skip test if bad quality
        qualityFlag = int(df[df['Test_ID'] == testID]['GoodQuality'])
        if not qualityFlag:
            print('skipping test %s due to bad quality' % testID)
            continue

        print(testID)
        # extract the mHVSR mean curve
        mHVSR_meanCurve_testID = np.array(df_mHVSR_meanCurve[df_mHVSR_meanCurve['testID'] == testID].values[0][1:]).astype(float)

        ax.semilogx(mHVSR_meanCurve_freqs, mHVSR_meanCurve_testID, label=testID)

    ax.set_title('Tests at SMS: %s' %SMS_ID)
    ax.set_ylabel('HVSR Mean Curve')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylim(0,10)
    ax.set_xlim(0.1, 25)
    ax.legend(loc='best', fancybox=True, prop={'size': 10}, framealpha=0.5)

    pdf.savefig(fig)
    plt.close(fig)

pdf.close()