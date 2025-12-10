import numpy as np 
import pandas as pd
import curie as ci
import matplotlib.pyplot as plt 
import os
import sys


sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
from nuclearanalysistools.Tendl import *

pathToSpectra = os.getcwd() + '/' # + '/experimentSpectra/'
pathToCalibrationFiles = os.getcwd() + '/generatedFiles/calibrationFiles/'
print(pathToCalibrationFiles)

class AnalyzeSpectrum:
    def __init__(self, detector, calibrationFile, spectrumFileName, listOfIsotopes):
        self.calibrationFile = calibrationFile
        self.detector = detector
        self.spectrumFileName = spectrumFileName
        if spectrumFileName is not None:
            self.spectrumName = spectrumFileName.replace(".Spe", "")
        self.listOfIsotopes = listOfIsotopes
        
    def analyze(self):
        try:
            cb = ci.Calibration(pathToCalibrationFiles + self.calibrationFile)
            sp = ci.Spectrum(pathToSpectra + self.detector + '/' + self.spectrumFileName)
            sp.cb = cb
            sp.isotopes = self.listOfIsotopes
            sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0}
            # a = sp.fit_peaks(xrays=True)
            # print(a)
            sp.saveas('./test/' + self.spectrumName + '_peak_data.csv', replace=False)
            sp.saveas('./test/' + self.spectrumName + '_peak_data.pdf')
            sp.plot()
        except AttributeError:
            print("Most likely no peaks...")

    def analyze_jobs(self, listOfJobSpectra, peakSummaryFilename):
        cb = ci.Calibration(pathToCalibrationFiles + self.calibrationFile)
        spectra = []
        for job in listOfJobSpectra:
            sp = ci.Spectrum(pathToSpectra + self.detector + '/' + job)
            spectra.append(sp)

        summed_spectrum = spectra[0]
        for i, spec in enumerate(spectra[1:], start=1):
            summed_spectrum += spec
            print(f"Added spectrum {i+1}/{len(spectra)}")

        # print(type(summed_spectrum))
        summed_spectrum.cb = cb
        summed_spectrum.isotopes = self.listOfIsotopes
        print(self.listOfIsotopes)
        summed_spectrum.fit_config = {'SNR_min':3.5}
        # summed_spectrum.plot()
        # print('generatedFiles/peakData/' + peakSummaryFilename + '_peak_data.csv')
        try:
            summed_spectrum.saveas('generatedFiles/peakData/' + peakSummaryFilename + '_peak_data.csv', replace=False)
            summed_spectrum.saveas('gener  atedFiles/peakData/' + peakSummaryFilename + '_peak_data.pdf')
        except AttributeError:
            print("Most likely no peaks...")
# idm_spes = os.getcwd() + '/IDM/idm_spes.txt' 
# det2_spes = os.getcwd() + '/DET2/det2_spes.txt'
# leps_spes = os.getcwd() + '/LEPS/leps_spes.txt'

# idm_filenames = np.genfromtxt(idm_spes, dtype=str, delimiter="\n")
# det2_filenames = np.genfromtxt(det2_spes, dtype=str, delimiter="\n")
# leps_filenames = np.genfromtxt(leps_spes, dtype=str, delimiter="\n")


def analyzeSeveralSpectra(filenames, detector, detectorfile, foil):
    numbers = range(1, 15)  # 1–14 inkludert
    listOfFoils = [f"{foil}{i:02d}" for i in numbers]
    for j in listOfFoils:
        # print(j)
        listOfIsotopes = getListOfIsotopesPerFoil(j)
        # print(listOfIsotopes)
        for file in filenames:
            # listOfIsotopes = getListOfIsotopesPerFoil(foil)
            # print(listOfIsotopes)
            if '_' + j in file:
                filename = file.replace(".Spe", "")
                # print(filename)
                expectedPeakDataFilename = 'generatedFiles/peakData/' + filename + '_peak_data.csv'
                if not os.path.exists(expectedPeakDataFilename):
                    print('trying to analyze file: ' + expectedPeakDataFilename)
                    AnalyzeSpectrum(detector, detectorfile, file, listOfIsotopes).analyze()


# listOfIsotopesNi = ['61CU', '58CO', '56CO', '64CU','62CU','60CU', '63NI', '57NI', '56NI', 'CO-62m','61CO', 'CO-58m', '57CO', '56CO','55CO','59FE', '55FE', '56MN', '54MN', '52MN', '51MN']
listOfIsotopesNi = [
    '61CUg', '62CU', '63CU', '64CU', '65CU',
    '56CO', '57CO', '58CO', '59CO', '60CO', '61CO', '62CO',
    '56NI', '57NI', '58NI', '59NI', '60NI', '61NI', '62NI', '63NI', '64NI',
    '55FE', '56FE', '57FE', '59FE',
    '54MN', '55MN', '56MN', '57MN',
    '52CR', '53CR', '54CR',
    'CO-58m', 'CO-60m'
]


listOfIsotopesCu = [
    '63ZN', '64ZN', '65ZN', '66ZN', '67ZN', '68ZN',
    '62CU', '63CU', '64CU', '65CU', '66CU', '67CU',
    '61NI', '62NI', '63NI', '64NI', '65NI',
    '60CO', '61CO', '62CO', '63CO',
    '58FE', '59FE'
]

listOfIsotopesSn = [
    '112IN', '113IN', '114IN', '115IN', '116IN', '117IN', '118IN', '119IN', '120IN', '121IN',
    '112SN', '113SN', '114SN', '115SN', '116SN', '117SN', '118SN', '119SN', '120SN', '121SN', '122SN', '123SN', '124SN',
    '110CD', '111CD', '112CD', '113CD', '114CD', '115CD', '116CD',
    '111INm', '113INm', '115INm', '117INm', '119INm',
    '115SB', '116SB', '117SB', '118SB', '119SB', '120SB', '121SB', '122SB',
    '120TE', '121TE', '122TE'
]

listOfIsotopesTa = [
    # W fra 181Ta(p,xn)
    '182W', '181W', '180W', '179W', '178W', '177W',

    # Ta fra 181Ta(p,pxn)
    '181TA', '180TA', '179TA', '178TA', '177TA',

    # Hf fra 181Ta(p,αxn)
    '179HF', '178HF', '177HF', '176HF', '175HF',

    # Svakere kanaler videre ned i Z (mulig ved 40–55 MeV, små tverrsnitt)
    '177LU', '176LU', '175LU'
]



def getListOfIsotopesPerFoil(foil, giveActivity=False):
    uniqueIsotopes_55 = pd.read_csv('isotopes_from_simulation/Unique_isotopes_55MeV.csv')
    uniqueIsotopes_30 = pd.read_csv('isotopes_from_simulation/Unique_isotopes_30MeV.csv')
    uniqueIsotopes = pd.concat([uniqueIsotopes_55, uniqueIsotopes_30])
    uniqueIsotopes = uniqueIsotopes[uniqueIsotopes['Activity foil'].astype(str).str.contains(foil, case=False, na=False)]
    if giveActivity:
        return uniqueIsotopes
    return list(uniqueIsotopes['isotope'])



# analyzeSeveralSpectra(det2_filenames, 'DET2', 'eng_cb_det2.json', listOfIsotopesNi,  'Ni')
# analyzeSeveralSpectra(idm_filenames, 'IDM', 'eng_cb_idm.json', listOfIsotopesNi,  'Ni')
# analyzeSeveralSpectra(leps_filenames, 'LEPS', 'eng_cb_leps.json', listOfIsotopesNi,  'Ni')


# analyzeSeveralSpectra(det2_filenames, 'DET2', 'eng_cb_det2.json', listOfIsotopesCu,  'Cu')
# analyzeSeveralSpectra(idm_filenames, 'IDM', 'eng_cb_idm.json', listOfIsotopesCu,  'Cu')
# analyzeSeveralSpectra(leps_filenames, 'LEPS', 'eng_cb_leps.json', listOfIsotopesCu,  'Cu')

# analyzeSeveralSpectra(det2_filenames, 'DET2', 'eng_cb_det2.json', listOfIsotopesTa,  'Ta')
# analyzeSeveralSpectra(idm_filenames, 'IDM', 'eng_cb_idm.json', listOfIsotopesTa,  'Ta')
# analyzeSeveralSpectra(leps_filenames, 'LEPS', 'eng_cb_leps.json', listOfIsotopesTa,  'Ta')

# analyzeSeveralSpectra(det2_filenames, 'DET2', 'eng_cb_det2.json', listOfIsotopesSn,  'Sn')
# analyzeSeveralSpectra(idm_filenames, 'IDM', 'eng_cb_idm.json', listOfIsotopesSn,  'Sn')
# analyzeSeveralSpectra(leps_filenames, 'LEPS', 'eng_cb_leps.json', listOfIsotopesSn,  'Sn')





def overviewPeaks(peakDataFiles, foil):
    uniqueIsotopes_55 = pd.read_csv('isotopes_from_simulation/Unique_isotopes_55MeV.csv')
    uniqueIsotopes_30 = pd.read_csv('isotopes_from_simulation/Unique_isotopes_30MeV.csv')
    uniqueIsotopes = pd.concat([uniqueIsotopes_55, uniqueIsotopes_30])
    uniqueIsotopes = uniqueIsotopes[uniqueIsotopes['Activity foil'].astype(str).str.contains(foil, case=False, na=False)]
    df_uniqueIsotopes = uniqueIsotopes[['isotope', 'A0 bequerel']]
    print(df_uniqueIsotopes)

    # summary = ", ".join(
    # f"{row['isotope']} ({row['A0 bequerel']:.1f} Bq)"
    # for _, row in uniqueIsotopes.iterrows())
    # # print(summary)

    
    filteredPeakDataFiles = []
    for file in peakDataFiles:
        if '_' + foil in file:
            print("_" + foil)
            filteredPeakDataFiles.append(file)
    print(filteredPeakDataFiles)

    
    newdata = []
    
    listOfIsotopes = getListOfIsotopesPerFoil(foil)
    for isotope in listOfIsotopes:
        for f in filteredPeakDataFiles:
            f = 'generatedFiles/peakData/' + f
            df = pd.read_csv(f)
            isotopes = df['isotope']; counts = df['counts']; filename = df['filename'];energy = df['energy'];intensity = df['intensity']
            live_time = df['live_time']
            start_time = df['start_time']
            timeSinceEob = getDelayTime(foil, start_time)[0]
            for idx in range(len(isotopes)):
                if isotope == isotopes.iloc[idx]:
                    match = uniqueIsotopes.loc[uniqueIsotopes['isotope'].eq(isotope), 'A0 bequerel']
                    activity_bq = int(round(match.iloc[0])) if not match.empty else None
                    formatted_livetime = format_time(live_time.iloc[idx])
                    half_life = ci.Isotope(isotope).half_life()
                    formatted_halflife = format_time(half_life)
                    newdata.append([foil, isotopes.iloc[idx], formatted_halflife, timeSinceEob, activity_bq, energy.iloc[idx], counts.iloc[idx], intensity.iloc[idx], formatted_livetime,  filename.iloc[idx]])
    df_isotope = pd.DataFrame(newdata, columns=['foil', 'isotope', 'half life', 'time since eob', 'expected A0 (Bq)', 'energy', 'counts', 'intensity', 'live_time', 'filename'])
    df_isotope.to_csv('peakSummaries/' + foil + '_summary.csv')
        

def format_time(t_seconds):
    if t_seconds < 60:
        return f"{t_seconds:.2f} s"
    elif t_seconds < 3600:
        return f"{t_seconds/60:.2f} min"
    elif t_seconds < 86400:
        return f"{t_seconds/3600:.2f} h"
    elif t_seconds < 31557600:  # ~1 year (365.25 days)
        return f"{t_seconds/86400:.2f} d"
    elif str(t_seconds) == 'nan':
        return "-"
    else:
        return f"{t_seconds/31557600:.2f} y"   
    
def getDelayTime(foil, start_time):
    # Velg riktig EOB-tid
    try:
        if int(foil[-2:]) <= 7:
            eob_time = pd.Timestamp("2025-09-24 14:32:54")
        else:
            eob_time = pd.Timestamp("2025-09-23 14:13:24")

        start_time = pd.to_datetime(start_time)
        delay = start_time - eob_time
        return delay
    except:
        return '-'


# ni_count_files = os.getcwd() + '/generatedFiles/peakData/ni_files.txt'
# ni_count_filenames = np.genfromtxt(ni_count_files, dtype=str, delimiter="\n")

# cu_count_files = os.getcwd() + '/generatedFiles/peakData/cu_files.txt'
# cu_count_filenames = np.genfromtxt(cu_count_files, dtype=str, delimiter="\n")

# ta_count_files = os.getcwd() + '/generatedFiles/peakData/ta_files.txt'
# ta_count_filenames = np.genfromtxt(ta_count_files, dtype=str, delimiter="\n")

# sn_count_files = os.getcwd() + '/generatedFiles/peakData/sn_files.txt'
# sn_count_filenames = np.genfromtxt(sn_count_files, dtype=str, delimiter="\n")


# overviewPeaks(ni_count_filenames, listOfIsotopesNi, 'Ni05')

def generatePeakOverviews(element, count_filenames):
    numbers = range(1, 15)  # 1–14 inkludert
    listOfFoils = [f"{element}{i:02d}" for i in numbers]
    
    for i in listOfFoils:
    # for i in ['Ni05', 'Ni06', 'Ni07', 'Ni08']:
        overviewPeaks(count_filenames, i)


# generatePeakOverviews('Ni', ni_count_filenames)
# generatePeakOverviews('Cu', cu_count_filenames)
# generatePeakOverviews('Ta', ta_count_filenames)
# generatePeakOverviews('Sn', sn_count_filenames)

# print(getListOfIsotopesPerFoil('Cu01', giveActivity=True))

# natTa = Tendl({'Ta181': 1.0}, 'proton')
# natNi = Tendl({'Ni64': 0.009256, 'Ni62': 0.036345, 'Ni61': 0.011399, 'Ni60': 0.262231, 'Ni58': 0.680769}, 'proton')
# natCu = Tendl({'Cu63': 0.6915, 'Cu65': 0.3085}, 'proton')
# W = '74'; Ta='73'; Hf='72'; Lu='71'; Yb='70'
# natTa.plotTendl23Unique(Hf, '173') #, isomerLevel='00')


# natNi.plotTendl23Unique('27', '55', color='red', label='55Co')
# natNi.plotTendl23Unique('28', '56', color='blue', label='56Ni')
# natNi.plotTendl23Unique('28', '57', color='green', label='57Ni')
# natNi.plotTendl23Unique('27', '58', color='yellow', label='58Co')
# natNi.plotTendl23Unique('27', '56', color='orange', label='56Co')
# natNi.plotTendl23Unique('27', '57', color='cyan', label='57Co')



# natCu.plotTendl23Unique('30', '65', color='red', label='Zn65')
# natCu.plotTendl23Unique('27', '56', color='blue', label='56Co')
# natCu.plotTendl23Unique('27', '58', color='green', label='58Co')
# natCu.plotTendl23Unique('30', '63', color='yellow', label='63Zn')
# natCu.plotTendl23Unique('30', '62', color='orange', label='62Zn')
# natCu.plotTendl23Unique('29', '61', color='cyan', label='61Cu')
# plt.title('Cu monitors')
# plt.legend()
# plt.show()