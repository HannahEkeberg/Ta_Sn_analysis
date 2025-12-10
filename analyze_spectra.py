import numpy as np 
import pandas as pd
import curie as ci
import matplotlib.pyplot as plt 
import os
import sys


sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
from nuclearanalysistools.Tendl import *

pathToSpectra = os.getcwd() + '/spectra/'
pathToCalibrationFiles = os.getcwd() + '/generatedfiles/calibrationfiles/'
pathToPeakData = os.getcwd() + '/generatedfiles/peakdata/'

class AnalyzeSpectrum:
    def __init__(self, detector, calibrationFile):
        self.calibrationFile = calibrationFile
        self.detector = detector
        
    def analyze(self, spectrumFileName, listOfIsotopes):
        spectrumFileName = spectrumFileName
        if spectrumFileName is not None:
            spectrumName = spectrumFileName.replace(".Spe", "")
        self.listOfIsotopes = listOfIsotopes
        try:
            cb = ci.Calibration(pathToCalibrationFiles + self.calibrationFile)
            sp = ci.Spectrum(pathToSpectra + self.detector + '/' + spectrumFileName)
            sp.cb = cb
            sp.isotopes = listOfIsotopes
            sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0}
            sp.saveas(pathToPeakData + '/data/' + spectrumName + '_peak_data.csv', replace=False)
            sp.saveas(pathToPeakData + '/figures/' + spectrumName + '_peak_data.pdf')
            sp.plot()
        except AttributeError:
            print("Most likely no peaks...")

    def analyze_jobs(self, listOfJobSpectra, peakSummaryFilename, listOfIsotopes):
        if peakSummaryFilename == None:
            peakSummaryFilename = listOfJobSpectra[0].replace("000.Spe", "job")
        cb = ci.Calibration(pathToCalibrationFiles + self.calibrationFile)
        spectra = []
        for job in listOfJobSpectra:
            sp = ci.Spectrum(pathToSpectra + self.detector + '/' + job)
            spectra.append(sp)

        summed_spectrum = spectra[0]
        for i, spec in enumerate(spectra[1:], start=1):
            summed_spectrum += spec
            print(f"Added spectrum {i+1}/{len(spectra)}")

        summed_spectrum.cb = cb
        summed_spectrum.isotopes = listOfIsotopes
        summed_spectrum.fit_config = {'SNR_min':3.5}
        summed_spectrum.plot()
        try:
            summed_spectrum.saveas('generatedFiles/peakData/data/' + peakSummaryFilename + '_peak_data.csv', replace=False)
            summed_spectrum.saveas('generatedFiles/peakData/figures/' + peakSummaryFilename + '_peak_data.pdf')
        except AttributeError:
            print("Most likely no peaks...")

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

def getListOfIsotopesPerFoil(foil=None, giveActivity=False):
    uniqueIsotopes_55 = pd.read_csv('isotopes_from_simulation/Unique_isotopes_55MeV.csv')
    uniqueIsotopes_30 = pd.read_csv('isotopes_from_simulation/Unique_isotopes_30MeV.csv')
    uniqueIsotopes = pd.concat([uniqueIsotopes_55, uniqueIsotopes_30])
    uniqueIsotopes = uniqueIsotopes[uniqueIsotopes['Activity foil'].astype(str).str.contains(foil, case=False, na=False)]
    if giveActivity:
        return uniqueIsotopes
    return list(uniqueIsotopes['isotope'])

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

def generatePeakOverviews(element, count_filenames):
    numbers = range(1, 15)  # 1–14 inkludert
    listOfFoils = [f"{element}{i:02d}" for i in numbers]
    
    for i in listOfFoils:
    # for i in ['Ni05', 'Ni06', 'Ni07', 'Ni08']:
        overviewPeaks(count_filenames, i)

def get_spectra(detector_folder, distance, foil, file_ending='.Spe'):
    root = os.getcwd() + '/spectra/' + detector_folder
    filenames = []
    for filename in os.listdir(root):
        if distance in filename and foil in filename and file_ending in filename:
            filenames.append(filename)

    d = {}
    for f in sorted(filenames):
        d.setdefault(f[:2], []).append(f)
    return [v[0] if len(v) == 1 else v for _, v in sorted(d.items())]