import sys
import numpy as np 
import curie as ci 
import pandas as pd
import matplotlib.pyplot as plt
import os

sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
pathToCalibrationSpectra = os.getcwd() + '/Calibration'
print(pathToCalibrationSpectra)

class Calibration:

    def __init__(self, pathToCalibrationSpectrum, detector, spec_Ba133, spec_Cs137, spec_Eu152, detectorCalibrationName):
        self.pathToCalibrationSpectrum = (pathToCalibrationSpectrum + '/' + detector + '/')
        self.spec_Ba133 = spec_Ba133; self.spec_Cs137 = spec_Cs137; self.spec_Eu152 = spec_Eu152#; self.spec_Am141 = spec_Am141
        self.detectorCalibrationName = detectorCalibrationName

    def sources(self):
        sources = [
            {'isotope':'133BA', 'A0':3.989E4, 'ref_date':'01/01/2009 12:00:00'},
            {'isotope':'137CS', 'A0':3.855E4, 'ref_date':'01/01/2009 12:00:00'},
            {'isotope':'152EU', 'A0':3.3929E4, 'ref_date':'01/01/2009 12:00:00'}
            # {'isotope':'152EU', 'A0':3.3929E4, 'ref_date':'01/01/2009 12:00:00'}
        ]
        return pd.DataFrame(sources)
    
    def sp_Ba133(self):
        sp_Ba133 = ci.Spectrum(self.pathToCalibrationSpectrum + self.spec_Ba133)
        sp_Ba133.isotopes = ['133BA']
        return sp_Ba133
    
    def sp_Cs137(self):
        sp_Cs137 = ci.Spectrum(self.pathToCalibrationSpectrum + self.spec_Cs137)
        sp_Cs137.isotopes = ['137CS']
        return sp_Cs137
    
    def sp_Eu152(self):
        sp_Eu152 = ci.Spectrum(self.pathToCalibrationSpectrum + self.spec_Eu152)
        sp_Eu152.isotopes = ['152EU']
        return sp_Eu152
    
    def calibrate(self, showPlot=False, detector=None):
        cb = ci.Calibration()
        if detector == 'IDM':
            cb.calibrate([self.sp_Ba133(), self.sp_Cs137(), self.sp_Eu152()], sources=self.sources())
            # cb.calibrate([self.sp_Ba133(), self.sp_Eu152()], sources=self.sources()) # only for energy cal det2
            # 
        if detector == 'DET2':
            # cb.calibrate([self.sp_Ba133(), self.sp_Cs137(), self.sp_Eu152()], sources=self.sources())
            cb.calibrate([self.sp_Ba133(), self.sp_Eu152()], sources=self.sources()) # only for energy cal det2 
        if detector=='LEPS':
            # cb.calibrate([self.sp_Ba133(), self.sp_Cs137(), self.sp_Eu152()], sources=self.sources())
            cb.calibrate([self.sp_Ba133()], sources=self.sources()) # only for energy cal det2 
        if showPlot:
            cb.plot(show=True, saveas = os.getcwd() + '/generatedFiles/calibrationFiles/' + self.detectorCalibrationName + 'pdf')
            
        # sp = ci.Spectrum(self.sp_Ba133)
        # sp.cb = cb
        # # sp.isotopes = self.listOfIsotopes
        # sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0}
        # # a = sp.fit_peaks(xrays=True)
        # # print(a)
        # # sp.saveas('generatedFiles/peakData/' + self.spectrumName + '_peak_data.csv', replace=False)
        # # sp.saveas('generatedFiles/peakData/' + self.spectrumName + '_peak_data.png')
        # sp.plot()
        cb.saveas(os.getcwd() + '/generatedFiles/calibrationFiles/' + self.detectorCalibrationName + '.json')


# cb_energy_idm = Calibration(pathToCalibrationSpectra, 'IDM', 'AB20250922_IDM_Ba133_30cm.Spe', 'AJ20250923_IDM_Cs137_35cm.Spe', 'AE20250922_IDM_Eu152_36cm.Spe', 'eng_cb_idm')
# cb_energy_det2 = Calibration(pathToCalibrationSpectra, 'DET2', 'AN20250923_Det2_Eu152_60cm.Spe', 'AF20250922_Det2_Ba133_60cm.Spe', 'AA20250922_Det2_Eu152_18cm.Spe', 'eng_cb_det2')
# cb_energy_leps = Calibration(pathToCalibrationSpectra, 'LEPS', 'AM20250923_LEPS_Ba133_30cm.Spe', 'AJ20250923_IDM_Cs137_35cm.Spe', 'AE20250922_IDM_Eu152_36cm.Spe', 'eng_cb_leps')

# cb_energy_idm.calibrate()
# cb_energy_idm.calibrate('IDM')
# cb_energy_det2.calibrate('DET2')
# cb_energy_leps.calibrate('LEPS')




def calibration_idm_30cm():
    cb = ci.Calibration()

    sp_Ba133 = ci.Spectrum('./Calibration/AB20250922_IDM_Ba133_30cm.Spe')
    sp_Ba133.isotopes = ['133BA'] 
    # sp_Ba133.plot()

    sp_Cs137 = ci.Spectrum('./Calibration/AK20250923_IDM_Cs137_30cm.Spe')
    sp_Cs137.isotopes = ['137CS'] 
    # sp_Cs137_2.plot()

    sp_Eu152 = ci.Spectrum('./Calibration/CD20251018_IDM_Eu152_30cm.Spe')
    sp_Eu152.isotopes = ['152EU']
    # sp_Eu152.plot()

    sources = [{'isotope':'133BA', 'A0':3.989E4, 'ref_date':'01/01/2009 12:00:00'},
               {'isotope':'137CS', 'A0':3.855E4, 'ref_date':'01/01/2009 12:00:00'},
               {'isotope':'152EU', 'A0':370000, 'ref_date':'11/01/1984 12:00:00'}]
    sources = pd.DataFrame(sources)

    cb.calibrate([sp_Ba133, sp_Cs137, sp_Eu152], sources=sources)
    cb.plot(show=True)
    # cb.saveas(os.getcwd() + '/generatedFiles/calibrationFiles/' + 'eng_cb_idm' + '.json')
    # cb.plot(show=True, saveas = './MyGeneratedFiles/Calibration/Figures/calibration_plots_10cm.pdf')

def calibration_det2_10cm():
    cb = ci.Calibration()

    sp_Ba133 = ci.Spectrum('./Calibration/BV20251017_Det2_Ba133_10cm.Spe')
    sp_Ba133.isotopes = ['133BA'] 
    # sp_Ba133.plot()

    sp_Cs137 = ci.Spectrum('./Calibration/BA20251017_Det2_Cs137_10cm.Spe')
    sp_Cs137.isotopes = ['137CS'] 
    # sp_Cs137_2.plot()

    sp_Eu152 = ci.Spectrum('./Calibration/CL20251020_Det2_Eu152_10cm.Spe')
    sp_Eu152.isotopes = ['152EU']
    # sp_Eu152.plot()

    sources = [{'isotope':'133BA', 'A0':3.989E4, 'ref_date':'01/01/2009 12:00:00'},
               {'isotope':'137CS', 'A0':3.855E4, 'ref_date':'01/01/2009 12:00:00'},
               {'isotope':'152EU', 'A0':370000, 'ref_date':'11/01/1984 12:00:00'}]
    sources = pd.DataFrame(sources)

    cb.calibrate([sp_Ba133, sp_Cs137, sp_Eu152], sources=sources)
    cb.plot(show=True)
    cb.saveas(os.getcwd() + '/generatedFiles/calibrationFiles/' + 'eng_cb_det2' + '.json')
    # cb.plot(show=True, saveas = './MyGeneratedFiles/Calibration/Figures/calibration_plots_10cm.pdf')

def calibration_leps_30cm():
    cb = ci.Calibration()

    sp_Ba133 = ci.Spectrum('./Calibration/AM20250923_LEPS_Ba133_30cm.Spe')
    sp_Ba133.isotopes = ['133BA'] 
    # sp_Ba133.plot()

    # sp_Cs137 = ci.Spectrum('./Calibration/BA20251017_Det2_Cs137_10cm.Spe')
    # sp_Cs137.isotopes = ['137CS'] 
    # # sp_Cs137_2.plot()

    # sp_Eu152 = ci.Spectrum('./Calibration/CL20251020_Det2_Eu152_10cm.Spe')
    # sp_Eu152.isotopes = ['152EU']
    # sp_Eu152.plot()

    sources = [{'isotope':'133BA', 'A0':3.989E4, 'ref_date':'01/01/2009 12:00:00'}]
            #    {'isotope':'137CS', 'A0':3.855E4, 'ref_date':'01/01/2009 12:00:00'},
            #    {'isotope':'152EU', 'A0':370000, 'ref_date':'11/01/1984 12:00:00'}]
    sources = pd.DataFrame(sources)

    cb.calibrate([sp_Ba133], sources=sources)
    cb.plot(show=True)
    cb.saveas(os.getcwd() + '/generatedFiles/calibrationFiles/' + 'eng_cb_leps' + '.json')
    # cb.plot(show=True, saveas = './MyGeneratedFiles/Calibration/Figures/calibration_plots_10cm.pdf')

def calibration_leps_30cm():
    cb = ci.Calibration()

    sp_Ba133 = ci.Spectrum('./Calibration/AM20250923_LEPS_Ba133_30cm.Spe')
    sp_Ba133.isotopes = ['133BA'] 
    # sp_Ba133.plot()

    # sp_Cs137 = ci.Spectrum('./Calibration/BA20251017_Det2_Cs137_10cm.Spe')
    # sp_Cs137.isotopes = ['137CS'] 
    # # sp_Cs137_2.plot()

    # sp_Eu152 = ci.Spectrum('./Calibration/CL20251020_Det2_Eu152_10cm.Spe')
    # sp_Eu152.isotopes = ['152EU']
    # sp_Eu152.plot()

    sources = [{'isotope':'133BA', 'A0':3.989E4, 'ref_date':'01/01/2009 12:00:00'}]
            #    {'isotope':'137CS', 'A0':3.855E4, 'ref_date':'01/01/2009 12:00:00'},
            #    {'isotope':'152EU', 'A0':370000, 'ref_date':'11/01/1984 12:00:00'}]
    sources = pd.DataFrame(sources)

    cb.calibrate([sp_Ba133], sources=sources)
    cb.plot(show=True)
    cb.saveas(os.getcwd() + '/generatedFiles/calibrationFiles/' + 'eng_cb_leps' + '.json')
    # cb.plot(show=True, saveas = './MyGeneratedFiles/Calibration/Figures/calibration_plots_10cm.pdf')


def calibration_xray_position4():
    cb = ci.Calibration()

    sp_Co57 = ci.Spectrum('./Calibration/CE10182025_Xray_Co57_HC6928_p4_converted.Spe')
    sp_Co57.isotopes = ['57CO'] 
    sp_Co57.plot()

    sp_Am241 = ci.Spectrum('./Calibration/CH10202025_Xray_Am241_p4_converted.Spe') 
    sp_Am241.isotopes = ['241AM'] 
    sp_Am241.plot()

    sources = [{'isotope':'57CO', 'A0':4.30E4, 'ref_date':'11/01/2023 13:27:02'},
               {'isotope':'241AM', 'A0':3.752E4, 'ref_date':'03/01/2019 12:00:00'}
               ]
    sources = pd.DataFrame(sources)

    cb.calibrate([sp_Co57, sp_Am241], sources=sources)
    cb.plot(show=True)
    cb.saveas(os.getcwd() + '/generatedFiles/calibrationFiles/' + 'eng_cb_leps' + '.json')
    
    
    # cb.plot(show=True, saveas = './MyGeneratedFiles/Calibration/Figures/calibration_plots_10cm.pdf')

# calibration_leps_30cm()
# calibration_det2_10cm()
# calibration_idm_30cm()

calibration_xray_position4()

