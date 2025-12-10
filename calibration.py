import sys
import numpy as np 
import curie as ci 
import pandas as pd
import matplotlib.pyplot as plt
import os

sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
pathToCalibrationSpectra = os.getcwd()
print(pathToCalibrationSpectra)

class Calibration:

    def __init__(self, detector_posititon, spectra, detectorCalibrationName):
        self.pathToCalibrationSpectrum = ('Spectra/Calibration' + '/' + detector_posititon + '/')
        # self.spec_Ba133 = spec_Ba133; self.spec_Cs137 = spec_Cs137; self.spec_Eu152 = spec_Eu152#; self.spec_Am141 = spec_Am141
        self.spectra = spectra
        self.detectorCalibrationName = detectorCalibrationName
        calibrationSpectra = []
        sources = []
        for isotope, spec in self.spectra.items():
            calibrationSpectra.append(self.spectrum(isotope, spec))
            isotope_callable = getattr(self, self.iso_to_method(isotope))
            sources.append(isotope_callable())

        cb = ci.Calibration()
        cb.calibrate(calibrationSpectra, sources=sources)
        # cb.plot(show=True)

        cb.saveas(os.getcwd() + '/generatedfiles/calibrationfiles/' + self.detectorCalibrationName + '.json')
        cb.plot(show=True, saveas = './generatedfiles/calibrationfiles/figures/' + self.detectorCalibrationName +'.pdf')

    def spectrum(self, isotope, spec): # isotope = '133BA'
        if type(spec) is list:
            sp = ci.Spectrum(self.pathToCalibrationSpectrum + spec[0]) + ci.Spectrum(self.pathToCalibrationSpectrum + spec[1])
        else:
            sp = ci.Spectrum(self.pathToCalibrationSpectrum + spec)
        sp.isotopes = [isotope]
        return sp
    
    def iso_to_method(self,name: str) -> str:
        digits = ''.join(ch for ch in name if ch.isdigit())
        letters = ''.join(ch for ch in name if ch.isalpha())
        return letters.capitalize() + digits

    def Ba133(self):
        return {'isotope':'133BA', 'A0':3.3859E4, 'ref_date':'03/01/2019 12:00:00'}
    
    def Cs137(self):
        return {'isotope':'137CS', 'A0':3.67055E4, 'ref_date':'03/01/2019 12:00:00'}
    
    def Eu152(self):
        return {'isotope':'152EU', 'A0':3.822E4, 'ref_date':'03/01/2019 12:00:00'}
    
    def Am241(self):
        return {'isotope':'241AM', 'A0':3.752E4, 'ref_date':'03/01/2019 12:00:00'}
    
    def Cd109(self):
        return {'isotope':'109CD', 'A0':3.774E4, 'ref_date':'03/01/2019 12:00:00'}
    
    def Co57(self):
        return {'isotope':'57CO', 'A0': 5.54E5, 'ref_date': '11/01/2023 13:27:02'}


#det 2:
# spectra_det2_10 = {'137CS': 'BA20251017_Det2_Cs137_10cm.Spe', '133BA': 'BV20251017_Det2_Ba133_10cm.Spe', '152EU': 'CL20251020_Det2_Eu152_10cm.Spe'}
# Calibration('Det2_10', spectra_det2_10, 'calibration_det2_10')
# spectra_det2_18 = {'137CS': 'BC20251017_Det2_Cs137_18cm.Spe', '133BA': 'BX20251017_Det2_Ba133_18cm.Spe', '152EU': 'CK20251020_Det2_Eu152_18cm.Spe'}
# Calibration('Det2_18', spectra_det2_18, 'calibration_det2_18')
# spectra_det2_24 = {'137CS': 'BD20251017_Det2_Cs137_24cm.Spe', '133BA': 'BZ20251017_Det2_Ba133_24cm.Spe', '152EU': 'CM20251020_Det2_Eu152_24cm.Spe'}
# Calibration('Det2_24', spectra_det2_24, 'calibration_det2_24')
# spectra_det2_30 = {'137CS': 'BF20251017_Det2_Cs137_30cm.Spe', '133BA': 'CA20251018_Det2_Ba133_30cm.Spe', '152EU': 'CR20251028_Det2_Eu152_30cm.Spe'}
# Calibration('Det2_30', spectra_det2_30, 'calibration_det2_30')
# spectra_det2_40 = {'137CS': 'BG20251017_Det2_Cs137_40cm.Spe', '133BA': 'CC20251018_Det2_Ba133_40cm.Spe', '152EU': 'CQ20251028_Det2_Eu152_40cm.Spe'}
# Calibration('Det2_40', spectra_det2_40, 'calibration_det2_40')
# spectra_det2_50 = {'137CS': 'BJ20251017_Det2_Cs137_50cm.Spe', '133BA': 'AD20250922_Det2_Ba133_50cm.Spe', '152EU': 'AL20250923_Det2_Eu152_50cm.Spe'}
# Calibration('Det2_50', spectra_det2_50, 'calibration_det2_50')
# spectra_det2_60 = {'137CS': 'BL20251017_Det2_Cs137_60cm.Spe', '133BA': 'AF20250922_Det2_Ba133_60cm.Spe', '152EU': 'CU20251029_Det2_Eu152_60cm.Spe'}
# Calibration('Det2_60', spectra_det2_60, 'calibration_det2_60')
# spectra_det2_70 = {'137CS': 'BM20251017_Det2_Cs137_70cm.Spe', '133BA': 'CG20251018_Det2_Ba133_70cm.Spe', '152EU': 'CY20251030_Det2_Eu152_70cm.Spe'}
# Calibration('Det2_70', spectra_det2_70, 'calibration_det2_70')
# spectra_det2_80 = {'137CS': 'CW20251029_Det2_Cs137_80cm.Spe', '133BA': ['CI20251020_Det2_Ba133_80cm_000.Spe', 'CJ20251020_Det2_Ba133_80cm_001.Spe'], '152EU': 'BU20251017_Det2_Eu152_80cm.Spe'}
# Calibration('Det2_80', spectra_det2_80, 'calibration_det2_80')

# IDM:
# spectra_idm_10 = {'137CS': 'AD20251016_IDM_Cs137_10cm.Spe', '133BA': 'AI20251016_IDM_Ba133_10cm.Spe', '152EU': 'AE20251016_IDM_Eu152_10cm.Spe'}
# Calibration('IDM_10', spectra_idm_10, 'calibration_idm_10')
# spectra_idm_15 = {'137CS': 'BP20251017_IDM_Cs137_15cm.Spe', '133BA': 'BB20251017_IDM_Ba133_15cm.Spe', '152EU': 'BW20251017_IDM_Eu152_15cm.Spe'}
# Calibration('IDM_15', spectra_idm_15, 'calibration_idm_15')
# spectra_idm_20 = {'137CS': 'BQ20251017_IDM_Cs137_20cm.Spe', '133BA': 'BE20251017_IDM_Ba133_20cm.Spe', '152EU': 'BY20251017_IDM_Eu152_20cm.Spe'}
# Calibration('IDM_20', spectra_idm_20, 'calibration_idm_20')
# spectra_idm_25 = {'137CS': 'BR20251017_IDM_Cs137_25cm.Spe', '133BA': 'BH20251017_IDM_Ba133_25cm.Spe', '152EU': 'CB20251018_IDM_Eu152_25cm.Spe'}
# Calibration('IDM_25', spectra_idm_25, 'calibration_idm_25')
# spectra_idm_30 = {'137CS': 'AK20250923_IDM_Cs137_30cm.Spe', '133BA': 'AB20250922_IDM_Ba133_30cm.Spe', '152EU': 'CD20251018_IDM_Eu152_30cm.Spe'}
# Calibration('IDM_30', spectra_idm_30, 'calibration_idm_30')
# spectra_idm_40 = {'137CS': 'AI20250923_IDM_Cs137_40cm.Spe', '133BA': 'BO20251017_IDM_Ba133_40cm.Spe', '152EU': 'CF20251018_IDM_Eu152_40cm.Spe'}
# Calibration('IDM_40', spectra_idm_40, 'calibration_idm_40')
# spectra_idm_45 = {'137CS': 'BS20251017_IDM_Cs137_45cm.Spe', '133BA': 'BN20251017_IDM_Ba133_45cm.Spe', '152EU': 'CO20251028_IDM_Eu152_45cm.Spe'}
# Calibration('IDM_45', spectra_idm_45, 'calibration_idm_45')
# spectra_idm_52 = {'137CS': 'BT20251017_IDM_Cs137_52cm.Spe', '133BA': 'BK20251017_IDM_Ba133_52cm.Spe', '152EU': 'CH20251020_IDM_Eu152_52cm.Spe'}
# Calibration('IDM_52', spectra_idm_52, 'calibration_idm_52')

# LEPS:
# spectra_leps_5 = {'241AM': 'AH20251016_LEPS_Am241_5cm.Spe'}
# Calibration('LEPS_5', spectra_leps_5, 'calibration_leps_5')
# spectra_leps_9 = {'241AM': 'DC20251030_LEPS_Am241_9cm.Spe', '133BA': 'DP20251103_LEPS_Ba133_9cm.Spe', '57CO': 'DO20251103_LEPS_Co57_9cm.Spe'}
# Calibration('LEPS_9', spectra_leps_9, 'calibration_leps_9')
# spectra_leps_10 = {'241AM': 'DB20251030_LEPS_Am241_10cm.Spe', '133BA': 'DA20251030_LEPS_Ba133_10cm.Spe', '57CO': 'DI20251030_LEPS_Co57_10cm.Spe'} #, '109CD': 'DQ20251103_LEPS_Cd109_10cm.Spe'}
# Calibration('LEPS_10', spectra_leps_10, 'calibration_leps_10')
# spectra_leps_15 = {'241AM': 'DD20251030_LEPS_Am241_15cm.Spe', '133BA': 'CX20251029_LEPS_Ba133_15cm.Spe', '57CO': 'DN20251103_LEPS_Co57_15cm.Spe'}
# Calibration('LEPS_15', spectra_leps_15, 'calibration_leps_15')
# spectra_leps_20 = {'241AM': 'DE20251030_LEPS_Am241_20cm.Spe', '133BA': 'CV20251029_LEPS_Ba133_20cm.Spe', '57CO': 'DL20251031_LEPS_Co57_20cm.Spe', '109CD': 'DS20251104_LEPS_Cd109_20cm.Spe'}
# Calibration('LEPS_20', spectra_leps_20, 'calibration_leps_20')
# spectra_leps_30 = {'241AM': 'DF20251030_LEPS_Am241_30cm.Spe', '133BA': 'AM20250923_LEPS_Ba133_30cm.Spe', '57CO': 'DK20251031_LEPS_Co57_30cm.Spe'} #'133BA': 'CZ20251030_LEPS_Ba133_30cm.Spe'}
# Calibration('LEPS_30', spectra_leps_30, 'calibration_leps_30')
# spectra_leps_40 = {'241AM': 'DH20251030_LEPS_Am241_40cm.Spe', '133BA': 'CT20251029_LEPS_Ba133_40cm.Spe', '57CO': 'DM20251103_LEPS_Co57_40cm.Spe', '109CD': 'DR20251104_LEPS_Cd109_40cm.Spe'}
# Calibration('LEPS_40', spectra_leps_40, 'calibration_leps_40')
# spectra_leps_60 = {'241AM': 'DG20251030_LEPS_Am241_60cm.Spe', '133BA': 'CP20251028_LEPS_Ba133_60cm.Spe', '57CO': 'DJ20251031_LEPS_Co57_60cm.Spe'}
# Calibration('LEPS_60', spectra_leps_60, 'calibration_leps_60')