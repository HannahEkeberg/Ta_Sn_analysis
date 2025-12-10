from analyze_spectra import *
            
def assemble(detector, distance, foil, calibrationfile):
    listOfIsotopes = getListOfIsotopesPerFoil(foil)
    analyze_spec = AnalyzeSpectrum(detector=detector, calibrationFile=calibrationfile)
    spectra =  get_spectra(detector, distance, foil)
    if not spectra:
        print('No files for for detector: ' + detector + ', foil: ' + foil + ' and distance: ' + distance )
    for spec in spectra:
        if type(spec) is str:
            analyze_spec.analyze(spec, listOfIsotopes)
        elif type(spec) is list:
            analyze_spec.analyze_jobs(spec, peakSummaryFilename=None, listOfIsotopes=listOfIsotopes)
        else:
            raise ValueError('Unexpected spec type: ' + spec)


def faster(foil):
    assemble('Det2', '10cm', foil, 'calibration_det2_10.json')
    assemble('Det2', '18cm', foil, 'calibration_det2_18.json')
    assemble('Det2', '24cm', foil, 'calibration_det2_24.json')
    assemble('Det2', '30cm', foil, 'calibration_det2_30.json')
    assemble('Det2', '40cm', foil, 'calibration_det2_40.json')
    assemble('Det2', '50cm', foil, 'calibration_det2_50.json')
    assemble('Det2', '60cm', foil, 'calibration_det2_60.json')
    assemble('Det2', '70cm', foil, 'calibration_det2_70.json')
    assemble('Det2', '80cm', foil, 'calibration_det2_80.json')
    assemble('IDM', '10cm', foil, 'calibration_idm_10.json')
    assemble('IDM', '15cm', foil, 'calibration_idm_15.json')
    assemble('IDM', '20cm', foil, 'calibration_idm_20.json')
    assemble('IDM', '25cm', foil, 'calibration_idm_25.json')
    assemble('IDM', '30cm', foil, 'calibration_idm_30.json')
    assemble('IDM', '40cm', foil, 'calibration_idm_40.json')
    assemble('IDM', '45cm', foil, 'calibration_idm_45.json')
    assemble('IDM', '52cm', foil, 'calibration_idm_52.json')
    assemble('LEPS', '9cm', foil, 'calibration_leps_9.json')
    assemble('LEPS', '10cm', foil, 'calibration_leps_10.json')
    assemble('LEPS', '15cm', foil, 'calibration_leps_15.json')
    assemble('LEPS', '20cm', foil, 'calibration_leps_20.json')
    assemble('LEPS', '25cm', foil, 'calibration_leps_25.json')
    assemble('LEPS', '30cm', foil, 'calibration_leps_30.json')
    assemble('LEPS', '40cm', foil, 'calibration_leps_40.json')
    assemble('LEPS', '60cm', foil, 'calibration_leps_60.json')


faster('Cu14')







