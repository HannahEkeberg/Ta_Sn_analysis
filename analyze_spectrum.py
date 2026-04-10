from analyze_spectra import *
            
def assemble(detector, distance, foil, calibrationfile, listOfIsotopes=None):
    print('Analyzing detector ' + detector + ' at distance ' + distance + ' foil: ' + foil)
    if listOfIsotopes:
        listOfIsotopes
    else:
        listOfIsotopes = getListOfIsotopesPerFoil(foil)
    # print(listOfIsotopes)
    analyze_spec = AnalyzeSpectrum(detector=detector, calibrationFile=calibrationfile)
    spectra =  get_spectra(detector, distance, foil)
    if not spectra:
        print('No files for for detector: ' + detector + ', foil: ' + foil + ' and distance: ' + distance )
    for spec in spectra:
        if type(spec) is str:
            analyze_spec.analyze(spec, listOfIsotopes, x_rays=True)
        elif type(spec) is list:
            analyze_spec.analyze_jobs(spec, peakSummaryFilename=None, listOfIsotopes=listOfIsotopes)
        else:
            raise ValueError('Unexpected spec type: ' + spec)


def faster(foil, listOfIsotopes=None):
    assemble('Det2', '10cm', foil, 'calibration_det2_10.json', listOfIsotopes)
    assemble('Det2', '18cm', foil, 'calibration_det2_18.json', listOfIsotopes)
    assemble('Det2', '24cm', foil, 'calibration_det2_24.json', listOfIsotopes)
    assemble('Det2', '30cm', foil, 'calibration_det2_30.json', listOfIsotopes)
    assemble('Det2', '40cm', foil, 'calibration_det2_40.json', listOfIsotopes)
    assemble('Det2', '50cm', foil, 'calibration_det2_50.json', listOfIsotopes)
    assemble('Det2', '60cm', foil, 'calibration_det2_60.json', listOfIsotopes)
    assemble('Det2', '70cm', foil, 'calibration_det2_70.json', listOfIsotopes)
    assemble('Det2', '80cm', foil, 'calibration_det2_80.json', listOfIsotopes)
    assemble('IDM', '10cm', foil, 'calibration_idm_10.json', listOfIsotopes)
    assemble('IDM', '15cm', foil, 'calibration_idm_15.json', listOfIsotopes)
    assemble('IDM', '20cm', foil, 'calibration_idm_20.json', listOfIsotopes)
    assemble('IDM', '25cm', foil, 'calibration_idm_25.json', listOfIsotopes)
    assemble('IDM', '30cm', foil, 'calibration_idm_30.json', listOfIsotopes)
    assemble('IDM', '40cm', foil, 'calibration_idm_40.json', listOfIsotopes)
    assemble('IDM', '45cm', foil, 'calibration_idm_45.json', listOfIsotopes)
    assemble('IDM', '52cm', foil, 'calibration_idm_52.json', listOfIsotopes)
    assemble('LEPS', '9cm', foil, 'calibration_leps_9.json', listOfIsotopes)
    assemble('LEPS', '10cm', foil, 'calibration_leps_10.json', listOfIsotopes)
    assemble('LEPS', '15cm', foil, 'calibration_leps_15.json', listOfIsotopes)
    assemble('LEPS', '20cm', foil, 'calibration_leps_20.json', listOfIsotopes)
    assemble('LEPS', '25cm', foil, 'calibration_leps_25.json', listOfIsotopes)
    assemble('LEPS', '30cm', foil, 'calibration_leps_30.json', listOfIsotopes)
    assemble('LEPS', '40cm', foil, 'calibration_leps_40.json', listOfIsotopes)
    assemble('LEPS', '60cm', foil, 'calibration_leps_60.json', listOfIsotopes)

# isotopes = ['172HFg', '173HFg', '175HFg', '177HFm2', '178HFm2', '179HFm2', '180HFm1',
#             '179LUg', '178LUg', '178LUm1', '177LUg', '177LUm1', '176LUm1', '174LUg', '174LUm1', '173LU', '172LUg',
#             '181Wg', '179Wg', '179Wm', '178Wg', '177Wg', '176Wg', '175Wg', '174Wg', '173Wg',
#             '180TAg', '179TAg', '178TAg', '178TAm2', '178TAm1', '177TAg', '176TAg', '175TAg', '174TAg', '173TAg', '172TAg', '171TAg'
#             ]


# isotopes = ['178Wg', '178TAg', '178TAm']
# faster('Ta01', listOfIsotopes=isotopes)
# faster('Ta02', listOfIsotopes=isotopes)
# faster('Ta03', listOfIsotopes=isotopes)
# faster('Ta04', listOfIsotopes=isotopes)
# faster('Ta05', listOfIsotopes=isotopes)
# faster('Ta06', listOfIsotopes=isotopes)
# faster('Ta07', listOfIsotopes=isotopes)
# faster('Ta08', listOfIsotopes=isotopes)
# faster('Ta09', listOfIsotopes=isotopes)
# faster('Ta10', listOfIsotopes=isotopes)
# faster('Ta11', listOfIsotopes=isotopes)
# faster('Ta12', listOfIsotopes=isotopes)
# faster('Ta13', listOfIsotopes=isotopes)
# faster('Ta14', listOfIsotopes=isotopes)

# isotopes = ['119SBg']
# faster('Sn01', listOfIsotopes=isotopes)
# faster('Sn02', listOfIsotopes=isotopes)
# faster('Sn03', listOfIsotopes=isotopes)
# faster('Sn04', listOfIsotopes=isotopes)
# faster('Sn05', listOfIsotopes=isotopes)
# faster('Sn06', listOfIsotopes=isotopes)
# faster('Sn07', listOfIsotopes=isotopes)
# faster('Sn08', listOfIsotopes=isotopes)
# faster('Sn09', listOfIsotopes=isotopes)
# faster('Sn10', listOfIsotopes=isotopes)
# faster('Sn11', listOfIsotopes=isotopes)
# faster('Sn12', listOfIsotopes=isotopes)
# faster('Sn13', listOfIsotopes=isotopes)
# faster('Sn14', listOfIsotopes=isotopes)





# for f in ['Cu01','Cu02','Cu03','Cu04','Cu05','Cu06','Cu07','Cu08','Cu09','Cu10','Cu11','Cu12','Cu13','Cu14']:
# for f in ['Ni01','Ni02','Ni03','Ni04','Ni05','Ni06','Ni07','Ni08','Ni09','Ni10','Ni11','Ni12','Ni13','Ni14']:
# for f in ['Sn01','Sn02','Sn03','Sn04','Sn05','Sn06','Sn07','Sn08','Sn09','Sn10','Sn11','Sn12','Sn13','Sn14']:
#     print(f)
#     print(getListOfIsotopesPerFoil(f))
#     print("****")



# listOfIsotopes = getListOfIsotopesPerFoil('Ta01')
# print(listOfIsotopes)

# faster('Ta02')
# faster('Ta03')
# faster('Ta04')
# faster('Ta05')
# faster('Ta06')

# faster('Ta07')
# faster('Ta08')
# faster('Ta09')
# faster('Ta10')
# faster('Ta11')
# faster('Ta12')
# faster('Ta13')
# faster('Ta14')


# faster('Sn01')
# faster('Sn02')
# faster('Sn03')
# faster('Sn04')
# faster('Sn05')
# faster('Sn06')
# faster('Sn07')
# faster('Sn08')
# faster('Sn09')
# faster('Sn10')
# faster('Sn11')
# faster('Sn12')
# faster('Sn13')
# faster('Sn14')









