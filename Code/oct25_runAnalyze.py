from analyze_spectra import *



def list_jobs(prefix, detector):
    pathToSpectra = os.getcwd() + '/' + detector
    # path = "/path/to/your/folder"
    # files = [f for f in os.listdir(pathToSpectra) if f.startswith(prefix) and os.path.isfile(os.path.join(pathToSpectra, f))]
    files = [
    f for f in os.listdir(pathToSpectra)
    if f.startswith(prefix) and f.endswith(".Spe")]
    return files


# JOBS

# fm_idm_Ta12_jobs = list_jobs('FM', 'IDM')
# # AnalyzeSpectrum(detector='IDM', calibrationFile='eng_cb_idm.json', spectrumFileName=None, listOfIsotopes=getListOfIsotopesPerFoil('Ta12')).analyze_jobs(listOfJobSpectra=fm_idm_Ta12_jobs, peakSummaryFilename='FM_IDM_tot')
# fx_idm_Ta03_jobs = list_jobs('FX', 'IDM')
# # AnalyzeSpectrum(detector='IDM', calibrationFile='eng_cb_idm.json', spectrumFileName=None, listOfIsotopes=getListOfIsotopesPerFoil('Ta12')).analyze_jobs(listOfJobSpectra=fm_idm_Ta12_jobs, peakSummaryFilename='FM_IDM_tot')

# fy_idm_Ta04_jobs = list_jobs('FY', 'IDM')
# gd_idm_Ta08_jobs = list_jobs('GD', 'IDM')
# ge_idm_Ta09_jobs = list_jobs('GE', 'IDM')
# gf_idm_Ta10_jobs = list_jobs('GF', 'IDM')
# gg_idm_Ta11_jobs = list_jobs('GG', 'IDM')
# gh_idm_Ta12_jobs = list_jobs('GH', 'IDM')
# gi_idm_Ta13_jobs = list_jobs('GI', 'IDM')
# gj_idm_Ta14_jobs = list_jobs('GJ', 'IDM')
# gk_idm_Cu13_jobs = list_jobs('GK', 'IDM')
# gl_idm_Cu14_jobs = list_jobs('GL', 'IDM')
# gm_idm_Cu12_jobs = list_jobs('GM', 'IDM')
# gn_idm_Cu10_jobs = list_jobs('GN', 'IDM')
# go_idm_Cu08_jobs = list_jobs('GO', 'IDM')
# gp_idm_Cu04_jobs = list_jobs('GP', 'IDM')
# gq_idm_Cu02_jobs = list_jobs('GQ', 'IDM')


# fb_det2_Ta08_jobs = list_jobs('FB', 'Det2')
# fj_det2_Ta11_jobs = list_jobs('FJ', 'Det2')
# fp_det2_Ta01_jobs = list_jobs('FP', 'Det2')
# ft_det2_Sn02_jobs = list_jobs('FT', 'Det2')
# fu_det2_Sn03_jobs = list_jobs('FU', 'Det2')
# gc_det2_Sn11_jobs = list_jobs('GC', 'Det2')
# gd_det2_Sn12_jobs = list_jobs('GD', 'Det2')
# ge_det2_Sn13_jobs = list_jobs('GE', 'Det2')
# gf_det2_Ni14_jobs = list_jobs('GF', 'Det2')
# gg_det2_Ni13_jobs = list_jobs('GG', 'Det2')
# gh_det2_Ni12_jobs = list_jobs('GH', 'Det2')
# gi_det2_Ni11_jobs = list_jobs('GI', 'Det2')
# gj_det2_Sn14_jobs = list_jobs('GJ', 'Det2')
# gk_det2_Ni10_jobs = list_jobs('GK', 'Det2')
# gl_det2_Ni09_jobs = list_jobs('GL', 'Det2')
# gm_det2_Ni08_jobs = list_jobs('GM', 'Det2')
# go_det2_Ni05_jobs = list_jobs('GO', 'Det2')
# gp_det2_Ni06_jobs = list_jobs('GP', 'Det2')
# gq_det2_Ni04_jobs = list_jobs('GQ', 'Det2')
# gr_det2_Ni03_jobs = list_jobs('GR', 'Det2')
# gs_det2_Ni02_jobs = list_jobs('GS', 'Det2')
# gt_det2_Ni01_jobs = list_jobs('GT', 'Det2')
# gu_det2_Cu11_jobs = list_jobs('GU', 'Det2')
# gv_det2_Cu09_jobs = list_jobs('GV', 'Det2')
# gw_det2_Cu07_jobs = list_jobs('GW', 'Det2')
# gx_det2_Cu06_jobs = list_jobs('GX', 'Det2')
# gy_det2_Cu05_jobs = list_jobs('GY', 'Det2')


# ez_LEPS_Ta10_jobs = list_jobs('EZ', 'LEPS')
# fg_LEPS_Ta13_jobs = list_jobs('FG', 'LEPS')
# fk_LEPS_Ta03_jobs = list_jobs('FK', 'LEPS')


# jobs_idm = [fm_idm_Ta12_jobs, fx_idm_Ta03_jobs, fy_idm_Ta04_jobs, gd_idm_Ta08_jobs, 
#         ge_idm_Ta09_jobs, gf_idm_Ta10_jobs, gg_idm_Ta11_jobs, gh_idm_Ta12_jobs,gi_idm_Ta13_jobs,gj_idm_Ta14_jobs,
#         gk_idm_Cu13_jobs, gl_idm_Cu14_jobs, gm_idm_Cu12_jobs,gn_idm_Cu10_jobs, go_idm_Cu08_jobs,gp_idm_Cu04_jobs,gq_idm_Cu02_jobs]

# jobs_det2 = [fb_det2_Ta08_jobs,fj_det2_Ta11_jobs,fp_det2_Ta01_jobs,ft_det2_Sn02_jobs,fu_det2_Sn03_jobs,gc_det2_Sn11_jobs,gd_det2_Sn12_jobs,ge_det2_Sn13_jobs,gf_det2_Ni14_jobs,gg_det2_Ni13_jobs,gh_det2_Ni12_jobs,gi_det2_Ni11_jobs,gj_det2_Sn14_jobs,gk_det2_Ni10_jobs,gl_det2_Ni09_jobs,gm_det2_Ni08_jobs,go_det2_Ni05_jobs,gp_det2_Ni06_jobs,gq_det2_Ni04_jobs,gr_det2_Ni03_jobs,gs_det2_Ni02_jobs,gt_det2_Ni01_jobs,gu_det2_Cu11_jobs,gv_det2_Cu09_jobs,gw_det2_Cu07_jobs,gx_det2_Cu06_jobs,gy_det2_Cu05_jobs]

# jobs_leps = [ez_LEPS_Ta10_jobs,fg_LEPS_Ta13_jobs,fk_LEPS_Ta03_jobs]

# for j in jobs_idm:
#     parts = j[0].split("_")
#     prefix = ''.join([c for c in parts[0] if c.isalpha()])
#     isotope = parts[2]
#     AnalyzeSpectrum(detector='IDM', calibrationFile='eng_cb_idm.json', spectrumFileName=None, listOfIsotopes=getListOfIsotopesPerFoil(isotope)).analyze_jobs(listOfJobSpectra=j, peakSummaryFilename=prefix + '_' + isotope + '_IDM_tot')

# for j in jobs_det2:
#     parts = j[0].split("_")
#     prefix = ''.join([c for c in parts[0] if c.isalpha()])
#     isotope = parts[2]
#     AnalyzeSpectrum(detector='DET2', calibrationFile='eng_cb_det2.json', spectrumFileName=None, listOfIsotopes=getListOfIsotopesPerFoil(isotope)).analyze_jobs(listOfJobSpectra=j, peakSummaryFilename=prefix + '_' + isotope +  '_DET2_tot')

# for j in jobs_leps:
#     parts = j[0].split("_")
#     prefix = ''.join([c for c in parts[0] if c.isalpha()])
#     isotope = parts[2]
#     AnalyzeSpectrum(detector='LEPS', calibrationFile='eng_cb_leps.json', spectrumFileName=None, listOfIsotopes=getListOfIsotopesPerFoil(isotope)).analyze_jobs(listOfJobSpectra=j, peakSummaryFilename=prefix + '_' + isotope + '_LEPS_tot')

# All:

# idm_spes = os.getcwd() + '/IDM/idm_specs.txt' 
# det2_spes = os.getcwd() + '/DET2/det2_specs.txt'
# leps_spes = os.getcwd() + '/LEPS/leps_specs.txt'


# idm_filenames = np.genfromtxt(idm_spes, dtype=str, delimiter="\n")
# det2_filenames = np.genfromtxt(det2_spes, dtype=str, delimiter="\n")
# leps_filenames = np.genfromtxt(leps_spes, dtype=str, delimiter="\n")

# print(idm_filenames)



def getListOfIsotopesPerFoil(foil, giveActivity=False):
    uniqueIsotopes_55 = pd.read_csv('isotopes_from_simulation/Unique_isotopes_55MeV.csv')
    uniqueIsotopes_30 = pd.read_csv('isotopes_from_simulation/Unique_isotopes_30MeV.csv')
    uniqueIsotopes = pd.concat([uniqueIsotopes_55, uniqueIsotopes_30])
    uniqueIsotopes = uniqueIsotopes[uniqueIsotopes['Activity foil'].astype(str).str.contains(foil, case=False, na=False)]
    if giveActivity:
        return uniqueIsotopes
    return list(uniqueIsotopes['isotope'])

# analyzeSeveralSpectra(det2_filenames, 'DET2', 'eng_cb_det2.json', 'Ni')
# analyzeSeveralSpectra(idm_filenames, 'IDM', 'eng_cb_idm.json', 'Ni')
# analyzeSeveralSpectra(leps_filenames, 'LEPS', 'eng_cb_leps.json', 'Ni')


# analyzeSeveralSpectra(det2_filenames, 'DET2', 'eng_cb_det2.json', 'Cu')
# analyzeSeveralSpectra(idm_filenames, 'IDM', 'eng_cb_idm.json', 'Cu')
# analyzeSeveralSpectra(leps_filenames, 'LEPS', 'eng_cb_leps.json', 'Cu')

# analyzeSeveralSpectra(det2_filenames, 'DET2', 'eng_cb_det2.json', 'Ta')
# analyzeSeveralSpectra(idm_filenames, 'IDM', 'eng_cb_idm.json', 'Ta')
# analyzeSeveralSpectra(leps_filenames, 'LEPS', 'eng_cb_leps.json', 'Ta')

# analyzeSeveralSpectra(det2_filenames, 'DET2', 'eng_cb_det2.json', 'Sn')
# analyzeSeveralSpectra(idm_filenames, 'IDM', 'eng_cb_idm.json', 'Sn')
# analyzeSeveralSpectra(leps_filenames, 'LEPS', 'eng_cb_leps.json', 'Sn')



# Get peak summaries in each foil

# peakDataFiles = os.getcwd() + '/generatedFiles/peakData/peakfiles.txt'
# peakDataFileNames = np.genfromtxt(peakDataFiles, dtype=str, delimiter="\n")
# print(peakDataFileNames)

# generatePeakOverviews(element='Sn', count_filenames=peakDataFileNames)




# Check calibration sources


# - sjekke om det har vært gain shift: 80 cm 152Eu kalibrering. 
# - Co57 x-ray detektor --> er det gode nok stats? 
# - Ba133 70 cm detektor 2 --> gode nok stats? 


#  AnalyzeSpectrum(detector='DET2', calibrationFile, spectrumFileName, listOfIsotopes)
AnalyzeSpectrum(detector='DET2', calibrationFile='eng_cb_det2.json', spectrumFileName='../Calibration/CN20251020_Det2_Eu152_80cm.Spe', listOfIsotopes=['152EU']).analyze()
AnalyzeSpectrum(detector='DET2', calibrationFile='eng_cb_det2.json', spectrumFileName='../Calibration/CG20251018_Det2_Ba133_70cm.Spe', listOfIsotopes=['133BA']).analyze()
# AnalyzeSpectrum(detector='DET2', calibrationFile='eng_cb_det2.json', spectrumFileName='../Calibration/CG20251018_Det2_Ba133_70cm.Spe', listOfIsotopes=['133BA']).analyze()



AnalyzeSpectrum(detector='LEPS', calibrationFile='eng_cb_leps.json', spectrumFileName='../Calibration/CN20251020_Det2_Eu152_80cm.Spe', listOfIsotopes=['152EU']).analyze()