from nuclearanalysistools.findGammas import *
import pandas as pd
from analyze_spectra import * #list of isotopes
import curie as ci
from datetime import datetime, timedelta

nuclei_cupper = ['65ZN', '63ZN', '62ZN', '60ZN', 
                 '64CU', '59CU', '62CU', '61CU', '60CU',
                 '63NI', '56NI', '57NI', 
                 '61CO', '60CO', '60COm1', '58CO', '57CO','56CO', '55CO',
                 '59FE', '55FE', 
                 '56MN',
                 '54CR']
nuclei_nickel = ['64CU', '59CU', '62CU', '61CU', '60CU', '63NI', '56NI', '57NI', '61CO', '60CO', '58CO', '57CO','56CO', '55CO','59FE', '55FE', '56MN', '54MN', '52MN', '52FE']
nuclei_ta = ['181W', '179Wg', '179Wm', '178Wg', '177W', '176W', '175W', '179TA', '175TA', '176TA', '177TA', '178TA','175HF']



def check_matching_gammas(nuclei, isotope, xrays):
    print(nuclei)
    ag = AnalyzeGammas(nuclei)
    gammas = ag.findGammasSpecificIsotope(isotope, xrays=xrays)['Energy']
    # gammas = ag.findGammasSpecificIsotope(isotope)['Energy']
    # print(gammas)
    for g in gammas:
        data = ag.matchByGamma(gammaLine=g, gammaLineTolerance=1.0, minIntensity=0.01, xrays=xrays)
        print("****" + str(g) + "****")
        print(data)

def match_gamma_peak_data_line(nuclei, gamma, filename):
    file = os.getcwd() + '/generatedfiles/peakdata/data_isotope/' + filename
    ag = AnalyzeGammas(nuclei)
    print(gamma)
    if gamma:
        ag.matchByGamma(gammaLine=gamma, gammaLineTolerance=2.5, minIntensity=None, xrays=False)
    else:
        df = pd.read_csv(file)
        gammas = df['energy']
        for g in gammas:
            print("****" + str(g) + "****")
            ag.matchByGamma(gammaLine=g, gammaLineTolerance=2.5, minIntensity=None, xrays=False)
            print("****")



def peak_data_filter(foil, peakdatafile, nuclei=None, xrays=False, getOnlyGammasUsed=False):
    if nuclei:
        nuclei=nuclei
    else:
        nuclei = getListOfIsotopesPerFoil(foil)
    # nuclei.append('178Wg')
    # nuclei.append('178TAm1')
    file = os.getcwd() + '/generatedfiles/peakdata/data_isotope/' + peakdatafile
    df = pd.read_csv(file).sort_values('intensity', ascending=False)
    df = add_delay_time(df, foil)
    ag = AnalyzeGammas(nuclei)
    independent_gammas = []
    non_independent_gammas = []
    filtered_dfs = []
    if not getOnlyGammasUsed:
        gammas = []
        for i, row in df.iterrows():
            g = row['energy']
            I = row['intensity']
            if not g in gammas:
                data = ag.matchByGamma(gammaLine=g, gammaLineTolerance=1.0, minIntensity=None, xrays=xrays)
                gammas.append(g)
                if len(data)==1:
                    independent_gammas.append(g)
                else:
                    non_independent_gammas.append(g)
                print("****" + str(g) + "****")
                print(data)
                df_filtered = df[df['energy']==g]
                print("--------------------")
                print(df_filtered[['filename', 'isotope', 'energy', 'intensity', 'delay_time_h']]) 
                # filtered_dfs.append(df_filtered)
                print("--------------------")

                # for _, t in data.iterrows():
                    # string = f"{t['Isotope']}: ({t['Energy']}, {t['Half life']}, {t['Intensity']}%)"
                    # strings_to_print.append(string)
                    # new_df = 
        # print(gammas_df)
        print("******")
        print("@@@@@@@ independent gammas @@@@@@@")
        print(independent_gammas)
        print("@@@@@@@ NON independent gammas @@@@@@@")
        print(non_independent_gammas)
    if getOnlyGammasUsed:
        print("------ spectra taken with delay time ---------")
        df_spectra = df[['filename', 'delay_time_h']].drop_duplicates()
        print(df_spectra)
        print("------- gammas used ---------")
        df_used = df[~df['Unnamed: 0'].astype(str).str.startswith('#')]
        print(df_used[['isotope', 'energy', 'intensity', 'delay_time_h', 'filename']])


        # print("------- gammas shared with... ---------")
        # energy = df_used['energy']
        
        # for e in strings_to_print:
        #     if 



    # Want:
    # Ta01, gamma (103), contaminants (181W, 119 d, 175Ta 10 h)
    # Ta02, gamma (103), contaminants (181W, 119 d, 175Ta 10 h)
    # Ta01, gamma (103), contaminants (181W, 119 d, 175Ta 10 h)







def peak_data_match_gammas(foil, peakdatafile, nuclei=None):
    if nuclei:
        nuclei=nuclei
    else:
        nuclei = getListOfIsotopesPerFoil(foil)
    # nuclei.append('178Wg')
    # nuclei.append('178TAm1')
    file = os.getcwd() + '/generatedfiles/peakdata/data_isotope/' + peakdatafile
    df = pd.read_csv(file).sort_values('intensity', ascending=False)
    df = add_delay_time(df, foil)
    # print(df)
    # gammas = df['energy']
    ag = AnalyzeGammas(nuclei)
    rows_to_remove_not_independent = []
    rows_to_keep_half_life = []
    row_to_remove = []
    rows_to_keep_independent = []
    intense_independent_gammas = []
    gammas_non_independent_but_decayed_out = []

    gammas_df = []
    for i, row in df.iterrows():
        g = row['energy']
        I = row['intensity']
        data = ag.matchByGamma(gammaLine=g, gammaLineTolerance=1.0, minIntensity=None, xrays=False)
        if not g in gammas_df:
            gammas_df.append(g)
            print("****" + str(g) + "****")
            print(data)
        
    #     isotope = row['isotope']
    #     try:
    #         contaminants = data[data['Isotope'] != isotope]
    #         row_numb = row['Unnamed: 0']
    #         if len(contaminants)>0:
    #             row_to_remove.append(row_numb)
    #         else:
    #             rows_to_keep_independent.append(row_numb)
    #             if not g in intense_independent_gammas:
    #                 intense_independent_gammas.append(g)

    #         if (row['delay_time'] >= 10 * contaminants['Half life (s)']).all():
    #             # print(f"Gamma can be used for {isotope}: all competing isotopes have decayed out.")
    #             rows_to_keep_half_life.append(row_numb)
    #             if not g in gammas_non_independent_but_decayed_out and not g in intense_independent_gammas:
    #                 gammas_non_independent_but_decayed_out.append(g)
    #         else:
    #             # print(f"Gamma is still ambiguous for {isotope}: some competing isotopes remain.")
    #             row_to_remove.append(row_numb)
    #     except:
    #         pass
    
    # print("Remove row")
    # print(list(set(row_to_remove)))
    # print("Rows that are not independent, but contaminant have decayed")
    # print(list(set(rows_to_keep_half_life)))
    # print("Rows that are independent")
    # print(list(set(rows_to_keep_independent)))
    # print("Independent gammas highest intensity")
    # print(intense_independent_gammas)
    # print("Non independent gammas but not fed")
    # print(gammas_non_independent_but_decayed_out)
    # TODO print all gammas used (either independent, or decayed out gammas)
    return rows_to_keep_independent, rows_to_keep_half_life

def add_delay_time(data, foil):
    eob = pd.to_datetime(getEob(foil))
    data['start_time'] = pd.to_datetime(data['start_time'])
    data['delay_time'] = (data['start_time'] - eob).dt.total_seconds()
    data['delay_time_h'] = (data['start_time'] - eob).dt.total_seconds() / 3600
    return data
    
def getEob(foil):
    element = foil[0:2]
    stack1_numbs = ['01', '02', '03', '04', '05', '06', '07']
    stack2_numbs = ['08', '09','10', '11', '12', '13', '14']

    if foil in [element + number for number in stack1_numbs]:
        return '09/24/2025 15:45:32'
    elif foil in [element + number for number in stack2_numbs]:
        return '09/23/2025 18:40:05'
    else:
        raise ValueError("No valid foil: " + foil)  
    
def find_time_where_decayed(foil, isotope):
    half_life = ci.Isotope(isotope).half_life()
    eob = getEob(foil)
        # Convert string → datetime
    eob_dt = datetime.strptime(eob, "%m/%d/%Y %H:%M:%S")

    # Add 10 half-lives
    decay_time = eob_dt + timedelta(seconds=10 * half_life)

    return decay_time

# print(find_time_where_decayed('Ta01', '177W'))

# peak_data_match_gammas_for_EC(foil='Ta01', isotope='177TAg', nuclei=None)


def match(foil, isotope, xrays=False):
    # peak_data_file = foil + '_' + isotope + '_gammas.csv'
    peak_data_file = foil + '_' + isotope + '_gammas.csv'
    print("Checking peak data file:")
    print(peak_data_file)
    listOfIsotopes = getListOfIsotopesPerFoil(foil)
    # listOfIsotopes.append('178Wg')
    # listOfIsotopes.append('178HFg')
    print(listOfIsotopes)
    peak_data_filter(foil, peak_data_file, nuclei=None, xrays=xrays, getOnlyGammasUsed=False)
    peak_data_filter(foil, peak_data_file, nuclei=None, xrays=xrays, getOnlyGammasUsed=True)


isotopes = ['172HFg', '173HFg', '175HFg', '177HFm2', '178HFm2', '179HFm2', '180HFm1','180HFm2',
            '179LUg', '178LUg', '178LUm1', '177LUg', '177LUm1', '176LUm1', '174LUg', '174LUm1', '173LU', '172LUg',
            '181Wg', '179Wg', '179Wm', '178Wg', '177Wg', '176Wg', '175Wg', '174Wg', '173Wg',
            '180TAg', '179TAg', '178TAg', '178TAm2', '178TAm1', '177TAg', '176TAg', '175TAg', '174TAg', '173TAg', '172TAg', '171TAg'
            ]

def foils(element):
    stack_numbs = ['01', '02', '03', '04', '05', '06', '07', '08', '09','10', '11', '12', '13', '14']
    foils = [element + number for number in stack_numbs]
    return foils


def run_single(foil, isotope, xrays=True):
    foil = foil
    match(foil, isotope, xrays=xrays)

    return foils
def run_all(isotope, xrays=True):
    foils_ = foils('Ta')
    for f in foils_:
        print(f)
        peak_data_file = f + '_' + isotope + '_gammas.csv'
        # print(peak_data_file)
    
        listOfIsotopes = getListOfIsotopesPerFoil(f)
        try:
            peak_data_filter(f, peak_data_file, nuclei=listOfIsotopes, xrays=xrays, getOnlyGammasUsed=True)
        except:
            pass
        # # listOfIsotopes.append('188Wg')
        # # print("#####")
        # # print(listOfIsotopes)
        # # print("#####")
        # ag = AnalyzeGammas(listOfIsotopes)
        # try:
        #     # match(f, isotope, xrays=True)
        #     peak_data_filter(f, peak_data_file, nuclei=None, xrays=xrays, getOnlyGammasUsed=True)
        # except:
        #     print(f + ' did not have any data')

# run_all('181Wg')
# run_single('Ta01', '178HFm2', xrays=True)
# run_single('Ta01', '177HFm2', xrays=True)
run_single('Ta01', '179LUg', xrays=True)
# run_single('Sn01', '109CD', xrays=True)

# run_all('175TAg')
#


# check_matching_gammas(getListOfIsotopesPerFoil('Sn01'), isotope='123SNg', xrays=True)

# print(getListOfIsotopesPerFoil('Ta01'))

# try:
#     run_single('Ta01', '179TAg', xrays=True)
#     run_single('Ta03', '179TAg', xrays=True)
#     run_single('Ta04', '179TAg', xrays=True)
#     run_single('Ta05', '179TAg', xrays=True)
#     run_single('Ta06', '179TAg', xrays=True)
#     run_single('Ta07', '179TAg', xrays=True)
#     run_single('Ta08', '179TAg', xrays=True)
#     run_single('Ta09', '179TAg', xrays=True)
# except:
#     pass

# # print(ci.Isotope('178TAg').half_life()/60)
# print(ci.Isotope('178TAm').half_life()/3600)



# data = AnalyzeGammas(getListOfIsotopesPerFoil('Sn01')).findGammasSpecificIsotope('109CD', xrays=True, minIntensity=1.0)
# print(data)
# data = AnalyzeGammas(getListOfIsotopesPerFoil('Ta01')).matchByGamma(426.4)
# energy = data['Energy'].values
# # energy = [332.274]
# for e in energy:
#     d = AnalyzeGammas(getListOfIsotopesPerFoil('Ta01')).matchByGamma(e)
#     print("::: " + str(e) + " :::")
#     print(d)

# print(AnalyzeGammas(getListOfIsotopesPerFoil('Ta01')).matchByGamma(332.274))




# check_matching_gammas(['181Wg', '178Wg'], '178Wg', xrays=True)