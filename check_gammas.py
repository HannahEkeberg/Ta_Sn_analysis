import sys
sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
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
                 '56MN', '54MN',
                 '54CR']
nuclei_nickel = ['64CU', '59CU', '62CU', '61CU', '60CU', '63NI', '56NI', '57NI', '61CO', '60CO', '58CO', '57CO','56CO', '55CO','59FE', '55FE', '56MN', '54MN', '52MN', '52FE', '51CR']
nuclei_ta = ['181W', '179Wg', '179Wm', '178Wg', '177W', '176W', '175W', '179TA', '175TA', '176TA', '177TA', '178TA','175HF']



def check_matching_gammas(nuclei, isotope, xrays):
    print(nuclei)
    ag = AnalyzeGammas(nuclei)
    gammas = ag.findGammasSpecificIsotope(isotope, xrays=xrays)['Energy']
    # gammas = ag.findGammasSpecificIsotope(isotope)['Energy']
    # print(gammas)
    for g in gammas:
        print("****" + str(g) + "****")
        data = ag.matchByGamma(gammaLine=g, gammaLineTolerance=1.5, minIntensity=1.0, xrays=xrays)
        print(data)
        print("****" + "****")

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
    # nuclei.append('179Wg')
    # nuclei.append('178Wg')
    # nuclei.append('178TAm1')
    file = os.getcwd() + '/generatedfiles/peakdata/data_isotope/' + peakdatafile
    df = pd.read_csv(file).sort_values('intensity', ascending=False)
    df = add_delay_time(df, foil)
    df['filename'] = df['filename'].str.split('/').str[-1] # fjern om ikke jobber med tinn
    # print(df)
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



        # UNCOMMENT LATER
        print("******")
        print("@@@@@@@ independent gammas @@@@@@@")
        print(independent_gammas)
        print("@@@@@@@ NON independent gammas @@@@@@@")
        print(non_independent_gammas)
    if getOnlyGammasUsed:
        # print("------ spectra taken with delay time ---------")
        # df_spectra = df[['filename', 'delay_time_h']].drop_duplicates()
        # print(df_spectra)
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



def get_all_spectra():
    list_of_files =  os.getcwd() + '/generatedfiles/peakdata/data'
    filenames = [f for f in os.listdir(list_of_files) if f.endswith(".csv")]
    for file in filenames:
        # Filnavn
        filename = os.path.basename(file)
        foil = os.path.basename(os.path.dirname(file))
        print(filename, foil)


# Parent folder (Det2)
        df = pd.read_csv(file).sort_values('intensity', ascending=False)
        df = add_delay_time(df, foil)

    # print(csv_files)
    # file = os.getcwd() + '/generatedfiles/peakdata/data' + peakdatafile
    # df = pd.read_csv(file).sort_values('intensity', ascending=False)
    # df = add_delay_time(df, foil)


# get_all_spectra()
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

def comment_out_gammas(foil, isotope, gamma):
    peak_data_file = foil + '_' + isotope + '_gammas.csv'
    file = os.getcwd() + '/generatedfiles/peakdata/data_isotope/' + peak_data_file
    df = pd.read_csv(file)
    for index, line in df.iterrows():
        if line['energy'] == gamma:
            print('Comment out gamma line ' + str(gamma))
            with open(file, 'r') as f:
                lines = f.readlines()
            # +1 because CSV line 0 is the header
            if not lines[index + 1].startswith('#'):
                lines[index + 1] = '#' + lines[index + 1]
            # lines[index + 1] = '#' + lines[index + 1]
            with open(file, 'w') as f:
                f.writelines(lines)
        # else: 
        #     print('no lines with gamma energy ' + str(gamma))

def comment_in_gammas(foil, isotope, gamma):
    peak_data_file = foil + '_' + isotope + '_gammas.csv'
    file = os.getcwd() + '/generatedfiles/peakdata/data_isotope/' + peak_data_file
    df = pd.read_csv(file)
    for index, line in df.iterrows():
        if line['energy'] == gamma:
            print('Comment in gamma line ' + str(gamma))
            with open(file, 'r') as f:
                lines = f.readlines()
            # +1 because CSV line 0 is the header
            lines[index + 1] = lines[index + 1].lstrip('#')
            with open(file, 'w') as f:
                f.writelines(lines)
        # else: 
        #     print('no lines with gamma energy ' + str(gamma))


def keep_selected_gammas(foil, isotope, gammas, unmark_first=False):
    peak_data_file = foil + '_' + isotope + '_gammas.csv'
    # print(gammas)
    # print(peak_data_file)
    file = os.getcwd() + '/generatedfiles/peakdata/data_isotope/' + peak_data_file
    df = pd.read_csv(file)
    marked_lines = []
    
    if unmark_first:
        print('unmarking all lines for file: ' + peak_data_file)
        unmark(file)
    if len(gammas) > 0:
        for index, line in df.iterrows():
            # print(line['energy'], str(gammas))
            # if index==0:
            #     print(line)
            #     continue
            idx = index + 1
            if line['energy'] not in gammas:
                if line['energy'] not in marked_lines:
                    marked_lines.append(line['energy'])
                with open(file, 'r') as f:
                    lines = f.readlines()
                if not lines[idx].startswith('#'):
                    lines[idx] = '#' + lines[idx]
                # lines[index + 1] = '#' + lines[index + 1]
                with open(file, 'w') as f:
                    f.writelines(lines)
    print("marked lines")
    print(marked_lines)
    print('unmarked lines')
    print(gammas)

        #     print(str(line['energy']) + 'not in gammas: ' + str(gammas)) 
        # else:
        #     print(line['energy'])
            # with open(file, 'r') as f:
            #     lines = f.readlines()
            # # +1 because CSV line 0 is the header
            # if '#' in lines[idx]:
            #     print(lines[idx])
            #     lines[idx] = lines[idx].lstrip('#')
            # #if not lines[idx].startswith('#'):
            # # lines[idx] = '#' + lines[idx]
            # # lines[0] = lines[0].lstrip('#')
            # # # lines[index + 1] = '#' + lines[index + 1]
            # with open(file, 'w') as f:
            #     f.writelines(lines)

def unmark(file):
    with open(file, 'r') as f:
        lines = f.readlines()

        for i, line in enumerate(lines):
            if '#' in line:
                lines[i] = line.lstrip('#')

    with open(file, 'w') as f:
        f.writelines(lines)

def get_all_used_gammalines(element, isotope):
    # element_foils = foils(element)
    unique_gammas_all_foils = []
    for foil in foils(element):
        print(foil)
        print('half life (min, h, d)')
        print(ci.Isotope(isotope).half_life('min'), ci.Isotope(isotope).half_life('h'), ci.Isotope(isotope).half_life('d'))
        peak_data_file = foil + '_' + isotope + '_gammas.csv'
        try:
            file = os.getcwd() + '/generatedfiles/peakdata/data_isotope/' + peak_data_file
            df = pd.read_csv(file, comment='#')
            print(df[['energy', 'intensity', 'unc_intensity']].drop_duplicates(subset='energy'))
            unique_gammas_all_foils.extend(g for g in df['energy'].unique() if g not in unique_gammas_all_foils)
            print("----------------------------")
        except:
            print("no peak data file found: " + peak_data_file)
    print("Gammas used for isotope: " + isotope)
    sorted_unique = np.sort(unique_gammas_all_foils)
    print(sorted_unique)
    return sorted_unique

def match(foil, isotope, xrays=False):
    # peak_data_file = foil + '_' + isotope + '_gammas.csv'
    peak_data_file = foil + '_' + isotope + '_gammas.csv'
    print("Checking peak data file:")
    print(peak_data_file)
    listOfIsotopes = getListOfIsotopesPerFoil(foil)
    # listOfIsotopes = ['181Wg', '180TAg', '179TAg', '180HFm', '179HFm2', '178HFm2', '176LUm']
    listOfIsotopes.append('178TAm1')
    # listOfIsotopes.append('178HFg')
    print(listOfIsotopes)
    peak_data_filter(foil, peak_data_file, nuclei=None, xrays=xrays, getOnlyGammasUsed=False)
    peak_data_filter(foil, peak_data_file, nuclei=None, xrays=xrays, getOnlyGammasUsed=True)


# isotopes = ['172HFg', '173HFg', '175HFg', '177HFm2', '178HFm2', '179HFm2', '180HFm1','180HFm2',
#             '179LUg', '178LUg', '178LUm1', '177LUg', '177LUm1', '176LUm1', '174LUg', '174LUm1', '173LU', '172LUg',
#             '181Wg', '179Wg', '179Wm', '178Wg', '177Wg', '176Wg', '175Wg', '174Wg', '173Wg',
#             '180TAg', '179TAg', '178TAg', '178TAm2', '178TAm1', '177TAg', '176TAg', '175TAg', '174TAg', '173TAg', '172TAg', '171TAg'
#             ]

def foils(element):
    stack_numbs = ['01', '02', '03', '04', '05', '06', '07', '08', '09','10', '11', '12', '13', '14']
    foils = [element + number for number in stack_numbs]
    return foils


def run_single(foil, isotope, xrays=True):
    foil = foil
    match(foil, isotope, xrays=xrays)

def run_all(isotope, element='Ta', xrays=True):
    # foils_ = foils(element)
    foils_ = foils(element)
    for f in foils_:
        print(f)
        peak_data_file = f + '_' + isotope + '_gammas.csv'
        listOfIsotopes = getListOfIsotopesPerFoil(f)

        try:
            peak_data_filter(f, peak_data_file, nuclei=listOfIsotopes, xrays=xrays, getOnlyGammasUsed=True)
        except:
            pass
        print("--------")



# run_single('Ni01', '58COg', xrays=False)



# run_single('Sn01', '109SNg', xrays=False)
# run_single('Sn01', '111SNg', xrays=False)
# run_single('Ta01', '178TAg', xrays=True)
# run_single('Ta12', '181Wg', xrays=True)
# run_single('Ta11', '181Wg', xrays=True)


# run_single('Ta01', '177Wg', xrays=False)

# comment_out_gammas(foil='Ta12', isotope='178TAm1', gamma=88.857)
# comment_out_gammas(foil='Ta12', isotope='178TAm1', gamma=426.355)
# comment_out_gammas(foil='Ta12', isotope='178TAm1', gamma=93.179)





# a = AnalyzeGammas([]).findGammasSpecificIsotope('178TAg', xrays=False)#['Energy']
# b = AnalyzeGammas([]).findGammasSpecificIsotope('178TAm', xrays=False)#['Energy']
# print(a)
# print(b)
# run_all('58COg', 'Ni', False)

# listOfIsotopes_all = [
#                         '181Wg', '179Wm', '179Wg', '178Wg', '177Wg', '176Wg', '175Wg',
#                         '180TAg', '179TAg', '178TAg', '178TAm', '177TAg', '176TAg', '165TAg', '174TAg', '173TAg',
#                         '180HFm', '179HFm2', '178HFm2', '177HFm2', '175HFg', '173HFg','172HFg',
#                         '179LUg', '178LUm', '178LUg', '177LUm', '177LUg', '176LUg', '174LUm', '174LUg', '173LUg', '171LUg'
#                     ]


# data = AnalyzeGammas(listOfIsotopes_all).findAllGammas(minIntensity=None, xrays=False)#matchByGamma(152.535, gammaLineTolerance=2.5)
# print(data)




# ag = AnalyzeGammas([]).findGammasSpecificIsotope(iso='116SBm', minIntensity=0.01, xrays=None)
# print(ag)


# isotoeps = ['157DY', '159DY', '163TB', '162TB', '161TB', '160TB', '156TB', '155TB']
# ag = AnalyzeGammas(isotoeps).findAllGammas(minIntensity=10)
# print(ag)

# xrays_178W = [56.28, 57.535, 64.948, 65.222, 66.982]
# keep_selected_gammas(foil='Ta01', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta02', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta01', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta04', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta03', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta04', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta05', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta06', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta07', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta08', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta09', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta10', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta11', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta12', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta13', isotope='178Wg', gammas=xrays_178W, unmark_first=True)
# keep_selected_gammas(foil='Ta14', isotope='178Wg', gammas=xrays_178W, unmark_first=True)



# keep_selected_gammas(foil='Ta01', isotope='177Wg', gammas=[417.16, 611.8, 1014.9, 367.52, 271.02, 377.35, 1036.4,  1066.9, 1182.5 ], unmark_first=True)
# keep_selected_gammas(foil='Ta02', isotope='177Wg', gammas=[417.16, 611.8, 1014.9, 367.52, 271.02, 377.35, 1036.4,  1066.9, 1182.5 ], unmark_first=True)
# keep_selected_gammas(foil='Ta01', isotope='177Wg', gammas=[271.02, 115.65], unmark_first=True)
# keep_selected_gammas(foil='Ta04', isotope='177Wg', gammas=[271.02, 115.65], unmark_first=True)
# keep_selected_gammas(foil='Ta03', isotope='177Wg', gammas=[271.02, 115.65], unmark_first=True)
# keep_selected_gammas(foil='Ta04', isotope='177Wg', gammas=[271.02, 115.65], unmark_first=True)
# keep_selected_gammas(foil='Ta05', isotope='177Wg', gammas=[271.02, 115.65], unmark_first=True)
# keep_selected_gammas(foil='Ta06', isotope='177Wg', gammas=[271.02, 115.65], unmark_first=True)
# keep_selected_gammas(foil='Ta07', isotope='177Wg', gammas=[271.02, 115.65], unmark_first=True)
# keep_selected_gammas(foil='Ta08', isotope='177Wg', gammas=[271.02, 115.65], unmark_first=True)

# keep_selected_gammas(foil='Ta01', isotope='177Wg', gammas=[271.02, 417.16, 186.2], unmark_first=True)
# keep_selected_gammas(foil='Ta02', isotope='177Wg', gammas=[271.02, 417.16, 186.2], unmark_first=True)
# keep_selected_gammas(foil='Ta03', isotope='177Wg', gammas=[271.02, 115.65], unmark_first=True)
# keep_selected_gammas(foil='Ta04', isotope='177Wg', gammas=[271.02, 115.65], unmark_first=True)
# keep_selected_gammas(foil='Ta05', isotope='177Wg', gammas=[271.02, 417.16, 186.2], unmark_first=True)
# keep_selected_gammas(foil='Ta06', isotope='177Wg', gammas=[271.02, 417.16, 186.2], unmark_first=True)
# keep_selected_gammas(foil='Ta07', isotope='177Wg', gammas=[271.02, 417.16, 186.2], unmark_first=True)
# keep_selected_gammas(foil='Ta08', isotope='177Wg', gammas=[271.02, 417.16, 186.2], unmark_first=True)




# 176TA
#keep_selected_gammas(foil='Ta01', isotope='176TAg', gammas=[1159.3, 201.84, 710.5, 1584.02, 1696.55, 1823.7], unmark_first=True)
#keep_selected_gammas(foil='Ta02', isotope='176TAg', gammas=[1159.3, 201.84, 710.5, 1584.02, 1696.55, 1823.7], unmark_first=True)
#keep_selected_gammas(foil='Ta03', isotope='176TAg', gammas=[1159.3, 201.84, 710.5, 1584.02, 1696.55, 1823.7], unmark_first=True)
#keep_selected_gammas(foil='Ta04', isotope='176TAg', gammas=[1159.3, 201.84, 710.5, 1584.02, 1696.55, 1823.7], unmark_first=True)
#keep_selected_gammas(foil='Ta05', isotope='176TAg', gammas=[1159.3, 201.84, 710.5, 1584.02, 1696.55, 1823.7], unmark_first=True)
#keep_selected_gammas(foil='Ta06', isotope='176TAg', gammas=[1159.3, 201.84, 710.5, 1584.02, 1696.55, 1823.7], unmark_first=True)


# 180TA 103.6,93.324
gammas_180Ta = [103.6, 93.324]
# keep_selected_gammas(foil='Ta01', isotope='180Tag', gammas=[103.6], unmark_first=True)
# keep_selected_gammas(foil='Ta01', isotope='180Tag', gammas=gammas_180Ta, unmark_first=True)
# keep_selected_gammas(foil='Ta02', isotope='180TAg', gammas=gammas_180Ta, unmark_first=True)
# keep_selected_gammas(foil='Ta03', isotope='180TAg', gammas=gammas_180Ta, unmark_first=True)
# keep_selected_gammas(foil='Ta04', isotope='180TAg', gammas=gammas_180Ta, unmark_first=True)
# keep_selected_gammas(foil='Ta05', isotope='180TAg', gammas=gammas_180Ta, unmark_first=True)
# keep_selected_gammas(foil='Ta06', isotope='180TAg', gammas=gammas_180Ta, unmark_first=True)
# keep_selected_gammas(foil='Ta07', isotope='180TAg', gammas=gammas_180Ta, unmark_first=True)
# keep_selected_gammas(foil='Ta08', isotope='180TAg', gammas=gammas_180Ta, unmark_first=True)
# keep_selected_gammas(foil='Ta09', isotope='180TAg', gammas=gammas_180Ta, unmark_first=True)
# keep_selected_gammas(foil='Ta10', isotope='180TAg', gammas=gammas_180Ta, unmark_first=True)
# keep_selected_gammas(foil='Ta11', isotope='180TAg', gammas=gammas_180Ta, unmark_first=True)
# keep_selected_gammas(foil='Ta12', isotope='180TAg', gammas=gammas_180Ta, unmark_first=True)
# keep_selected_gammas(foil='Ta13', isotope='180TAg', gammas=gammas_180Ta, unmark_first=True)
# keep_selected_gammas(foil='Ta14', isotope='180TAg', gammas=gammas_180Ta, unmark_first=True)


nuclei_ta = ['181Wg', '179Wg', '179Wm', '178Wg', '177Wg', '176Wg', '175Wg', 
          '180TAg', '179TAg', '178TAg', '178TAm1', '177TAg','176TAg', '175TAg', '174TAg', '173TAg',
          '180HFm', '179HFm2', '178HFm2', '177HFm2', '175HFg', '173HFg', 
          '179LUg', '178LUg', '177LUg', '176LUm', ]
# check_matching_gammas(nuclei_ta, '178HFm2', False)
# run_single('Ta01', '179HFm2', xrays=False)



# 500.697, 332.274, 443.162,215.426
# keep_selected_gammas(foil='Ta01', isotope='180HFm', gammas=[443.162,215.426], unmark_first=True)
# keep_selected_gammas(foil='Ta02', isotope='180HFm', gammas=[443.162,215.426], unmark_first=True)
# keep_selected_gammas(foil='Ta03', isotope='180HFm', gammas=[443.162,215.426], unmark_first=True)
# keep_selected_gammas(foil='Ta04', isotope='180HFm', gammas=[443.162,215.426], unmark_first=True)
# keep_selected_gammas(foil='Ta05', isotope='180HFm', gammas=[443.162,215.426], unmark_first=True)
# keep_selected_gammas(foil='Ta06', isotope='180HFm', gammas=[443.162,215.426], unmark_first=True)
# keep_selected_gammas(foil='Ta07', isotope='180HFm', gammas=[443.162,215.426], unmark_first=True)
# keep_selected_gammas(foil='Ta08', isotope='180HFm', gammas=[443.162,215.426], unmark_first=True)
# keep_selected_gammas(foil='Ta09', isotope='180HFm', gammas=[443.162,215.426], unmark_first=True)
# keep_selected_gammas(foil='Ta10', isotope='180HFm', gammas=[443.162,215.426], unmark_first=True)




# keep_selected_gammas(foil='Ta01', isotope='177TAg', gammas=[745.9, 944.8, 112.9, 208.4, 1057.8], unmark_first=True)
# keep_selected_gammas(foil='Ta02', isotope='177TAg', gammas=[745.9, 944.8, 112.9, 208.4, 1057.8], unmark_first=True)
# keep_selected_gammas(foil='Ta03', isotope='177TAg', gammas=[745.9, 944.8, 112.9, 208.4, 1057.8], unmark_first=True)
# keep_selected_gammas(foil='Ta04', isotope='177TAg', gammas=[745.9, 944.8, 112.9, 208.4, 1057.8], unmark_first=True)
# keep_selected_gammas(foil='Ta05', isotope='177TAg', gammas=[745.9, 944.8, 112.9, 208.4, 1057.8], unmark_first=True)
# keep_selected_gammas(foil='Ta06', isotope='177TAg', gammas=[745.9, 944.8, 112.9, 208.4, 1057.8], unmark_first=True)
# keep_selected_gammas(foil='Ta07', isotope='177TAg', gammas=[745.9, 944.8, 112.9, 208.4, 1057.8], unmark_first=True)
# keep_selected_gammas(foil='Ta08', isotope='177TAg', gammas=[745.9, 944.8, 112.9, 208.4, 1057.8], unmark_first=True)


# keep_selected_gammas(foil='Ta01', isotope='175TAg', gammas=[348.5, 436.4, 998.3, 1793.1], unmark_first=True)
# keep_selected_gammas(foil='Ta02', isotope='175TAg', gammas=[348.5, 436.4, 998.3, 1793.1], unmark_first=True)
# keep_selected_gammas(foil='Ta03', isotope='175TAg', gammas=[348.5, 436.4, 998.3, 1793.1], unmark_first=True)
# keep_selected_gammas(foil='Ta04', isotope='175TAg', gammas=[348.5, 436.4, 998.3, 1793.1], unmark_first=True)
# keep_selected_gammas(foil='Ta05', isotope='175TAg', gammas=[348.5, 436.4, 998.3, 1793.1], unmark_first=True)



# AnalyzeGammas(nuclei)

# 207.4, 266.9, 348.5,  436.4, 475.0, 998.3, 1793.1
# AnalyzeGammas(nuclei).matchByGamma(gammaLine=207.4, gammaLineTolerance=1.5, minIntensity=None, xrays=False)
# AnalyzeGammas(nuclei).matchByGamma(gammaLine=266.9, gammaLineTolerance=1.5, minIntensity=None, xrays=False)
# AnalyzeGammas(nuclei).matchByGamma(gammaLine=475.0, gammaLineTolerance=1.5, minIntensity=None, xrays=False)
# AnalyzeGammas(nuclei).matchByGamma(gammaLine=348.5, gammaLineTolerance=1.5, minIntensity=None, xrays=False)
# AnalyzeGammas(nuclei).matchByGamma(gammaLine=436.4, gammaLineTolerance=1.5, minIntensity=None, xrays=False)
# AnalyzeGammas(nuclei).matchByGamma(gammaLine=998.3, gammaLineTolerance=1.5, minIntensity=None, xrays=False)
# AnalyzeGammas(nuclei).matchByGamma(gammaLine=133.9, gammaLineTolerance=1, minIntensity=None, xrays=False)
# AnalyzeGammas(nuclei).matchByGamma(gammaLine=30.7, gammaLineTolerance=1, minIntensity=None, xrays=True)


# gammas_173Hf = [123.675, 296.974, 139.635, 306.568, 134.965]

# keep_selected_gammas(foil='Ta01', isotope='173HFg', gammas=gammas_173Hf, unmark_first=True)
# keep_selected_gammas(foil='Ta02', isotope='173HFg', gammas=gammas_173Hf, unmark_first=True)
# keep_selected_gammas(foil='Ta03', isotope='173HFg', gammas=gammas_173Hf, unmark_first=True)
# keep_selected_gammas(foil='Ta04', isotope='173HFg', gammas=gammas_173Hf, unmark_first=True)
# keep_selected_gammas(foil='Ta05', isotope='173HFg', gammas=gammas_173Hf, unmark_first=True)
# keep_selected_gammas(foil='Ta06', isotope='173HFg', gammas=gammas_173Hf, unmark_first=True)

# gammas_179Hfm2 = [362.6, 453.6]
# keep_selected_gammas(foil='Ta01', isotope='179HFm2', gammas=gammas_179Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta02', isotope='179HFm2', gammas=gammas_179Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta03', isotope='179HFm2', gammas=gammas_179Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta04', isotope='179HFm2', gammas=gammas_179Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta05', isotope='179HFm2', gammas=gammas_179Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta06', isotope='179HFm2', gammas=gammas_179Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta07', isotope='179HFm2', gammas=gammas_179Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta08', isotope='179HFm2', gammas=gammas_179Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta09', isotope='179HFm2', gammas=gammas_179Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta10', isotope='179HFm2', gammas=gammas_179Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta11', isotope='179HFm2', gammas=gammas_179Hfm2, unmark_first=True)

# gammas_177Hfm2=[606.5, 572.4, 326.7, 295.1, 277.3]
# gammas_177Hfm2=[606.5, 295.1, 277.3]
# keep_selected_gammas(foil='Ta01', isotope='177HFm2', gammas=gammas_177Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta02', isotope='177HFm2', gammas=gammas_177Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta03', isotope='177HFm2', gammas=gammas_177Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta04', isotope='177HFm2', gammas=gammas_177Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta05', isotope='177HFm2', gammas=gammas_177Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta06', isotope='177HFm2', gammas=gammas_177Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta07', isotope='177HFm2', gammas=gammas_177Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta08', isotope='177HFm2', gammas=gammas_177Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta09', isotope='177HFm2', gammas=gammas_177Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta10', isotope='177HFm2', gammas=gammas_177Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta11', isotope='177HFm2', gammas=gammas_177Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta12', isotope='177HFm2', gammas=gammas_177Hfm2, unmark_first=True)
# # keep_selected_gammas(foil='Ta13', isotope='177HFm2', gammas=gammas_177Hfm2, unmark_first=True)
# keep_selected_gammas(foil='Ta14', isotope='177HFm2', gammas=gammas_177Hfm2, unmark_first=True)




# lst = get_all_used_gammalines(element='Ta', isotope='179Wg')
ag = AnalyzeGammas(nuclei_nickel)
# print(ag.findGammasSpecificIsotope('177HFm2'))

def go_through(lst, ag, write=False, xrays=False):
    for l in lst:
        a = ag.matchByGamma(gammaLine=l, gammaLineTolerance=1.5, minIntensity=None, xrays=xrays)
        if write:
            for _, row in a.iterrows():
                s = (
                    f"& & \\makecell{{{row['Energy']}}} & "
                    f"\\makecell{{{row['Intensity']} ({row['Unc intensity']})}} &"
                )
                print(s)

# lst = get_all_used_gammalines(element='Ni', isotope='51CRg')
# go_through(lst, ag, write=True, xrays=False) 




# a = AnalyzeGammas(nuclei_nickel).matchByGamma(gammaLine=320.0824, gammaLineTolerance=2.5, minIntensity=None, xrays=False)
