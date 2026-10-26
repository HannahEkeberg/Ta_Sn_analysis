import pandas as pd
import numpy as np
import curie as ci 
import os 
from activity_new import *

activity = Activity_manual()

def search(foil, nuclide, energy):
    print('Searching for ' + nuclide + ' in ' + foil)
    folder = os.getcwd() + '/generatedfiles/peakdata/data/'
    # foil = 'Ta02'
    # nuclide = '180TAg'
    # energy = 103.6
    # energy = 93.324
    energy_tol = 1
    files = [f for f in os.listdir(folder) if foil in f and f.endswith(".csv")]
    df_full = pd.DataFrame()
    # print(df_full)
    for f in files:
        # print(f)
        df = pd.read_csv(folder + f)
        # print(df.columns)
        # df_sorted = df[df['isotope'] == nuclide]
        matches = df[df["energy"].between(energy - energy_tol, energy + energy_tol)]
        if not matches.empty:
            filename  = matches['filename'].str.split("/").str[-1]
            # print(filename)
            # print(matches[['isotope', 'energy']])
            df_full = pd.concat([df_full, matches], ignore_index=True)
            # print("-----")
    print(df_full)
    return df_full

def search_in_isotopes(foil, energy):
    folder = os.getcwd() + '/generatedfiles/peakdata/data_isotope/'
    energy_tol = 1
    files = [f for f in os.listdir(folder) if foil in f and f.endswith(".csv")]
    print(files)
    df_full = pd.DataFrame()
    for f in files:
        df = pd.read_csv(folder + f)
        matches = df[df["energy"].between(energy - energy_tol, energy + energy_tol)]
        if not matches.empty:
            filename  = matches['filename'].str.split("/").str[-1]
            df_full = pd.concat([df_full, matches], ignore_index=True)
    print(df_full[['filename', 'isotope', 'energy', 'intensity']])
    # print(df_full.columns)
    return df_full


def get_vals_per_isotope(df, isotope, eob):
    peak_data = df[df['isotope'] == isotope]
    E = peak_data['energy']
    filename = peak_data['filename'].apply(os.path.basename)
    counts = peak_data['counts']
    unc_counts = peak_data['unc_counts']
    intensity = peak_data['intensity']
    unc_intensity = peak_data['unc_intensity']
    efficiency = peak_data['efficiency']
    unc_efficiency = peak_data['unc_efficiency']
    count_time = peak_data['live_time']
    iso = peak_data['isotope']
    eob = pd.to_datetime(eob)
    peak_data['start_time'] = pd.to_datetime(peak_data['start_time'])
    delay_time = (peak_data['start_time'] - eob).dt.total_seconds() 

    return filename, E, counts, intensity, efficiency, count_time, delay_time

def plot_decay_rate(foil, peak_data, isotope1, isotope2):
    eob = activity.getEob(foil)
    # print(df['isotope'])
    df = peak_data.sort_values(by='start_time', ascending=True).reset_index(drop=True)
    
    peak_data = df[df['isotope'].isin([isotope1, isotope2])]
    # print(peak_data)
    E = peak_data['energy']
    filename = peak_data['filename'].apply(os.path.basename)
    counts = peak_data['counts']
    unc_counts = peak_data['unc_counts']
    intensity = peak_data['intensity']
    unc_intensity = peak_data['unc_intensity']
    efficiency = peak_data['efficiency']
    unc_efficiency = peak_data['unc_efficiency']
    count_time = peak_data['live_time']
    iso = peak_data['isotope']
    eob = pd.to_datetime(eob)
    peak_data['start_time'] = pd.to_datetime(peak_data['start_time'])
    delay_time = (peak_data['start_time'] - eob).dt.total_seconds()
    # activity = 
    
    df['delay_time'] = delay_time/3600
    print(df[['filename', 'isotope', 'energy', 'intensity', 'counts', 'delay_time']])
    
    dc_isotope1 = ci.Isotope(isotope1).decay_const()
    half_life = ci.Isotope(isotope1).half_life()
    dc_isotope2 = ci.Isotope(isotope2).decay_const()
    # dc_isotope2 = ci.Isotope('178Wg').decay_const()


    def decay(t, A1, A2):
        return (
            A1 * np.exp(-dc_isotope1 * t) +
            A2 * np.exp(-dc_isotope2 * t)
        )

    popt, pcov = curve_fit( decay, delay_time, peak_data['counts'], sigma=peak_data['unc_counts'], absolute_sigma=True, p0=[1000, 1000], bounds=(0, np.inf))

    A1, A2 = popt

    # print(A1, A2)
    t = np.linspace(delay_time.min(), delay_time.max(), 500) 

    component1 = A1 * np.exp(-dc_isotope1 * t)
    component2 = A2 * np.exp(-dc_isotope2 * t)

    eff1 = peak_data[peak_data['isotope'] == isotope1]['efficiency']
    eff2 = peak_data[peak_data['isotope'] == isotope2]['efficiency']
    
    energy1 = peak_data[peak_data['isotope'] == isotope1]['energy']
    energy2 = peak_data[peak_data['isotope'] == isotope2]['energy']
    
    I_gamma1 = peak_data[peak_data['isotope'] == isotope1]['intensity']
    I_gamma2 = peak_data[peak_data['isotope'] == isotope2]['intensity']

    live_time1 = peak_data[peak_data['isotope'] == isotope1]['live_time']
    live_time2 = peak_data[peak_data['isotope'] == isotope2]['live_time']

    A1_Bq = A1 / (efficiency * I_gamma1 * live_time1)
    A2_Bq = A2 / (efficiency * I_gamma2 * live_time2)

    # print(live_time1, filename, A1_Bq)
    # print(A2_Bq)

    for a, b, c in zip(live_time1, filename, A1_Bq):
        print(f"{a/3600}.2f,{b},{c}")

    plt.xlim(0,10*half_life/3600)
    plt.plot(delay_time/3600, A1_Bq*1e-6, 'o', label=isotope1)
    plt.plot(delay_time/3600, A2_Bq*1e-6, 'o', label=isotope2)
    plt.legend()
    plt.xlabel('Time since eob hours')
    plt.ylabel('Activity MBq')
    plt.show()

def try_to_calculate_subtracted_A0(isotope_feeding, foil, independent, df_data):

    end = 'ind' if independent else 'cum'
    filename = foil + '_' + isotope_feeding + '_' + end + '.csv'
    activity_file = os.getcwd() + '/generatedfiles/activity/data_isotope2/' + filename
    df_activity = pd.read_csv(activity_file) # A0, dA0
    print(df_activity)

    activity_file_180Ta = os.getcwd() + '/generatedfiles/activity/data_isotope2/' + foil + '_' + '180TA' + '_ind' + '.csv'

    df  = pd.read_csv(activity_file_180Ta)
    df_new = df.copy()

    A0_180Ta_subtract = df['fit'] - df_activity['fit']

    df_new['fit'] = A0_180Ta_subtract

    print(df_new)

    df_new.to_csv(os.getcwd () + '/generatedfiles/activity/data_isotope2/' + foil + '_' + '180TA_corrected' + '_ind' + '.csv')


    # # FOR FEEDING ISOTOPE
    # dc = ci.Isotope(isotope_feeding).decay_const()
    # eob = activity.getEob(foil + '01')
    # end = 'ind' if independent else 'cum'
    # filename = foil + '_' + isotope_feeding + '_' + end + '.csv'
    # df_activity = pd.read_csv(activity_file) # A0, dA0
    # activity_file = os.getcwd() + '/generatedfiles/activity/data_isotope2/' + filename
    

    # eob = pd.to_datetime(eob)
    # df_data['start_time'] = pd.to_datetime(df_data['start_time'])
    # delay_time = (df_data['start_time'] - eob).dt.total_seconds()
    # df_data['delay_time'] = delay_time
    # # df_filter = df_data[df_data['isotope'].isin(['180TAg', '178TAm1'])]
    # df_filter = df_data[df_data['isotope'].isin(['178TAm1'])]
    # df_sorted = df_filter[['filename', 'isotope', 'energy', 'intensity', 'counts', 'delay_time']].sort_values(['delay_time', 'energy'])
    # # print(df_sorted)


    # N_A = df_activity['fit'][0] * np.exp(-dc*df_sorted['delay_time']) * df_sorted['intensity'] * df_sorted['efficiency'] * 
    # print(len(N_A))
    # print(N_A)

    # E_gamma = 93.324


    

# Incident where production of 181W, 179W and 180Ta seem to increase for the highest energies
#


# df = search('Ta01', '180TAg', 103.6)
# df = search('Ta09', '180TAg', 93.324)
# plot_decay_rate('Ta09', df, '180TAg', '180HFm')
# df = search_in_isotopes('Ta01', 332.274)
# plot_decay_rate('Ta01', df, '180HFm', '178TAm1')



df = search_in_isotopes('Ta01', 93.324)
# plot_decay_rate('Ta01', df, '180TAg', '178TAm1')
# df = search_in_isotopes('Ta03', 103.6)
# plot_decay_rate('Ta03', df, '180TAg', '175TAg')

try_to_calculate_subtracted_A0('178TAm1', 'Ta', False, df)



# df = search('Ta01', '180TAg', 103.6)
# plot_decay_rate('Ta01', df, '180TAg', '180HFm')
# df = search('Ta01', '180TAg', 103.6)
# plot_decay_rate('Ta04', df, '180TAg', '180HFm')

# df = search('Ta05', '180TAg', 93.324)
# plot_decay_rate('Ta05', df, '180TAg', '180HFm')