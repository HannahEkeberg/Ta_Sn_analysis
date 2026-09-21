import pandas as pd
import os 
import numpy as np

path = os.getcwd() + '/generatedfiles/peakdata/data_isotope/'

def concat_peak_data(foil):
    files = [f for f in os.listdir(path) if foil in f]
        # files = [f for f in files if '179W' in f]
    dfs = []
    for f in files:
        df = pd.read_csv(path + f)
        dfs.append(df)
    dfs = [df for df in dfs if not df.empty]
    total_df = pd.concat(dfs, ignore_index=True).sort_values(['filename', 'energy'])
    return total_df

def get_peaks(foil, energy, tolerance):
    total_df = concat_peak_data(foil)
    E_lim = ( energy - tolerance, energy + tolerance )
    matching_rows = total_df[total_df['energy'].between(E_lim[0], E_lim[1])]
    print(matching_rows)


# concat_peak_data('Ta01')
# get_peaks('Ta01', 103.600, 2.5)
get_peaks('Ta01', 93.324, 2.5)
# def correct_counts(energy):
#     A0_interference=1; efficiency=1; intensity=1; 
#     decay_const_interference=1; delay_time=1; irradiation_time=1
#     counts_E = (A0_interference * efficiency * intensity) / (decay_const_interference) * np.exp(-decay_const_interference*delay_time) * (1 - np.exp(-decay_const_interference*irradiation_time))
#     pass

# def correct_counts2(energy, df_corrected, df_interference):
#     energy_window_keV = 1
#     energy_matches_df = df_corrected[df_corrected['energy'] == energy]
#     print(energy_matches_df)
#     interference = df_interference[
#         df_interference['energy'].between(
#         energy - energy_window_keV,
#         energy + energy_window_keV
#         )]
    
#     print(interference)
#     matching_interference = interference[
#         interference['filename'].isin(energy_matches_df['filename'])
#         ]
#     corrected_rows = []
#     for i, interference_row in matching_interference.iterrows():
#         mask = (
#             (df_corrected['energy'] == energy) &
#             (df_corrected['filename'] == interference_row['filename'])
#             )
    
#         target_row = energy_matches_df[energy_matches_df['filename'] == interference_row['filename']].iloc[0]
#         # counts_corrected = target_row['counts'] - interference_row['counts'] * (
#         #     (target_row['efficiency'] * target_row['intensity']) /
#         #     (interference_row['efficiency'] * interference_row['intensity'])
#         # )
        
#         scale = (
#             (target_row['efficiency'] * target_row['intensity']) /
#             (interference_row['efficiency'] * interference_row['intensity'])
#         )
#         counts_corrected = target_row['counts'] - interference_row['counts'] * scale
#         sigma_corrected = np.sqrt(
#         target_row['counts'] + (scale**2) * interference_row['counts']
#         )
#         print(target_row['isotope'], target_row['filename'], target_row['counts'])
#         print(interference_row['isotope'], interference_row['filename'], interference_row['counts'])
#         print('counts corrected: ' + str(counts_corrected) + ' +/- ' + str (sigma_corrected))
#         if counts_corrected > 0:
#             print("****")
#             print(target_row['counts'], counts_corrected)
#             print(df_corrected[mask])
#             df_corrected.loc[mask, 'counts'] = counts_corrected
#             print(df_corrected[mask])
    


# def update(foil, isotope, isotope_interference, gamma):
#     # df_combined = pd.concat(
#     #     [pd.read_csv(os.path.join(path, file)) for file in os.listdir(path)if foil in file]
#     #     ,ignore_index=True).sort_values("energy")
#     # df_target = foil + '_' + isotope + '_gammas.csv'
#     # df_corrected = pd.read_csv(path + df_target).copy()

#     # energy_tolerance = 1
#     # interference = df_combined[
#     #     df_combined['energy'].between(
#     #     gamma - energy_tolerance,
#     #     gamma + energy_tolerance
#     #     )]

#     # print(interference)

#     peak_data_file = foil + '_' + isotope + '_gammas.csv'
#     peak_data_file_interference = foil + '_' + isotope_interference + '_gammas.csv'
#     df = pd.read_csv(path + peak_data_file)
#     df_corrected = df.copy()
#     df_interference = pd.read_csv(path + peak_data_file_interference)
#     correct_counts2(gamma, df_corrected, df_interference)

# # update('Ta03', '180TAg', '178TAg', 93.324)
# update('Ta04', '180TAg', '175TAg', 103.6)


# # foil = 'Ta01'; isotope='180TAg'; isotope_interference = '178TAg'
# # peak_data_file = foil + '_' + isotope + '_gammas.csv'
# # peak_data_file_interference = foil + '_' + isotope_interference + '_gammas.csv'

# # df = pd.read_csv(path + peak_data_file)
# # df_corrected = df.copy()
# # df_interference = pd.read_csv(path + peak_data_file_interference) # her kan man sikkert concate fra andre duplikater... 
# # unique_energies = df['energy'].unique()

# # gamma = 103.6; 
# # isotope_interference = '178TAg'

# # correct_counts2(93.324)
