import os
import pandas as pd
import numpy as np
from scipy.constants import elementary_charge
from flux_stack import *

def cross_sections(element, isotope, independent):
    root = os.getcwd() + '/generatedfiles/crossections/data_manual/'

    # print(type(element))
    if independent:
        filename = element + '_' + isotope + '_ind.csv'
    elif independent==False:
        filename = element + '_' + isotope + '_cum.csv'
    else:
        filename = element + '_' + isotope + '.csv'
    # print(root + filename)
    print(filename)
    df = pd.read_csv(root + filename, comment='#')
    print(df)
    cross_sections = df['cross section (mb)'].values
    unc_cross_sections = df['unc cross section (mb)'].values
    energy = df['energy (MeV)'].values
    unc_left = df['unc left (MeV)'].values
    unc_right = df['unc right (MeV)'].values
    return energy, unc_left, unc_right, cross_sections, unc_cross_sections

def eob_activity_manually(element, isotope, independent, stack=None):
    root = os.getcwd() + '/generatedfiles/activity/data_isotope2/'
    if independent==True:
        filename = element + '_' + isotope + '_ind.csv'
    elif independent==False:
        filename = element + '_' + isotope + '_cum.csv'
    else:
        filename = element + '_' + isotope + '.csv'
    # print(filename)
    df = pd.read_csv(root + filename)
    # print(df)
    eob_activity = df['fit'].values
    cov_eob_activity = df['cov'].values
    unc_eob_activity = np.sqrt(cov_eob_activity)
    std_eob_activity = np.where(
            eob_activity > 0,
            (unc_eob_activity / eob_activity)**2,
            0.0)    
    # print(eob_activity)
    # print("eob activities unx")
    # print(cov_eob_activity)
    # print(unc_eob_activity)
    # print(std_eob_activity)
    # print("eob activities unx")
    if stack:
        start_idx, end_idx = get_indexes_stack(stack)
        # return eob_activity[start_idx:end_idx], std_eob_activity[start_idx:end_idx]
        return eob_activity[start_idx:end_idx], cov_eob_activity[start_idx:end_idx]
    else:
        # return eob_activity, cov_eob_activity #std_eob_activity
        return eob_activity, std_eob_activity

def eob_activity_curie(element, isotope, stack=None, independent=None):
    root = os.getcwd() + '/generatedfiles/activity/data_isotope/'
    if independent==True:
        filename = element + '_' + isotope + '_ind.csv'
    elif independent==False:
        filename = element + '_' + isotope + '_cum.csv'
    else:
        filename = element + '_' + isotope + '.csv'
    df = pd.read_csv(root + filename)
    eob_activity = df['fit'].values
    cov_eob_activity = df['cov'].values
    unc_eob_activity = np.sqrt(cov_eob_activity)
    std_eob_activity = np.where(
            eob_activity > 0,
            (unc_eob_activity / eob_activity)**2,
            0.0)
    if stack:
        start_idx, end_idx = get_indexes_stack(stack)
        return eob_activity[start_idx:end_idx], std_eob_activity[start_idx:end_idx]
    else:
        return eob_activity, std_eob_activity

def eob_activity_from_files(foils, isotope):
    root = os.getcwd() + '/generatedfiles/activity/data/'
    eob_activity = np.zeros(len(foils)); cov_eob_activity = np.zeros(len(foils))
    for i, foil in enumerate(foils):
        for filename in os.listdir(root):
            if foil in filename and 'all_isotopes' in filename:
                df = pd.read_csv(root + filename)
                df_isotope = df[df['isotope'].astype(str).str.contains(isotope, case=False, na=False)]
                if df_isotope.empty:
                    eob_activity[i] = 0.0; cov_eob_activity[i]=0.0
                else:
                    # TODO add if test for array length > 1 
                    if len(df_isotope['fit'].values) > 1:
                        pass
                    else:
                        eob_activity[i] = df_isotope['fit'].values[0]
                        if eob_activity[i] == 0:
                            cov_eob_activity[i]  = 0
                        else:
                            cov_eob_activity[i] = df_isotope['cov'].values[0]
    return eob_activity, cov_eob_activity

def areal_density_from_files(element, stack=None):
    root = os.getcwd() + '/generatedfiles/arealdensity/'
    for filename in os.listdir(root):
        if element in filename:
                df = pd.read_csv(root + filename)
    areal_density = df['nuclei/cm2'].values
    unc_areal_density = df['unc nuclei/cm2'].values
    foils = df['foil'].values
    if stack:
        start_idx,end_idx = get_indexes_stack(stack)
        return foils[start_idx:end_idx], areal_density[start_idx:end_idx], unc_areal_density[start_idx:end_idx]
    return foils, areal_density, unc_areal_density

def weighted_average_beam_energy(element, stack=None):
    # run when flux stack. 
    root = os.getcwd() + '/generatedfiles/fluxweightedaverageenergy/'
    if element == 'Ni':
        stack_30 = 'stack_30_MeV_Ni_weighted_average_beam_energy.csv'
        stack_55 = 'stack_55_MeV_Ni_weighted_average_beam_energy.csv'
    elif element == 'Cu':
        stack_30 = 'stack_30_MeV_Cu_weighted_average_beam_energy.csv'
        stack_55 = 'stack_55_MeV_Cu_weighted_average_beam_energy.csv'
    elif element == 'Ta':
        stack_30 = 'stack_30_MeV_Ta_weighted_average_beam_energy.csv'
        stack_55 = 'stack_55_MeV_Ta_weighted_average_beam_energy.csv'
    elif element == 'Sn':
        stack_30 = 'stack_30_MeV_Sn_weighted_average_beam_energy.csv'
        stack_55 = 'stack_55_MeV_Sn_weighted_average_beam_energy.csv'
    else: 
        raise Exception('Invalid element: ' + element)
    
    df_30 = pd.read_csv(root + stack_30)
    df_55 = pd.read_csv(root + stack_55)
    if stack == 'stack_30_MeV':
        df = df_30
    elif stack == 'stack_55_MeV':
        df = df_55
    else:
        df = pd.concat([df_55, df_30])
    energy = df['wabe'].values
    unc_left = df['unc wabe left'].values
    unc_right = df['unc wabe right'].values
    return energy, unc_left, unc_right

def weighted_average_beam_current(file='IAEA25'):
    # run when plotting....
    print("*****")
    root = os.getcwd() + '/generatedfiles/beamcurrent/true_weighted_average_beam_current/'
    # IAEA 17
    # IAEA 25 - old monitors
    # IAEA 25 - all values.... 
    # IAEA 25 - Selected parameters 
    print(root)
    if file=='IAEA17':
        # generatedfiles/beamcurrent/true_weighted_average_beam_current/old_old_params_30.csv
        print("OLD IAEA DATA...")
        # stack_30 = 'old_old_params_30.csv'    # old iaea data
        # stack_55 = 'old_old_params_55.csv'    # old iaea data
        stack_30 = 'iaea17_30MeV_dp1.010.csv'    # old iaea data
        stack_55 = 'iaea17_55MeV_dp1.020.csv'    # old iaea data


        # iaea17_30MeV_dp1.006.csv
    elif file=='IAEA17_not_weighted':
        # generatedfiles/beamcurrent/true_weighted_average_beam_current/old_old_params_30.csv
        print("OLD IAEA DATA")
        stack_30 = 'beam_current_stack_30_MeV_not_weighted.csv'    # old iaea data
        stack_55 = 'beam_current_stack_55_MeV_not_weighted.csv'    # old iaea data
    elif file == 'IAEA25':
        stack_30 = 'iaea25_30MeV_dp1.006.csv'  # new iaea data
        stack_55 = 'iaea25_55MeV_dp1.016.csv'  # new iaea dat
    else:
        stack_30 = 'iaea25_30MeV_dp1.006.csv'  # new iaea data
        stack_55 = 'iaea25_55MeV_dp1.016.csv'  # new iaea dat
    df_30 = pd.read_csv(root + stack_30)
    df_55 = pd.read_csv(root + stack_55)
    df = pd.concat([df_55, df_30])
    # print(df)
    # print(file)
    # print(df)
    beam_current = df['Average Current (nA)'].values.copy()
    unc_beam_current = df['Uncertainty in Average Current (nA)'].values.copy()
    # print(beam_current)
    #convert to protons/s:
    beam_current *= 1e-9/elementary_charge
    unc_beam_current *= 1e-9/elementary_charge
    return beam_current, unc_beam_current

def get_indexes_stack(stack):
    if stack == 'stack_55_MeV':
        start_idx = 0; end_idx = 7
    elif stack == 'stack_30_MeV':
        start_idx = 7; end_idx = 14
    else:
        print(stack)
        raise Exception('Not a valid stack: ' + stack + '. Must be stack_30_MeV or stack_55_MeV')
    return start_idx, end_idx

def get_wa_from_stack_files():
    stack_fluxes_30 = 'TaSn_stack_30MeV_dp_1.010%_fluxes.csv'; full_stack_30 = 'TaSn_stack_30MeV_dp_1.010%.csv'
    stack_fluxes_55 = 'TaSn_stack_55MeV_dp_1.020%_fluxes.csv'; full_stack_55 = 'TaSn_stack_55MeV_dp_1.010%.csv'
    wa_55 = WeightedAverageFlux(stack_fluxes = stack_fluxes_55, full_stack = full_stack_55, stack = 'stack_55_MeV')
    wa_30 = WeightedAverageFlux(stack_fluxes = stack_fluxes_30, full_stack = full_stack_30, stack = 'stack_30_MeV')
    return wa_55, wa_30

def update_weighted_average_beam_current():
    pass
def update_weighted_average_beam_energy():
    pass

    
def get_irradiation_time(stack=None):
    # TODO fix in bc also...
    print(stack)
    stack_30_start = pd.to_datetime('09/23/2025 17:36:59')
    stack_30_end   = pd.to_datetime('09/23/2025 18:40:05')
    stack_55_start = pd.to_datetime('09/24/2025 14:45:05')
    stack_55_end   = pd.to_datetime('09/24/2025 15:45:32')
    unc = 3
    irr_time_55 = (stack_55_end - stack_55_start).total_seconds()
    irr_time_30 = (stack_30_end - stack_30_start).total_seconds()
    if stack == 'stack_55_MeV':
        irr_time = np.full(7, irr_time_55)
        unc_irr_time = np.full(7,unc)
    elif stack == 'stack_30_MeV':
        irr_time = np.full(7, irr_time_30)
        unc_irr_time = np.full(7,unc)
    else:
        irr_time = np.concatenate([np.full(7, irr_time_55), np.full(7, irr_time_30)])  # → 14 elements
        unc_irr_time = np.full(14,unc)
    return irr_time, unc_irr_time


def foils(element):
    stack_numbs = ['01', '02', '03', '04', '05', '06', '07', '08', '09','10', '11', '12', '13', '14']
    foils = [element + number for number in stack_numbs]
    return foils