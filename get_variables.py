import os
import pandas as pd
import numpy as np
from scipy.constants import elementary_charge
from flux_stack import *

def cross_sections(element, isotope, independent):
    root = os.getcwd() + '/generatedfiles/crossections/data_manual/'
    if independent:
        filename = filename = element + '_' + isotope + '_ind.csv'
    elif independent==False:
        filename = filename = element + '_' + isotope + '_cum.csv'
    else:
        filename = filename = element + '_' + isotope + '.csv'
    df = pd.read_csv(root + filename, comment='#')
    cross_sections = df['cross_section'].values
    unc_cross_sections = df['unc_cross_section'].values
    energy = df['energy'].values
    unc_left = df['unc_left'].values
    unc_right = df['unc_right'].values
    return energy, unc_left, unc_right, cross_sections, unc_cross_sections

def eob_activity_manually(element, isotope, independent, stack=None):
    root = os.getcwd() + '/generatedfiles/activity/data_isotope2/'
    if independent==True:
        filename = element + '_' + isotope + '_ind.csv'
    elif independent==False:
        filename = element + '_' + isotope + '_cum.csv'
    else:
        filename = element + '_' + isotope + '.csv'
    print(filename)
    df = pd.read_csv(root + filename)
    # print(df)
    eob_activity = df['fit'].values
    cov_eob_activity = df['cov'].values
    unc_eob_activity = np.sqrt(cov_eob_activity)
    std_eob_activity = np.where(
            eob_activity > 0,
            (unc_eob_activity / eob_activity)**2,
            0.0)
    print(stack)
    if stack:
        start_idx, end_idx = get_indexes_stack(stack)
        return eob_activity[start_idx:end_idx], std_eob_activity[start_idx:end_idx]
    else:
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

def weighted_average_beam_energy(element):
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

    df = pd.concat([df_55, df_30])
    energy = df['wabe'].values
    unc_left = df['unc wabe left'].values
    unc_right = df['unc wabe right'].values
    return energy, unc_left, unc_right

def weighted_average_beam_current():
    # run when plotting....
    root = os.getcwd() + '/generatedfiles/beamcurrent/'
    stack_30 = 'beam_current_stack_30_MeV.csv'
    stack_55 = 'beam_current_stack_55_MeV.csv'
    df_30 = pd.read_csv(root + stack_30)
    df_55 = pd.read_csv(root + stack_55)
    df = pd.concat([df_55, df_30])
    beam_current = df['beam current (nA)'].values
    unc_beam_current = df['unc beam current (nA)'].values
    print(beam_current)
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

# i, di = weighted_average_beam_current()
# print(i)