import os
import pandas as pd
import numpy as np


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

def areal_density_from_files(element):
    root = os.getcwd() + '/generatedfiles/arealdensity/'
    for filename in os.listdir(root):
        if element in filename:
                df = pd.read_csv(root + filename)
    areal_density = df['nuclei/cm2'].values
    unc_areal_density = df['unc nuclei/cm2'].values
    foils = df['foil'].values
    return foils, areal_density, unc_areal_density

def beam_current_from_files(stack_files=None):
    pass # TODO must be via file
#     bc = BeamCurrent(stack_files)
#     return bc.weighted_average_beam_current()
