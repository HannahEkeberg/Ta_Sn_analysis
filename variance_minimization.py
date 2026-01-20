from flux_stack import WeightedAverageFlux
from beamcurrent import BeamCurrent
import os
import numpy as np
import matplotlib.pyplot as plt

# wa = WeightedAverageFlux('TaSn_stack_55MeV_fluxes.csv', 'TaSn_stack_30MeV_fluxes.csv', 'TaSn_stack_55MeV.csv','TaSn_stack_30MeV.csv')


dp_array = np.arange(0.90, 1.10, 0.01)

"""
for dp in dp_array[:3]:
    percentage = f'{dp:.3f}'
    stack_55_fluxes = 'TaSn_stack_55MeV_dp_' + percentage  + '%_fluxes.csv'
    stack_30_fluxes = 'TaSn_stack_30MeV_dp_' + percentage  + '%_fluxes.csv'
    stack_55 = 'TaSn_stack_55MeV_dp_'+ percentage + '%.csv'
    stack_30 = 'TaSn_stack_30MeV_dp_'+ percentage + '%.csv'
    wa = WeightedAverageFlux(stack_55_fluxes, stack_30_fluxes, stack_55, stack_30)
    print(wa)
    bc = BeamCurrent(stack_files=[stack_55_fluxes, stack_30_fluxes, stack_55, stack_30])
    I = bc.beam_current(element='Cu', isotope='63ZN')
    # print(I)
    # bc.cross_sections(element='Cu', isotope='63ZN')
    # print(wa)



    # print(I)

"""

class VarianceMinimization:

    def __init__(self, dp_array):
        self.dp_array = dp_array
        pass

    def get_files(self, dp, stack):
        percentage = f'{dp:.3f}'
        if stack == '55_MeV':
            full_stack = 'TaSn_stack_55MeV_dp_'+ percentage + '%.csv'
            stack_fluxes = 'TaSn_stack_55MeV_dp_' + percentage  + '%_fluxes.csv'
        elif stack == '30_MeV':
            full_stack = 'TaSn_stack_30MeV_dp_'+ percentage + '%.csv'
            stack_fluxes = 'TaSn_stack_30MeV_dp_' + percentage  + '%_fluxes.csv'
        else:
            raise Exception('Not a valid stack: ' + stack + '. Must be 55_MeV or 30_MeV')
        return full_stack, stack_fluxes

    def p0(self, x, a):
        a = np.ones(len(x))
        return a

    def beam_current(self, element, isotope):
        flux_weighted_average_energy, unc_energy_left, unc_energy_right, beam_current, unc_beam_current, protons_per_second, unc_protons_per_second = bc.beam_current(element, isotope)

    def chi_squared(self, expected, observed, unc_observed):
        observed = np.array(observed)
        expected = np.array(expected)
        diff = observed-expected 
        chi2 = np.sum(np.multiply(diff, diff)/np.mcultiply(unc_observed,unc_observed))
        return chi2
    
    def plot_chi_squared(self, element, isotope, stack):
        for dp in self.dp_array:
            # flux, stack = self.get_files(dp, stack)
            full_stack, flux_stack = self.get_files(dp, stack)
            wa = WeightedAverageFlux(stack_55_fluxes, stack_30_fluxes, stack_55, stack_30)
            print(wa)
            # bc = BeamCurrent(stack_files=[stack_55_fluxes, stack_30_fluxes, stack_55, stack_30])
            # I = bc.beam_current(element='Cu', isotope='63ZN')

            
            
            
        pass

VarianceMinimization(dp_array[:3]).plot_chi_squared('Cu', '63ZN', '55_MeV')



"""
bc.plot_beam_current(element='Cu', isotope='63ZN', stack=None, color='blue', label='63Zn')
bc.plot_beam_current(element='Cu', isotope='65ZN', stack=None, color='red', label='65Zn')
bc.plot_beam_current(element='Cu', isotope='62ZN', stack=None, color='yellow', label='63Zn')
bc.plot_beam_current(element='Cu', isotope='56CO', stack=None, color='orange', label='56Co')
bc.plot_beam_current(element='Cu', isotope='58CO', stack=None, color='darkred', label='58Co')
bc.plot_beam_current(element='Ni', isotope='57NI', stack=None, color='cyan', label='Ni 57Ni')
"""
