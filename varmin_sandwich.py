import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
import curie as ci
import os
from weighted_average_updated import *

import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

class VarMin:
      
    def __init__(self, dp_array):
        self.dp_array = np.sort(dp_array)

    def get_files(self, dp, stack):
        percentage = f'{dp:.3f}'
        if stack == 'stack_55_MeV':
            full_stack = 'TaSn_stack_55MeV_dp_'+ percentage + '%.csv'
            stack_fluxes = 'TaSn_stack_55MeV_dp_' + percentage  + '%_fluxes.csv'
        elif stack == 'stack_30_MeV':
            full_stack = 'TaSn_stack_30MeV_dp_'+ percentage + '%.csv'
            stack_fluxes = 'TaSn_stack_30MeV_dp_' + percentage  + '%_fluxes.csv'
        else:
            raise Exception('Not a valid stack: ' + stack + '. Must be stack_55_MeV or stack_30_MeV')
        return full_stack, stack_fluxes

    def p0(self, x, a):
        return a
    
    def fit_p0(self, x_data, y_data, unc_data):
        x_data = np.asarray(x_data)
        y_data = np.asarray(y_data)
        unc_data = np.asarray(unc_data)

        mask = np.isfinite(y_data) & np.isfinite(unc_data) & (unc_data > 0)
        x_data = x_data[mask]
        y_data = y_data[mask]
        unc_data = unc_data[mask]

        popt, cov = curve_fit(
            self.p0,
            x_data,
            y_data,
            p0=[100],
            sigma=unc_data,
            absolute_sigma=True
        )
        return popt[0]

    def chi_squared(self, expected, observed, unc_observed):
        observed = np.array(observed)
        expected = np.array(expected)
        chi2 = np.sum((observed - expected)**2 / unc_observed**2)
        return chi2
    
    def get_data_chi_squared(self, x_data, y_data, unc_y_data):
        x = np.array(x_data)
        print('p0')
        # print(len(y_data))
        true = self.fit_p0(x_data, y_data, unc_y_data)
        true_array = np.zeros(len(y_data))
        true_array.fill(true)
        dgf = 1
        return x, y_data, true_array, unc_y_data, dgf

    def run_chi_squared(self, x_data, y_data, unc_y_data):
        x, y_data, true_array, unc_y_data, dgf = self.get_data_chi_squared(x_data, y_data, unc_y_data)
        chi2 = self.chi_squared(y_data, true_array, unc_y_data)
        red_chi2 = chi2/(len(y_data)-dgf)
        # if red_chi2<0.000001:
        return chi2, red_chi2
    
    def get_data(self, number_of_compartments, number_of_monitor_foils, monitor_reactions_per_foil, stack, list_of_isotopes):
        for dp in self.dp_array:
            print("dp: "+ str(dp))
            dp_str = '%.3f' %dp
            full_stack, flux_stack = self.get_files(dp, stack)
            wa = WeightedAverageFlux(flux_stack, full_stack, stack)
            wabc = WeightedAverageBeamCurrent(number_of_compartments, number_of_monitor_foils, monitor_reactions_per_foil, wa, stack, list_of_isotopes)
            bc, unc_bc = wabc.weighted_average_beam_current()
            energy = wabc.flux_weighted_average_energy('Ta')
            print(bc, energy)
            break

    

        indices = self.get_compartment_indices(foils, compartments) if compartments else slice(None)
        energy, flux = wa.get_flux_energy_stack(element)
        flux_weighted_average_energy, unc_energy_left, unc_energy_right = wa.flux_weighted_average_energy(energy, flux)
        flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section = wa.monitor_flux_weighted_average_cross_section(element, isotope)
        flux_weighted_average_energy, unc_energy_left, unc_energy_right, beam_current, unc_beam_current, protons_per_second, unc_protons_per_second = bc.beam_current(element, isotope)
        # print(isotope, beam_current[indices])
        return flux_weighted_average_energy[indices], beam_current[indices], unc_beam_current[indices]

    # def plot(self, compartment):
    #     self.plot_chi2(self.dp_array, red_chi2s,False, title_plot, method,'saddlebrown')
    #     plt.title('Compartment ' + compartments[0])
    #     name = 'p0' + '_comparment' + compartments[0] + '.pdf'
    #     plt.savefig(os.getcwd() + '/generatedfiles/varianceminimization/' + name)
    #     plt.show()


dp_array = np.array([0.990, 0.991, 1.003, 1.001, 1.021, 0.998, 1.013, 
                    1.024, 1.070, 0.996, 1.050, 1.014, 
                    0.900, 1.012, 1.110, 1.200, 1.007, 1.027, 1.040, 0.997, 1.008, 1.016, 1.022, 1.011, 1.025, 0.999, 0.994, 1.005, 1.009, 1.026, 1.080, 1.090, 0.970, 1.023, 1.029, 1.002, 1.019, 1.020, 1.006, 1.010, 0.993, 1.100, 0.950, 1.000, 1.015, 1.030, 0.992, 1.018, 0.800, 1.017, 1.028, 1.060, 1.004, 0.980, 0.960, 0.995])

params_55_stack = [
    'params_Cu_65Zn', 'params_Cu_63Zn', 'params_Cu_62Zn', 'params_Cu_61Cu', 'params_Cu_58Co',
    'params_Ni_57Ni', 'params_Ni_56Ni', 'params_Ni_58Co', 'params_Ni_57Co', 'params_Ni_56Co', 'params_Ni_55Co'
]

params_30_stack = [
    'params_Cu_65Zn', 'params_Cu_63Zn', 'params_Cu_62Zn',
    'params_Ni_56Ni', 'params_Ni_58Co', 'params_Ni_57Co', 'params_Ni_56Co', 'params_Ni_55Co'
]

# 55 MeV

VarMin(dp_array).get_data(7, 2, [6,6], 'stack_55_MeV', params_55_stack)