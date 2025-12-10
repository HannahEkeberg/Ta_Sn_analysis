import curie as ci
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.constants import elementary_charge
import sys
sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
from nuclearanalysistools.Tendl import *
from nuclearanalysistools.tools import *

class WeightedAverageFlux:

    def __init__(self, stack_55_fluxes, stack_30_fluxes, stack_55, stack_30):
        self.root = os.getcwd() + '/generatedfiles/stack/'
        self.stack_55_fluxes = pd.read_csv(self.root + stack_55_fluxes)
        self.stack_30_fluxes = pd.read_csv(self.root + stack_30_fluxes)
        self.stack_55       = pd.read_csv(self.root + stack_55)
        self.stack_30       = pd.read_csv(self.root + stack_30)

    def stack(self, element):
        stack_55 = self.stack_55[self.stack_55['compound'] == element]
        stack_30 = self.stack_30[self.stack_30['compound'] == element]
        E_55 = stack_55['mu_E'].values; E_30 = stack_30['mu_E'].values
        return np.concatenate((E_55, E_30))
    
    def get_flux_energy_stack(self, element):
        flux = []; energy = []
        foils = self.foils(element)
        total_stack = pd.concat([self.stack_55_fluxes, self.stack_30_fluxes], ignore_index=True)
        for foil in foils:
            filtered_stack = total_stack[total_stack['name'] == foil]
            E = filtered_stack['energy'].values; F = filtered_stack['flux'].values
            energy.append(E); flux.append(F)
        return energy, flux

    def plot_distributions(self, element):
        foils = self.foils(element)
        energy, flux = self.get_flux_energy_stack(element)
        colors = Tools().colors()
        mean_energy, unc_energy_left, unc_energy_right, fwhm, half_max = self.variables_for_flux_weighted_average(energy, flux)
        
        for i in range(len(energy)):
            plt.plot(energy[i], flux[i], color=colors[i], linewidth = 0.7)
            half_flux = np.max(flux[i])/2
            plt.plot(half_max[i], [half_flux, half_flux], color=colors[i], linewidth=0.8, label=foils[i]+ ' - {0:.2f}'.format(mean_energy[i])  + ' MeV (FWHM: {0:.2f})'.format(fwhm[i]))
            plt.vlines(mean_energy[i], ymin=0.0, ymax = np.max(flux[i]), linewidth=0.4, linestyle='--')#, label=r'$\mu=${}'.format(mu))
            plt.legend(fontsize='xx-small')
            plt.show()

    def foils(self,element):
        stack_numbs = ['01', '02', '03', '04', '05', '06', '07', '08', '09','10', '11', '12', '13', '14']
        foils = [element + number for number in stack_numbs]
        return foils

    def flux_weighted_average_energy(self, energy, flux, plot_distribution=False, foils=None):
        mean_energy, unc_energy_left, unc_energy_right, fwhm, half_max = self.variables_for_flux_weighted_average(energy, flux)
        return mean_energy, unc_energy_left, unc_energy_right
    
    def variables_for_flux_weighted_average(self, energy, flux):
        unc_energy_left = np.zeros(len(energy)); unc_energy_right = np.zeros(len(energy))
        fwhm = np.zeros(len(energy)); half_max = []; mean_energy = np.zeros(len(energy))
        for i in range(len(energy)):
            max_flux = np.max(flux[i])
           
            def line_interpolation(x, y, i, half):
                return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

            def half_max_flux_energy(E,F):
                half_max_flux = max(F)/2.0
                signs = np.sign(np.add(F, -half_max_flux))  # for each flux, if over half max +1, if under -1, else 0
                zero_crossings = (signs[0:-2] != signs[1:-1]) # find all the points where the flux is over next to something under.
                # gives multiple arrays of len 2, just take the first in case there are more than 1 crossing. 
                # Arrays with coordinates on the left and right side?
                zero_crossings_i = np.where(zero_crossings)[0] 
                return [line_interpolation(E, F, zero_crossings_i[0], half_max_flux),
                        line_interpolation(E, F, zero_crossings_i[1], half_max_flux)]
            
            half_max_flux_energy = half_max_flux_energy(energy[i], flux[i])
            half_max.append(half_max_flux_energy)
            mean_energy[i] = np.trapezoid(flux[i]*energy[i], energy[i])/np.trapezoid(flux[i],energy[i])
            fwhm[i] = half_max_flux_energy[1]-half_max_flux_energy[0]
            unc_energy_left[i] = mean_energy[i]-half_max_flux_energy[0]; unc_energy_right[i] = half_max_flux_energy[1]-mean_energy[i]   #left and right uncertainty in energy

        return mean_energy, unc_energy_left, unc_energy_right, fwhm, half_max
        
    
    def monitor_data(self, element, isotope):
        if element == 'Cu' and isotope =='62ZN':
            monitor_file = './monitordata/cup62znt/cup62znt.txt'
            useFile=True
        elif element == 'Cu' and isotope =='63ZN':
            monitor_file = './monitordata/cup63znt/cup63znt.txt'
            useFile=True
        elif element == 'Cu' and isotope =='65ZN':
            monitor_file = 'monitordata/cup65znt/cup65znt.txt'
            useFile=True
        elif element == 'Cu' and isotope =='58CO':
            monitor_file = 'monitordata/cup58cot/cup58cot.txt'
            useFile=True
        elif element == 'Cu' and isotope =='56CO':
            monitor_file = 'monitordata/cup56cot/cup56cot.txt'
            useFile=True
        elif element == 'Ni' and isotope =='57NI':
            monitor_file = 'monitordata/nip57nit/cup57nit.txt'
            useFile=True
        elif element == 'Ni' and 'isotope' == '55CO':
            useFile=False
        
        energy = np.loadtxt(monitor_file, usecols=[0], skiprows=6)
        cs = np.loadtxt(monitor_file, usecols=[1], skiprows=6)
        unc_cs = np.loadtxt(monitor_file, usecols=[2], skiprows=6)
        if cs[0]>0:
            energy_below_interpolation = np.linspace(0,energy[0]-0.5, 4)
            cs = np.concatenate((np.zeros(len(energy_below_interpolation)), cs))
            energy = np.concatenate((energy_below_interpolation, energy))
            unc_cs = np.concatenate((np.zeros(len(energy_below_interpolation)), unc_cs))

        tck = interpolate.splrep(energy, cs, s=0)
        sigma_tck = interpolate.splrep(energy, unc_cs, s=0)

        # TODO for remaining, fit experimental data to pade parameters from https://doi.org/10.1007/s10967-024-09513-7
        return energy, cs, unc_cs, tck, sigma_tck
    
    def plot_monitor_reaction(self, element, isotope, label='IAEA recommended data'):
        energy, cs, unc_cs, tck, sigma_tck = self.monitor_data(element, isotope)
        plt.plot(energy, cs, label=label)
        plt.fill_between(energy, cs+unc_cs, cs-unc_cs, color='blue', alpha=0.1)
        plt.legend()
        plt.show()


    def flux_weighted_average_cross_section(self, element, foils, isotope):
        mon_energy, mon_cs, mon_unc_cs, tck, sigma_tck = self.monitor_data(element, isotope)
        energy, flux = self.get_flux_energy_stack(foils)
        flux_weighted_average_cross_section = np.zeros(len(energy))
        unc_flux_weighted_average_cross_section = np.zeros(len(energy))
        interpolated_cs_list = []
        for i in range(len(energy)):
            interpolated_cs = interpolate.splev(energy[i], tck, der=0) *1e-27 # mb--> 1e-27 cm^2. #gives interpolated cross section
            interpolated_cs_list.append(interpolated_cs)
            interpolated_unc_cs = interpolate.splev(energy[i], sigma_tck, der=0) * 1e-27
            flux_weighted_average_cross_section[i] = np.trapezoid(flux[i]*interpolated_cs, energy[i])/np.trapezoid(flux[i],energy[i])
            unc_flux_weighted_average_cross_section[i] = np.trapezoid(flux[i] * interpolated_unc_cs, energy[i])/np.trapezoid(flux[i],energy[i])    
            # plt.plot(energy[i], interpolated_cs*1e27)
        if isotope=='63ZN':
            plt.plot(mon_energy, mon_cs, label='monitor')
        return flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section
    

wa = WeightedAverageFlux('TaSn_stack_55MeV_fluxes.csv', 'TaSn_stack_30MeV_fluxes.csv', 'TaSn_stack_55MeV.csv','TaSn_stack_30MeV.csv')
wa.plot_monitor_reaction(element='Cu', isotope='63ZN', label = 'recommended data 63Zn')
wa.plot_distributions('Cu')