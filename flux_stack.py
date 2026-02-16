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

    def __init__(self, stack_fluxes, full_stack, stack, dp=None):
        self.root = os.getcwd() + '/generatedfiles/stack/'
        self.stack_fluxes = pd.read_csv(self.root + stack_fluxes)
        # self.full_stack   = pd.read_csv(self.root + full_stack) # Not actually used
        self.stack        = stack # 55_MeV, 30_MeV
        self.dp           = dp # Whatever var min running... 

    def get_flux_energy_stack(self, element):
        flux = []; energy = []
        foils = self.foils(element, self.stack)
        # total_stack = pd.concat([self.stack_55_fluxes, self.stack_30_fluxes], ignore_index=True)
        
        for foil in foils:
            filtered_stack = self.stack_fluxes[self.stack_fluxes['name'] == foil]
            E = filtered_stack['energy'].values; F = filtered_stack['flux'].values
            energy.append(E); flux.append(F)

        max_index_for_zero_patting = 2
        for i in range(max_index_for_zero_patting):
            if flux[-1][i] > flux[-1][i+1]:
                flux[-1][i] = 0
        return energy, flux
    
    def foils(self, element, stack):
        if stack == 'stack_55_MeV':
            stack_numbs = ['01', '02', '03', '04', '05', '06', '07']
        elif stack == 'stack_30_MeV':
            stack_numbs = ['08', '09','10', '11', '12', '13', '14']
        else:
            raise Exception('Not a valid stack: ' + stack + '. Must be stack_30_MeV or stack_55_MeV')
        # stack_numbs = ['01', '02', '03', '04', '05', '06', '07', '08', '09','10', '11', '12', '13', '14']
        foils = [element + number for number in stack_numbs]
        return foils
    
    def plot_flux_distributions(self, element):
        foils = self.foils(element, self.stack)
        energy, flux = self.get_flux_energy_stack(element)
        colors = Tools().colors()
        mean_energy, unc_energy_left, unc_energy_right = self.flux_weighted_average_energy(energy, flux)
        for i in range(len(energy)):
            plt.plot(energy[i], flux[i], color=colors[i], linewidth = 0.7)
            half_max_flux = np.max(flux[i])/2
            fwhm = unc_energy_left[i] + unc_energy_right[i]
            fwhm_left = mean_energy[i] - unc_energy_left[i]; fwhm_right = mean_energy[i] + unc_energy_right[i]
            plt.vlines(mean_energy[i], ymin=0.0, ymax = np.max(flux[i]), linewidth=0.4, color=colors[i], linestyle='--')#, label=r'$\mu=${}'.format(mu))
            plt.plot([fwhm_left, fwhm_right], [half_max_flux, half_max_flux],  color=colors[i], linewidth=0.8, label=foils[i]+ ' - {0:.2f}'.format(mean_energy[i])  + ' MeV (fwhm: {0:.2f})'.format(fwhm))
            if self.dp:
                plt.title('Stack simulation for: ' + element + ', dp: ' + self.dp + ' - ' +self.stack)
            else:
                stack_split_string = self.stack.split("_")
                plt.title(element + ' - ' + stack_split_string[1] + ' ' + stack_split_string[2])
        # generatedfiles/fluxweightedaverageenergy/stack_30_MeV_Cu_weighted_average_beam_energy.csv
        plt.legend(fontsize='xx-small')
        plt.savefig(os.getcwd() + '/generatedfiles/fluxweightedaverageenergy/' + self.stack + '.pdf')
        plt.show()

    def flux_weighted_average_energy(self, energy, flux, element=None):
        unc_energy_left = np.zeros(len(energy)); unc_energy_right = np.zeros(len(energy))
        fwhm = np.zeros(len(energy)); half_max = []; mean_energy = np.zeros(len(energy))
        for i in range(len(energy)):
           
            def line_interpolation(x, y, i, half):
                return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

            def half_max_flux_energy(E,F):
                half_max_flux = max(F)/2.0
                signs = np.sign(np.add(F, -half_max_flux))  # for each flux, if over half max +1, if under -1, else 0
                zero_crossings = (signs[0:-2] != signs[1:-1]) # find all the points where the flux is over next to something under.
                zero_crossings_i = np.where(zero_crossings)[0]
                return [line_interpolation(E, F, zero_crossings_i[0], half_max_flux),
                        line_interpolation(E, F, zero_crossings_i[-1], half_max_flux)]
            
            half_max_flux_energy = half_max_flux_energy(energy[i], flux[i])
            half_max.append(half_max_flux_energy)
            mean_energy[i] = np.trapezoid(flux[i]*energy[i], energy[i])/np.trapezoid(flux[i],energy[i])
            fwhm[i] = half_max_flux_energy[1]-half_max_flux_energy[0]
            unc_energy_left[i] = mean_energy[i]-half_max_flux_energy[0]; unc_energy_right[i] = half_max_flux_energy[1]-mean_energy[i]   #left and right uncertainty in energy
        
        if element:
            data = []
            for i in range(len(mean_energy)):
                data.append([mean_energy[i], unc_energy_left[i], unc_energy_right[i]])
            df = pd.DataFrame(data, columns=['wabe', 'unc wabe left', 'unc wabe right'])
            
            df.to_csv(os.getcwd () + '/generatedfiles/fluxweightedaverageenergy/' + self.stack + '_' + element + '_weighted_average_beam_energy.csv')

        return mean_energy, unc_energy_left, unc_energy_right
    
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
            monitor_file = 'monitordata/nip57nit/nip57nit.txt'
            useFile=True
        elif element == 'Ni' and 'isotope' == '55CO':
            useFile=False
        # print(element, isotope)
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
    
    def plot_monitor_reaction(self, element, isotope, label=None):
        energy, cs, unc_cs, tck, sigma_tck = self.monitor_data(element, isotope)
        plt.plot(energy, cs, label=label, color='royalblue')
        plt.fill_between(energy, cs+unc_cs, cs-unc_cs, color='blue', alpha=0.1)
        
        energy, flux = self.get_flux_energy_stack(element)
        mean_energy, unc_energy_left, unc_energy_right = self.flux_weighted_average_energy(energy, flux)
        flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section = self.monitor_flux_weighted_average_cross_section(element, isotope)
        """
        plt.errorbar(mean_energy, flux_weighted_average_cross_section*1e27, marker='P', color='peru',linewidth=0.0001,
            xerr=[unc_energy_left, unc_energy_right], yerr=unc_flux_weighted_average_cross_section, elinewidth=1.0, capthick=1.0, capsize=3.0,
            label=None, linestyle='none')
        """
            # label='flux weighted average monitor cross section', linestyle='none')
        plt.title(element + isotope)
        plt.title(element + '(p,x)' + isotope)
        plt.xlabel('Energy (MeV)')
        plt.ylabel('Cross section (mb)')
        plt.xlim(0,80)
        # plt.show()

    def monitor_flux_weighted_average_cross_section(self, element, isotope):
        mon_energy, mon_cs, mon_unc_cs, tck, sigma_tck = self.monitor_data(element, isotope)
        energy, flux = self.get_flux_energy_stack(element)
        flux_weighted_average_cross_section = np.zeros(len(energy))
        unc_flux_weighted_average_cross_section = np.zeros(len(energy))
        interpolated_cs_list = []
        for i in range(len(energy)):
            interpolated_cs = interpolate.splev(energy[i], tck, der=0) *1e-27 # mb--> 1e-27 cm^2. #gives interpolated cross section
            interpolated_cs_list.append(interpolated_cs)
            interpolated_unc_cs = interpolate.splev(energy[i], sigma_tck, der=0) * 1e-27
            flux_weighted_average_cross_section[i] = np.trapezoid(flux[i]*interpolated_cs, energy[i])/np.trapezoid(flux[i],energy[i])
            unc_flux_weighted_average_cross_section[i] = np.trapezoid(flux[i] * interpolated_unc_cs, energy[i])/np.trapezoid(flux[i],energy[i])
            
            
            minimum_cross_section = 1e-28

            flux_weighted_average_cross_section[i] = max(flux_weighted_average_cross_section[i], 0.0)
            unc_flux_weighted_average_cross_section[i] = max(unc_flux_weighted_average_cross_section[i], 0.0)
            if flux_weighted_average_cross_section[i] < minimum_cross_section:
                flux_weighted_average_cross_section[i] = 0.0
                unc_flux_weighted_average_cross_section[i] = 0.0
            
        return flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section
    
"""
wa = WeightedAverageFlux('TaSn_stack_30MeV_dp_1.000%_fluxes.csv', 'TaSn_stack_30MeV_dp_1.000%.csv', 'stack_30_MeV')
wa_2 = WeightedAverageFlux('TaSn_stack_55MeV_dp_1.000%_fluxes.csv', 'TaSn_stack_55MeV_dp_1.000%.csv', 'stack_55_MeV')
wa.plot_flux_distributions('Ni')
wa.plot_monitor_reaction(element='Cu', isotope='62ZN', label=None)
# plt.legend()
wa_2.plot_monitor_reaction(element='Cu', isotope='62ZN', label=None)
plt.show()
"""    