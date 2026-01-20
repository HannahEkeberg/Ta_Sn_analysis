import curie as ci
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.constants import elementary_charge
import sys
from get_variables import *
from flux_stack import WeightedAverageFlux
# from activity import Acitivity
# from areal_density import ArealDensity
sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
from nuclearanalysistools.Tendl import *
from nuclearanalysistools.tools import *

save_beam_current_to = './generatedfiles/beamcurrent/'

class BeamCurrent():
    
    def __init__(self, stack_files):
        self.wa = WeightedAverageFlux(stack_files[0], stack_files[1], stack_files[2],stack_files[3])
    
    def beam_current(self, element, isotope):
        foils, areal_density, unc_areal_density = areal_density_from_files(element)
        eob_activities, cov_eob_activities = eob_activity_from_files(foils, isotope)
        unc_eob_activity = np.sqrt(cov_eob_activities)
        energy, flux = self.wa.get_flux_energy_stack(element)
        flux_weighted_average_energy, unc_energy_left, unc_energy_right = self.wa.flux_weighted_average_energy(energy, flux)
        flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section = self.wa.monitor_flux_weighted_average_cross_section(element, isotope)
        decay_constant = ci.Isotope(isotope).decay_const()
        irradiation_time = 3600; unc_irradiation_time = 1  # s
        std_eob_activity = np.where(
            eob_activities > 0,
            (unc_eob_activity / eob_activities)**2,
            0.0)
        protons_per_second = eob_activities / (np.array(areal_density) * (1- np.exp (-decay_constant * irradiation_time))  * flux_weighted_average_cross_section) 
        unc_protons_per_second = protons_per_second * np.sqrt(
              std_eob_activity
            + (unc_areal_density/areal_density)**2 
            + (unc_irradiation_time/irradiation_time)**2 
            + (unc_flux_weighted_average_cross_section/flux_weighted_average_cross_section)**2)
        beam_current = protons_per_second * elementary_charge * 1e9 # nano ampere
        unc_beam_current = unc_protons_per_second * elementary_charge * 1e9
        data = []
        for i in range(len(beam_current)):
            data.append([foils[i], beam_current[i], protons_per_second[i], eob_activities[i], areal_density[i], flux_weighted_average_energy[i], [flux_weighted_average_cross_section[i]]])
        df = pd.DataFrame(data, columns=['foil', 'beam current (nA)', 'protons/s', 'eob activity (Bq)', 'areal density (p/cm2)', 'flux weighted average energy (MeV)', 'flux weighted average cross section (cm^2)'])
        # df.to_csv(save_beam_current_to + 'beam_current_' + element +  '_' + isotope + '.csv')
        return flux_weighted_average_energy, unc_energy_left, unc_energy_right, beam_current, unc_beam_current, protons_per_second, unc_protons_per_second

    def plot_beam_current(self, element, isotope, stack, color, label):
        flux_weighted_average_energy, unc_energy_left, unc_energy_right, beam_current, unc_beam_current, protons_per_second, unc_protons_per_second = self.beam_current(element, isotope)
        if stack == '55':
            start_idx = 0; end_idx = 6
        elif stack == '30':
            start_idx = 6; end_idx = 13
        else:
            plt.errorbar(flux_weighted_average_energy, beam_current, color=color, marker='.', linewidth=0.001, xerr=[unc_energy_left, unc_energy_right], yerr=np.abs(unc_beam_current), elinewidth=0.5, capthick=0.5, capsize=3.0,label=label)
        plt.xlabel('Energy MeV')
        plt.ylabel('Beam current nA')
        # plt.ylim(50,200)

    def weighted_average_beam_current(self, beam_current_files=[]):
        # weighted_average = np.zeros(14)
        root = os.getcwd() + '/generatedfiles/beamcurrent/'
        bc_cu_65Zn = self.get_current_from_csv('beam_current_Cu_65ZN.csv')
        bc_cu_63Zn = self.get_current_from_csv('beam_current_Cu_63ZN.csv')
        bc_cu_62Zn = self.get_current_from_csv('beam_current_Cu_62ZN.csv')
        bc_cu_58Co = self.get_current_from_csv('beam_current_Cu_58CO.csv')
        bc_cu_56Co = self.get_current_from_csv('beam_current_Cu_56CO.csv')
        bc_ni_57Ni = self.get_current_from_csv('beam_current_Ni_57NI.csv')
        
        average_beam_current = np.average([bc_cu_65Zn,bc_cu_63Zn,bc_cu_62Zn,bc_cu_58Co,bc_cu_56Co,bc_ni_57Ni],axis=0)
        pd.DataFrame(average_beam_current).to_csv(root + 'weighted_average_beam_current.csv')
        return average_beam_current

    def get_current_from_csv(self, beam_current_file):
        root = os.getcwd() + '/generatedfiles/beamcurrent/'
        df = pd.read_csv(root + beam_current_file)
        return df['protons/s']#, df['unc protons/s']

    def cross_sections(self, element, isotope, protons_second=None):
        foils, areal_density, unc_areal_density = areal_density_from_files(element)
        # foils, areal_density, unc_areal_density = self.areal_density(element)
        # eob_activities, cov_eob_activities = self.eob_activity(foils, isotope)
        eob_activities, cov_eob_activities = eob_activity_from_files(foils, isotope)
        # energy, flux = self.get_flux_energy_stack(foils)
        energy, flux = self.wa.get_flux_energy_stack(element)
        # flux_weighted_average_energy, unc_energy_left, unc_energy_right = self.flux_weighted_average_energy(energy, flux)
        flux_weighted_average_energy, unc_energy_left, unc_energy_right = self.wa.flux_weighted_average_energy(energy, flux)
        flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section = self.wa.monitor_flux_weighted_average_cross_section(element, isotope)
        # flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section = self.flux_weighted_average_cross_section(element, foils, isotope)
        beam_current = np.ones(len(energy))*100 #nA
        protons_per_second = beam_current/ elementary_charge * 1e-9 # nano ampere
        protons_per_second_calculated = self.beamcurrent(element, isotope)[2]
        lamb = ci.Isotope(isotope).decay_const()
        irradiation_time = 3600  # s
        cross_section = eob_activities / (np.array(areal_density) *  np.exp (-lamb * irradiation_time)) * 1 / protons_per_second
        cross_section_calculated = eob_activities / (np.array(areal_density) *  np.exp (-lamb * irradiation_time)) * 1 / protons_per_second_calculated
        # cross_section_averaged = eob_activities / (np.array(areal_density) *  np.exp (-lamb * irradiation_time)) * 1 / protons_second
        # print(cross_section*1e27)
        # print(flux_weighted_average_cross_section*1e27)
        # print(flux_weighted_average_energy)
        plt.plot(flux_weighted_average_energy, cross_section*1e27, '--', label='cross sections standard bc')
        plt.plot(flux_weighted_average_energy, cross_section*1e27, 'o', label='cross sections standard bc')
        plt.plot(flux_weighted_average_energy,flux_weighted_average_cross_section*1e27, 'o', label='fwacs')
        plt.plot(flux_weighted_average_energy,cross_section_calculated*1e27, 'o', label='with bc calculations')
        # plt.plot(flux_weighted_average_energy, cross_section_averaged*1e27, '*', label='with av bc calculations')



# colors = Tools().colors()
# bc = BeamCurrent(stack_files=None)
# bc.plot_beam_current('Cu','63ZN',  stack=None, color = colors[0], label= r'$^{nat}$Cu(p,x)$^{63}$Zn')
# bc.plot_beam_current('Cu', '65ZN', stack=None, color = colors[1], label= r'$^{nat}$Cu(p,x)$^{65}$Zn')
# bc.plot_beam_current('Cu', '62ZN', stack=None, color = colors[2], label= r'$^{nat}$Cu(p,x)$^{62}$Zn')
# bc.plot_beam_current('Cu', '58CO', stack=None, color = colors[3], label= r'$^{nat}$Cu(p,x)$^{58}$Co')
# bc.plot_beam_current('Cu', '56CO', stack=None, color = colors[4], label= r'$^{nat}$Cu(p,x)$^{56}$Co')
# bc.plot_beam_current('Ni', '57NI', stack=None, color = colors[5], label= r'$^{nat}$Ni(p,x)$^{57}$Ni')
# bc.beam_current('Cu','63ZN')
# bc.beam_current('Cu', '65ZN')
# bc.beam_current('Cu', '62ZN')
# bc.beam_current('Cu', '58CO')
# bc.beam_current('Cu', '56CO')
# bc.beam_current('Ni', '57NI')

# plt.legend()
# plt.show()