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

class BeamCurrent:

    def __init__(self, weighted_average_class, stack, monitor_reactions):
        self.stack = stack
        self.wa = weighted_average_class
        self.monitor_reactions = monitor_reactions if monitor_reactions else self.get_monitor_reactions()
        # if monitor_reactions:
        #     self.monitor_reactions = monitor_reactions
        # else:
        #     self.monitor_reactions = self.get_monitor_reactions()
    
    def get_monitor_reactions(self):
        return [
            ('Cu', '63ZN'),
            ('Cu', '65ZN'),
            ('Cu', '62ZN'),
            ('Cu', '58CO'),
            ('Cu', '56CO'),
            ('Ni', '57NI'),
            ]

    def beam_current(self, element, isotope):
        foils, areal_density, unc_areal_density = areal_density_from_files(element, self.stack)
        eob_activitiy, std_eob_activity = eob_activity_manually(element, isotope, independent=None, stack=self.stack)
        # eob_activitiy, std_eob_activity = eob_activity_curie(element, isotope, self.stack)
        energy, flux = self.wa.get_flux_energy_stack(element)
        flux_weighted_average_energy, unc_energy_left, unc_energy_right = self.wa.flux_weighted_average_energy(energy, flux)
        flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section = self.wa.monitor_flux_weighted_average_cross_section(element, isotope)
        decay_constant = ci.Isotope(isotope).decay_const()
        irradiation_time = 3600; unc_irradiation_time = 1  # s
        print(eob_activitiy.shape, areal_density.shape, flux_weighted_average_cross_section.shape)
        protons_per_second = eob_activitiy / (np.array(areal_density) * (1- np.exp (-decay_constant * irradiation_time))  * flux_weighted_average_cross_section)
        unc_protons_per_second = protons_per_second * np.sqrt(
              std_eob_activity
            + (unc_areal_density/areal_density)**2 
            + (unc_irradiation_time/irradiation_time)**2 
            + (unc_flux_weighted_average_cross_section/flux_weighted_average_cross_section)**2)
        protons_per_second = np.nan_to_num(protons_per_second, nan=0.0, posinf=0.0, neginf=0.0); unc_protons_per_second = np.nan_to_num(unc_protons_per_second, nan=0.0, posinf=0.0, neginf=0.0)
        beam_current = protons_per_second * elementary_charge * 1e9 # nano ampere
        unc_beam_current = unc_protons_per_second * elementary_charge * 1e9
        data = []
        # for i in range(len(beam_current)):
            # data.append([foils[i], beam_current[i], unc_beam_current[i], protons_per_second[i], unc_protons_per_second[i], eob_activity[i]*1e-6, areal_density[i], flux_weighted_average_energy[i], [flux_weighted_average_cross_section[i]*1e27]])
        # df = pd.DataFrame(data, columns=['foil', 'beam current (nA)', 'unc beam current (nA)', 'protons/s',' unc protons/s', 'eob activity (MBq)', 'areal density (nuclei/cm2)', 'flux weighted average energy (MeV)', 'flux weighted average cross section (mb)'])
        # df.to_csv(save_beam_current_to + 'beam_current_' + element +  '_' + isotope + '.csv')
        return flux_weighted_average_energy, unc_energy_left, unc_energy_right, beam_current, unc_beam_current, protons_per_second, unc_protons_per_second

    def test_eob_activity(self, lamb, areal_density=None, cs=None):
        protons_per_second = 100 / (elementary_charge * 1e9) #nA * 
        beam_current = protons_per_second * elementary_charge * 1e9 # nano ampere
        A0 = protons_per_second * cs * areal_density * (1-np.exp(-lamb * 3600))
        return A0, beam_current, protons_per_second, cs, areal_density

    def plot_beam_current_isotope(self, element, isotope, color, label, compartments=None, remove_zeros=True):
        flux_weighted_average_energy, unc_energy_left, unc_energy_right, beam_current, unc_beam_current, protons_per_second, unc_protons_per_second = self.beam_current(element, isotope)
        indices = self.get_compartment_indices(compartments) if compartments else slice(None)
        if compartments:
            label = element + ' ' + isotope + ' %.1f' %beam_current[indices][0]
        # if beam_current[indices]>0
        if remove_zeros:
            indices = beam_current[indices] > 0
        try:
            plt.errorbar(flux_weighted_average_energy[indices], beam_current[indices], color=color, marker='.', ls = '',linewidth=0.001, xerr=[unc_energy_left[indices], unc_energy_right[indices]], yerr=np.abs(unc_beam_current[indices]), elinewidth=0.5, capthick=0.5, capsize=3.0,label=label)
        except:
            print("not plotting bc for " + element +'_'+ isotope)
        plt.xlabel('Energy MeV')
        plt.ylabel('Beam current nA')

    def get_compartment_indices(self, compartments):
        if isinstance(compartments, str):
            compartments = [compartments]
        if self.stack == 'stack_30_MeV':
            start = '08'
            offset = int(start)
            return [int(c) - offset for c in compartments]
        elif self.stack == 'stack_55_MeV':
            return [int(c) - 1 for c in compartments]
        else:
            return None

    def reactions(self):
        return [
            ('Cu', '63ZN'),
            ('Cu', '62ZN'),
            ('Cu', '65ZN'),
            ('Cu', '58CO'),
            ('Cu', '56CO'),
            ('Ni', '57NI'),
            ]

    def average_beam_current(self, reactions=None, threshold=None):
        if reactions is None:
            reactions = self.reactions()
        beam_currents = []; protons_per_second = []
        for i in range(len(reactions)):
            reaction = reactions[i]; element = reaction[0]; isotope = reaction[1]
            flux_weighted_average_energy, unc_energy_left, unc_energy_right, beam_current, unc_beam_current, protons_per_second, unc_protons_per_second = self.beam_current(element, isotope)
            beam_currents.append(beam_current)
        if threshold is None:
            threshold = 0  # nA
        arr = np.array(beam_currents)
        average_beam_current = [
            col[col >= threshold].mean() if np.any(col >= threshold) else 0.0
            for col in arr.T
        ]
        energies = []
        for element in ['Ni', 'Cu', 'Ta', 'Sn']:
            energy, flux = self.wa.get_flux_energy_stack(element)
            flux_weighted_average_energy, unc_energy_left, unc_energy_right = self.wa.flux_weighted_average_energy(energy, flux)
            energies.append(flux_weighted_average_energy)
        arr = np.array(energies)
        average_energy = arr.mean(axis=0)
        average_protons_per_second = np.array(average_beam_current) / (elementary_charge * 1e9) #nA 
        average_unc_protons_per_second = np.ones(len(average_protons_per_second))*1e-9
        data = []
        """
        for i in range(len(average_protons_per_second)):
            data.append([average_energy[i], average_protons_per_second[i], average_unc_protons_per_second[i]])
        df = pd.DataFrame(data, columns=['wa energy', 'wa protons/s', 'wa unc proton/s'])
        df.to_csv(save_beam_current_to + self.stack + '_weighted_average_beam_current.csv')
        """
        # for i in range(len(average_beam_current)):

        return average_energy, average_beam_current, average_protons_per_second, average_unc_protons_per_second

    def plot_all(self, title=None):
        colors = Tools().colors()
        # self.plot_beam_current_isotope('Cu', '63ZN', color=colors[0], label=r'$^{63}$Zn')
        reactions=self.monitor_reactions
        for i in range(len(reactions)):
            reaction = reactions[i]; element = reaction[0]; isotope = reaction[1]
            if isotope=='62ZN': 
                self.plot_beam_current_isotope('Cu', '62ZN', color=colors[2], label=r'$^{62}$Zn')
            if isotope=='65ZN': 
                self.plot_beam_current_isotope('Cu', '65ZN', color=colors[1], label=r'$^{65}$Zn')
            if isotope=='63ZN': 
                self.plot_beam_current_isotope('Cu', '63ZN', color=colors[0], label=r'$^{63}$Zn')
            if isotope=='58CO': 
                self.plot_beam_current_isotope('Cu', '58CO', color=colors[3], label=r'$^{58}$Co')
            if isotope=='56CO': 
                self.plot_beam_current_isotope('Cu', '56CO', color=colors[4], label=r'$^{56}$Co')
            if isotope=='57NI': 
                self.plot_beam_current_isotope('Ni', '57NI', color=colors[5], label=r'$^{57}$Ni')
        wa_energy, wa_beamcurrent, wa_protons_per_second, wa_unc_protons_per_second = self.average_beam_current(reactions)
        weighted_bc, unc_weighted_bc = self.get_average_and_unc()
        plt.errorbar(wa_energy, weighted_bc, color='darkblue', marker='^', linewidth=0.001, yerr=np.abs(unc_weighted_bc), ls='', elinewidth=0.5, capthick=0.5, capsize=5.0, markersize=5,label='weighted average')
        # plt.errorbar(wa_energy, wa_beamcurrent, color='peru', marker='D', linewidth=0.001, yerr=np.abs(unc_weighted_bc), elinewidth=0.5, capthick=0.5, capsize=5.0, markersize=5,label='average')
        if title:
            plt.title(title)
        # plt.axhline(y=100, color=colors[-1], linestyle='--', linewidth=0.5, label = '100 nA')
        # plt.ylim(0,500)

    def average_value(self, values, unc_values):
        values = np.array(values)
        unc_values = np.array(unc_values)
        mask = (unc_values > 0) & np.isfinite(unc_values) & np.isfinite(values)
        values = values[mask]
        unc_values = unc_values[mask]
        weights = 1 / unc_values**2
        average = np.average(values, weights=weights)
        unc_average = np.sqrt(1 / np.sum(weights))
        return average, unc_average

    def get_average_and_unc(self):
        # flux_weighted_average_energy, unc_energy_left, unc_energy_right, beam_current, unc_beam_current, protons_per_second, unc_protons_per_second 
        reactions=self.monitor_reactions
        I = []; dI = []
        i_bc = 3; i_unc_bc = 4
        for i in range(len(reactions)):
            reaction = reactions[i]; element = reaction[0]; isotope = reaction[1]
            if isotope=='62ZN':
                variable = self.beam_current(element='Cu', isotope='62ZN')
                I.append(variable[i_bc]); dI.append(variable[i_unc_bc])
            if isotope=='65ZN': 
                variable = self.beam_current(element='Cu', isotope='65ZN')
                I.append(variable[i_bc]); dI.append(variable[i_unc_bc])
            if isotope=='63ZN': 
                variable = self.beam_current(element='Cu', isotope='63ZN')
                I.append(variable[i_bc]); dI.append(variable[i_unc_bc])
            if isotope=='58CO': 
                variable = self.beam_current(element='Cu', isotope='58CO')
                I.append(variable[i_bc]); dI.append(variable[i_unc_bc])
            if isotope=='56CO': 
                variable = self.beam_current(element='Cu', isotope='56CO')
                I.append(variable[i_bc]); dI.append(variable[i_unc_bc])
            if isotope=='57NI': 
                variable = self.beam_current(element='Ni', isotope='57NI')
                I.append(variable[i_bc]); dI.append(variable[i_unc_bc])
        I = np.array(I); dI = np.array(dI)
        I_average = []; dI_average = []; percentage_dI_average = []

        for i in range(I.shape[1]):
            av, unc_av = self.average_value(I[:, i], dI[:, i])
            I_average.append(av)
            dI_average.append(unc_av)
            percentage_dI_average.append(unc_av/av*100)
        df = pd.DataFrame({
            'beam current (nA)': I_average,
            'unc beam current (nA)': dI_average,
            'unc beam current (%)': percentage_dI_average
        })

        # df.to_csv(save_beam_current_to + 'beam_current_' + self.stack + '.csv')
        return I_average, dI_average


