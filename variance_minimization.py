from flux_stack import WeightedAverageFlux
from beamcurrent import BeamCurrent
import os
import numpy as np
from get_variables import *
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
import curie as ci
import sys
sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
from nuclearanalysistools.Tendl import *
from nuclearanalysistools.tools import *

import warnings

warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in divide",
    category=RuntimeWarning
)
warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in scalar divide",
    category=RuntimeWarning
)
warnings.filterwarnings(
    "ignore",
    message="divide by zero encountered in divide",
    category=RuntimeWarning
)

# warnings.filterwarnings("error", category=RuntimeWarning)


# dp_array = np.arange(0.80, 1.30, 0.1)
# dp_array = np.arange(0.80, 1.21, 0.001)
# dp_array = dp_array = np.arange(0.80, 1.21, 0.1)
# dp_array = np.arange(0.95, 1.11, 0.01).append(np.arange(0.80, 1.21, 0.1))
dp_array = np.union1d(
    np.arange(0.95, 1.11, 0.01),
    np.arange(0.80, 1.21, 0.1)
)
# dp_array = np.arange(0.80, 1.21, 0.1)
colors = Tools().colors()

class VarianceMinimization:

    def __init__(self, dp_array, monitor_reactions):
        self.dp_array = dp_array
        if monitor_reactions:
            self.monitor_reactions = monitor_reactions
        else:
            self.monitor_reactions = self.get_monitor_reactions()
     
    def get_monitor_reactions(self):
        return [
            ('Cu', '63ZN'),
            ('Cu', '65ZN'),
            ('Cu', '62ZN'),
            ('Cu', '58CO'),
            ('Cu', '56CO'),
            ('Ni', '57NI'),
            ]

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
    
    def p1(self, x, a, b):
        x = np.array(x)
        b_array = np.zeros(len(x))
        b_array.fill(b)
        return a*x + b_array
    
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

    def fit_p1(self, x_data, y_data, unc_data):
        y_data = np.array(y_data)
        y_data[np.isnan(y_data)] = 0
        popt, cov  = curve_fit(self.p1, x_data, y_data, p0=[0, 100], sigma=unc_data)
        return popt[0], popt[1]

    def chi_squared(self, expected, observed, unc_observed):
        observed = np.array(observed)
        expected = np.array(expected)
        chi2 = np.sum((observed - expected)**2 / unc_observed**2)
        return chi2
    
    def get_data_chi_squared(self, x_data, y_data, unc_y_data, method):
        x = np.array(x_data)
        if method == 'p0':
            print('p0')
            print(len(y_data))
            true = self.fit_p0(x_data, y_data, unc_y_data)
            true_array = np.zeros(len(y_data))
            true_array.fill(true)
            dgf = 1
        elif method == 'p1':
            print('p1')
            a,b = self.fit_p1(x_data, y_data, unc_y_data)
            b_array = np.zeros(len(x))
            b_array.fill(b)
            true_array = a*x+b_array
            dgf = 2
        return x, y_data, true_array, unc_y_data, dgf

    def run_chi_squared(self, x_data, y_data, unc_y_data, method):
        x, y_data, true_array, unc_y_data, dgf = self.get_data_chi_squared(x_data, y_data, unc_y_data, method)
        chi2 = self.chi_squared(y_data, true_array, unc_y_data)
        red_chi2 = chi2/(len(y_data)-dgf)
        return chi2, red_chi2

    def get_data_cross_section(self, element, isotope, stack, wa, bc, compartments, plot=False):
        foils, areal_density, unc_areal_density = areal_density_from_files(element, stack)
        indices = self.get_compartment_indices(foils, compartments) if compartments else slice(None)
        eob_activities, cov_eob_activities = eob_activity_from_files(foils, isotope)
        energy, flux = wa.get_flux_energy_stack(element)
        flux_weighted_average_energy, unc_energy_left, unc_energy_right = wa.flux_weighted_average_energy(energy, flux)
        flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section = wa.monitor_flux_weighted_average_cross_section(element, isotope)
        I_average, dI_average = bc.get_average_and_unc()
        protons_per_second_transformed = np.array(I_average)*1e-9/elementary_charge
        unc_protons_per_second_transformed = np.array(dI_average)*1e-9/elementary_charge
        lamb = ci.Isotope(isotope).decay_const()
        irradiation_time = 3600  # s
        cross_section = eob_activities / (np.array(areal_density) *  (1-np.exp (-lamb * irradiation_time))) * 1 / protons_per_second_transformed
        unc_irradiation_time = 0.3
        unc_eob_activity = np.sqrt(cov_eob_activities)
        std_eob_activity = np.where(
            eob_activities > 0,
            (unc_eob_activity / eob_activities)**2,
            0.0)
        unc_cross_section = cross_section * np.sqrt( 
                std_eob_activity
                + (unc_areal_density/areal_density)**2 
                + (unc_irradiation_time/irradiation_time)**2  
                + (unc_protons_per_second_transformed/protons_per_second_transformed)** 2
                # + add unc decay constant
            )
        monitor_cross_section, unc_monitor_cross_section = wa.monitor_flux_weighted_average_cross_section(element, isotope)
        
        if plot:
            cross_section*=1e27 ; flux_weighted_average_cross_section*=1e27
            ratio = flux_weighted_average_cross_section/cross_section
            # print('Ratio for ' + isotope + ' - ' + str(ratio))
            unc_cross_section = np.ones(len(cross_section))*1e-27  # TODO Must fix
            wa.plot_monitor_reaction(element, isotope, label='IAEA recommended data')
            plt.errorbar(flux_weighted_average_energy, cross_section, marker='P', color='forestgreen',linewidth=0.0001,
            xerr=[unc_energy_left, unc_energy_right], yerr=unc_cross_section, elinewidth=1.0, capthick=1.0, capsize=3.0,
            label='flux weighted average cross section', linestyle='none')
            plt.legend()
            plt.show()
        return flux_weighted_average_energy[indices], cross_section[indices], unc_cross_section[indices], monitor_cross_section[indices], unc_monitor_cross_section[indices]

    def get_compartment_indices(self, all_foils, compartments):
        indices = [i for i, foil in enumerate(all_foils) if any(foil.endswith(c) for c in compartments)]
        if not indices:
            raise Exception("Compartments and foils do not match, wrong stack?")
        return indices
    
    def float(self, x):
        return float(np.asarray(x).squeeze())

    def get_data_beam_current(self, element, isotope, stack, wa, bc, compartments):
        foils, areal_density, unc_areal_density = areal_density_from_files(element, stack)
        indices = self.get_compartment_indices(foils, compartments) if compartments else slice(None)
        energy, flux = wa.get_flux_energy_stack(element)
        flux_weighted_average_energy, unc_energy_left, unc_energy_right = wa.flux_weighted_average_energy(energy, flux)
        flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section = wa.monitor_flux_weighted_average_cross_section(element, isotope)
        flux_weighted_average_energy, unc_energy_left, unc_energy_right, beam_current, unc_beam_current, protons_per_second, unc_protons_per_second = bc.beam_current(element, isotope)
        print(isotope, beam_current[indices])
        return flux_weighted_average_energy[indices], beam_current[indices], unc_beam_current[indices]

    def plot_chi2(self, dp_array, chi2_array, smooth_curve=True, title=None, label=None, color='saddlebrown'):
        label = label if label else r'reduced $\chi^2$'
        if smooth_curve:
            dp_smooth = np.linspace(dp_array.min(), dp_array.max(), 500)
            tck = interpolate.splrep(dp_array, chi2_array, s=0)
            chi2_smooth_d_zero = interpolate.splev(dp_smooth, tck)
            plt.plot(dp_smooth, chi2_smooth_d_zero, label=label, color=color)
            dp_min = dp_smooth[np.argmin(chi2_smooth_d_zero)]
            print("Smooth curve minimum: "  + str(dp_min))
        else:
            plt.plot(dp_array, chi2_array, color=color, label=label, alpha=0.5)
        plt.plot(dp_array, chi2_array, 'o', color=color, label=r'reduced $\chi^2 data$')
        min_index_data = chi2_array.index(min(chi2_array))
        min_chi2_data = chi2_array[min_index_data]
        min_dp_data = dp_array[min_index_data]
        label_vertical_line = r'$\chi^2_\text{min}$=%.2f' %min_chi2_data +  ' at '+ str(min_dp_data) 
        plt.axvline(min_dp_data, color=color, linestyle='--', linewidth=0.5, label=label_vertical_line) 
        plt.legend()
        if title:
            plt.title(title)
        else:
            plt.title(r'Reduced $\chi^2$')
            # plt.show()
            # sigma_tck = interpolate.splrep(energy, unc_cs, s=0)

    def plot_chi2_energy(self, energies, chi2, label=None):
        label = label if label else r'reduced $\chi^2$'
        plt.plot(energies, chi2, color='lightpink', label=r'reduced $\chi^2$')
        plt.plot(energies, chi2, 'o', color='dodgerblue', label=r'reduced $\chi^2 data$')

    def get_beam_current_data(self, stack, wa, bc, compartments):
        energies = []; beam_currents = []; unc_beam_currents = []
        for isotope in self.monitor_reactions:
            element = isotope[0]; isotope=isotope[1]
            wabe, beam_current, unc_beam_current = self.get_data_beam_current(element, isotope, stack, wa, bc, compartments)
            energies.append(wabe)
            beam_currents.append(beam_current)
            unc_beam_currents.append(unc_beam_current)

        energies = np.asarray(energies)
        beam_currents = np.asarray(beam_currents)
        unc_beam_currents = np.asarray(unc_beam_currents)
        mask = (
            np.isfinite(energies) &
            np.isfinite(beam_currents) &
            np.isfinite(unc_beam_currents) &
            (unc_beam_currents > 0)
        )
        energies = energies[mask]
        beam_currents = beam_currents[mask]
        unc_beam_currents = unc_beam_currents[mask]
        return energies, beam_currents, unc_beam_currents

    def get_chi2_values(self, stack, wa, bc, compartments, method):
        energies, beam_currents, unc_beam_currents = self.get_beam_current_data(stack, wa, bc, compartments)
        # energies = []; beam_currents = []; unc_beam_currents = []
        # for isotope in self.monitor_reactions:
        #     element = isotope[0]; isotope=isotope[1]
        #     wabe, beam_current, unc_beam_current = self.get_data_beam_current(element, isotope, stack, wa, bc, compartments)
        #     energies.append(wabe)
        #     beam_currents.append(beam_current)
        #     unc_beam_currents.append(unc_beam_current)

        # energies = np.asarray(energies)
        # beam_currents = np.asarray(beam_currents)
        # unc_beam_currents = np.asarray(unc_beam_currents)
        # mask = (
        #     np.isfinite(energies) &
        #     np.isfinite(beam_currents) &
        #     np.isfinite(unc_beam_currents) &
        #     (unc_beam_currents > 0)
        # )
        # energies = energies[mask]
        # beam_currents = beam_currents[mask]
        # unc_beam_currents = unc_beam_currents[mask]
        chi2, red_chi2 = self.run_chi_squared(energies, beam_currents, unc_beam_currents, method)
        return np.average(energies), chi2, red_chi2

    def plot_beam_currents(self, stack, title_plot=''):
        for dp in self.dp_array:
            full_stack, flux_stack = self.get_files(dp, stack)
            wa = WeightedAverageFlux(flux_stack, full_stack, stack)
            bc = BeamCurrent(wa, stack, self.monitor_reactions)
            # print("dp " + str(dp))
            bc.plot_all()
            plt.show()
            for isotope in self.monitor_reactions:
                element = isotope[0]; isotope=isotope[1]
                self.get_data_cross_section(element, isotope, stack, wa, bc, compartments=None, plot=True)

    def plot_fitted_bc(self, stack, wa, bc, compartments, dp=None):
        energies, beam_currents, unc_beam_currents = self.get_beam_current_data(stack, wa, bc, compartments)
        x, y_data, true_array, unc_y_data, dgf = self.get_data_chi_squared(energies, beam_currents, unc_beam_currents, method='p0')
        x_array = np.linspace(x[0]-1, x[-1]+1, len(true_array))
        plt.plot(x_array, true_array, label='fit p0')
        plt.plot(x_array, true_array, label='fit p0')
        x, y_data, true_array, unc_y_data, dgf = self.get_data_chi_squared(energies, beam_currents, unc_beam_currents, method='p1')
        x_array = np.linspace(x[0]-1, x[-1]+1, len(true_array))
        plt.plot(x_array, true_array, label='fit p1')
        plt.plot(x_array, true_array, label='fit p1')
        color_idx=0
        for i in self.monitor_reactions:
            element = i[0]; isotope= i[1]
            color_idx += 1
            color = colors[color_idx]
            bc.plot_beam_current_isotope(element, isotope, color=color, label=(element + ' ' + isotope), compartments=compartments)
        plt.title('Compartment ' + str(compartments) + ' - dp: ' + str(dp))    
        plt.legend()
        plt.show()

    def plot_chi_squared(self, stack, method, compartments=None, title_plot=''):
        chi2s = []; red_chi2s = []; energies = []; chi2s_2=[]; red_chi2s_2=[]
        for dp in self.dp_array:
            print("dp: "+ str(dp))
            full_stack, flux_stack = self.get_files(dp, stack)
            wa = WeightedAverageFlux(flux_stack, full_stack, stack)
            bc = BeamCurrent(wa, stack, self.monitor_reactions)
            if method == 'p0' or method == 'p1':
                energy, chi2, red_chi2 = self.get_chi2_values(stack, wa, bc, compartments, method)
                chi2s.append(chi2)
                red_chi2s.append(red_chi2)
                energies.append(energy)
            elif method == 'p0 and p1':
                energy, chi2, red_chi2 = self.get_chi2_values(stack, wa, bc, compartments, method='p0')
                chi2s.append(chi2)
                red_chi2s.append(red_chi2)
                energy, chi2, red_chi2 = self.get_chi2_values(stack, wa, bc, compartments, method='p1')
                chi2s_2.append(chi2)
                red_chi2s_2.append(red_chi2)
            elif method =='global_xs':
                energies = []; cross_sections = []; unc_cross_sections = []
                for isotope in self.monitor_reactions:
                    element = isotope[0]; isotope=isotope[1]
                    wabe, cross_section, unc_cross_section, monitor_cross_section, unc_monitor_cross_section = self.get_data_cross_section(element, isotope, stack, wa, bc, compartments)
                    cross_section = cross_section
                    energies.append(wabe)
                    cross_sections.append(cross_section)
                    unc_cross_sections.append(unc_cross_section)
                energies = np.concatenate(energies)
                cross_sections = np.concatenate(cross_sections)
                unc_cross_sections = np.concatenate(unc_cross_sections)
                chi2 = self.chi_squared(monitor_cross_section, cross_section, unc_cross_section)
                dgf_cs = 1
                red_chi2 = chi2/(len(cross_sections)-dgf_cs) 
                chi2s.append(chi2)    
                red_chi2s.append(red_chi2)
            elif method == 'test':
                for isotope in self.monitor_reactions:
                    element = isotope[0]; isotope=isotope[1]
                    # print(dp, isotope) 
                    bc.beam_current(element, isotope)
            elif  method == 'test2':
                self.plot_fitted_bc(stack, wa, bc, compartments, dp)
                plt.show()
            elif method == 'test3':
                for isotope in self.monitor_reactions:
                    element = isotope[0]; isotope=isotope[1]
                    # bc.plot_all()
                    # plt.show()
                    self.get_data_cross_section(element, isotope, stack, wa, bc, compartments, plot=True)


        if method == 'p0 and p1':
            self.plot_chi2(self.dp_array, red_chi2s_2,False, title_plot, 'p1','dodgerblue')
            self.plot_chi2(self.dp_array, red_chi2s,False, title_plot, 'p0','saddlebrown')
            plt.show()
        elif method == 'p0' or method =='p1':
            self.plot_chi2(self.dp_array, red_chi2s,False, title_plot, method,'saddlebrown')
            plt.show()


zn65 = ('Cu', '65ZN'); zn63 = ('Cu', '63ZN'); zn62 = ('Cu', '62ZN')
co58 = ('Cu', '58CO'); co56 = ('Cu', '56CO'); ni57 = ('Ni', '57NI')
varmin = VarianceMinimization(dp_array, [zn65, zn62,zn63, co58, co56, ni57])
# varmin = VarianceMinimization(dp_array, [zn65, zn62,co58, co56, ni57])
# varmin.plot_chi_squared('stack_55_MeV', method='p0 and p1', compartments=['07'])
varmin.plot_chi_squared('stack_30_MeV', method='p0', compartments=['14'])
# varmin.plot_chi_squared('stack_55_MeV', method='test2', compartments=['07'])
# varmin.plot_chi_squared('stack_30_MeV', method='p0 and p1', compartments=['11'])
# varmin.plot_chi_squared('stack_55_MeV', method='p0', compartments=['07'])
# varmin.plot_chi_squared('stack_55_MeV', method='test3')
# varmin.plot_chi_squared('stack_30_MeV', method='test3')
plt.show()

