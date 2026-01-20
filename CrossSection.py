
import numpy as np
from get_variables import *
from flux_stack import WeightedAverageFlux
import curie as ci
import matplotlib.pyplot as plt

import sys
sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
from nuclearanalysistools.Tendl import *
from nuclearanalysistools.tools import *

class CrossSection:

    def __init__(self):
        self.wa = WeightedAverageFlux('TaSn_stack_55MeV_fluxes.csv', 'TaSn_stack_30MeV_fluxes.csv', 'TaSn_stack_55MeV.csv','TaSn_stack_30MeV.csv')

    def cross_section(self, element, isotope):
        foils, areal_density, unc_areal_density = areal_density_from_files(element)
        eob_activities, cov_eob_activities = eob_activity_from_files(foils, isotope)
        unc_eob_activity = np.sqrt(cov_eob_activities)
        energy, flux = self.wa.get_flux_energy_stack(element)
        flux_weighted_average_energy, unc_energy_left, unc_energy_right = self.wa.flux_weighted_average_energy(energy, flux)
        # flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section = self.wa.monitor_flux_weighted_average_cross_section(element, isotope)
        protons_per_second = beam_current_from_files()
        decay_constant = ci.Isotope(isotope).decay_const()
        irradiation_time = 3600; unc_irradiation_time = 1  # s
        std_eob_activity = np.where(
            eob_activities > 0,
            (unc_eob_activity / eob_activities)**2,
            0.0)
        eob_activities = eob_activity_from_files(foils, isotope)
        cross_section = eob_activities / (np.array(areal_density) * (1- np.exp (-decay_constant * irradiation_time))  * protons_per_second)
        plt.plot(flux_weighted_average_energy, cross_section[0]*1e27, 'o')
        cu_target = {'Cu63': 0.6915, 'Cu65':0.3085}
        # ta_target = {"Ta181": 1.0}
        # Tendl({"Ta181": 1.0}, 'proton').plotTendl23Unique(productZ='74', productA='177', isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
        Tendl(cu_target, 'proton').plotTendl23Unique(productZ='29', productA='63', isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
        # plt.show()

CrossSection().cross_section('Cu', '63ZN')
wa = WeightedAverageFlux('TaSn_stack_55MeV_fluxes.csv', 'TaSn_stack_30MeV_fluxes.csv', 'TaSn_stack_55MeV.csv','TaSn_stack_30MeV.csv')
wa.plot_monitor_reaction('Cu', '63ZN')
plt.show()