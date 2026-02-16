
import numpy as np
from get_variables import *
from flux_stack import WeightedAverageFlux
from beamcurrent import BeamCurrent
import curie as ci
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge

import sys
sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
from nuclearanalysistools.Tendl import *
from nuclearanalysistools.tools import *

class CrossSection:

    def __init__(self):
        pass

    def cross_section(self, element, isotope):
        foils, areal_density, unc_areal_density = areal_density_from_files(element)
        eob_activitiy, std_eob_activity = eob_activity(element, isotope)
        beam_current, unc_beam_current = weighted_average_beam_current()
        energy, unc_left, unc_right = weighted_average_beam_energy(element)
        decay_constant = ci.Isotope(isotope).decay_const()
        irradiation_time = 3600; unc_irradiation_time = 1  # s
        cross_section = eob_activitiy / (np.array(areal_density) * (1- np.exp (-decay_constant * irradiation_time))  * beam_current)
        unc_cross_section = cross_section * np.sqrt(
            std_eob_activity
            + (unc_areal_density/areal_density)**2 
            + (unc_irradiation_time/irradiation_time)**2 
            + (unc_beam_current/beam_current)**2
        )
        mask = cross_section>0
        return energy[mask], unc_left[mask], unc_right[mask], cross_section[mask]*1e27, unc_cross_section[mask]*1e27
    
    def cross_section_subtract(self, element, isotope_parent, isotope_daughter, branching_ratio):
        idx_energy = 0; idx_unc_left=1; idx_unc_right=2; idx_cs = 3; idx_unc_cs = 4; 
        daugher_cumulative = self.cross_section(element, isotope_daughter)
        parent_independent = self.cross_section(element, isotope_parent)
        energy = daugher_cumulative[idx_energy]
        unc_left = daugher_cumulative[idx_unc_left]; unc_right = daugher_cumulative[idx_unc_right]
        daughter_independent = daugher_cumulative[idx_cs] - parent_independent[idx_cs]*branching_ratio
        unc_daughter_independent = daugher_cumulative[idx_unc_cs] - parent_independent[idx_unc_cs]*branching_ratio
        mask = daughter_independent>0
        return energy[mask], unc_left[mask], unc_right[mask], daughter_independent[mask], unc_daughter_independent[mask]*1e27
         
    def substract_cross_section(self, element, isotope_product, isotope_feeder, branching_ratio):
        #energy, unc_left, unc_right, cross_section, unc_cross_section
        i_cs = 3; i_unc_cs = 4
        product_cumulative = self.cross_section(element, isotope_product)
        product_feeding = self.cross_section(element, isotope_feeder)
        cs_independent = product_cumulative[i_cs] - product_feeding[i_cs] * branching_ratio
        unc_cs_independent = product_cumulative[i_unc_cs] - product_feeding[i_unc_cs] * branching_ratio
        return 
        pass

    def monitor_cross_section(self,element,isotope):
        wa_55, wa_30 = get_wa_from_stack_files()
        wa_55.plot_monitor_reaction(element, isotope, label='IAEA recommended data')
        wa_30.plot_monitor_reaction(element, isotope)

    def plot(self, element, isotope):
        energy, unc_left, unc_right, cross_section, unc_cross_section = self.cross_section(element, isotope)
        plt.errorbar(energy, cross_section, marker='P', color='darkred',linewidth=0.0001,
        xerr=[unc_left, unc_right], yerr=unc_cross_section, elinewidth=1.0, capthick=1.0, capsize=3.0,
        label='flux weighted average cross section', linestyle='none')

    # def plot_subtracted(self, element, isotope_parent, isotope_daughter, branching_ratio):
        # energy[mask], unc_left[mask], unc_right[mask], daughter_independent[mask], unc_daughter_independent[mask]*1e27