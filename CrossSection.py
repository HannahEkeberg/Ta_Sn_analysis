
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

path = os.getcwd() + '/generatedfiles/crossections/data/'
path_manual = os.getcwd() + '/generatedfiles/crossections/data_manual/'

class CrossSection:

    def __init__(self):
        pass

    def cross_section(self, element, isotope, mask=True, independent=False):
        foils, areal_density, unc_areal_density = areal_density_from_files(element)
        # eob_activity, std_eob_activity = eob_activity_manually(element, isotope)
        eob_activity, std_eob_activity = eob_activity_curie(element, isotope, independent)
        beam_current, unc_beam_current = weighted_average_beam_current()
        energy, unc_left, unc_right = weighted_average_beam_energy(element)
        decay_constant = ci.Isotope(isotope).decay_const()
        irradiation_time = 3600; unc_irradiation_time = 1  # s
        cross_section = eob_activity / (np.array(areal_density) * (1- np.exp (-decay_constant * irradiation_time))  * beam_current)
        unc_cross_section = cross_section * np.sqrt(
            std_eob_activity
            + (unc_areal_density/areal_density)**2 
            + (unc_irradiation_time/irradiation_time)**2 
            + (unc_beam_current/beam_current)**2
        )
        df = pd.DataFrame({
            "cross_section": cross_section,
            "unc_cross_section": unc_cross_section  
        })
        if independent:
            str = path + element + '_' + isotope + '_ind.csv'
        else:
            str = path + element + '_' + isotope + '_cum.csv'
        df.to_csv(str)
        if mask:
            mask = cross_section>0
            return energy[mask], unc_left[mask], unc_right[mask], cross_section[mask]*1e27, unc_cross_section[mask]*1e27
        else:
            # print(mask)
            return energy, unc_left, unc_right, cross_section*1e27, unc_cross_section*1e27

    def cross_section_manual(self, element, isotope, independent=False):
        foils, areal_density, unc_areal_density = areal_density_from_files(element)
        eob_activity, std_eob_activity = eob_activity_manually(element, isotope, independent)
        # eob_activity, std_eob_activity = eob_activity_curie(element, isotope)
        beam_current, unc_beam_current = weighted_average_beam_current()
        energy, unc_left, unc_right = weighted_average_beam_energy(element)
        decay_constant = ci.Isotope(isotope).decay_const()
        irradiation_time = 3600; unc_irradiation_time = 1  # s
        # print(len(eob_activity), len(beam_current), len(areal_density))
        # print(eob_activity)
        cross_section = eob_activity / (np.array(areal_density) * (1- np.exp (-decay_constant * irradiation_time))  * beam_current)
        unc_cross_section = cross_section * np.sqrt(
            std_eob_activity
            + (unc_areal_density/areal_density)**2 
            + (unc_irradiation_time/irradiation_time)**2 
            + (unc_beam_current/beam_current)**2
        )
        df = pd.DataFrame({
            "energy": energy,
            "unc_left": unc_left,
            "unc_right": unc_right,
            "cross_section": cross_section,
            "unc_cross_section": unc_cross_section  
        })
        # print(cross_section)
        if independent:
            str = path_manual + element + '_' + isotope + '_ind.csv'
        elif independent == False:
            str = path_manual + element + '_' + isotope + '_cum.csv'
        else:
            str = path_manual + element + '_' + isotope + '.csv'
        # print(df)
        df.to_csv(str)
        # print(mask)

    def calculate_cross_section(self, eob_activity, areal_density, decay_constant, irradiation_time, beam_current,
                                std_eob_activity, unc_areal_density, unc_irradiation_time, unc_beam_current):
        cross_section = eob_activity / (np.array(areal_density) * (1- np.exp (-decay_constant * irradiation_time))  * beam_current)
        unc_cross_section = cross_section * np.sqrt(
            std_eob_activity
            + (unc_areal_density/areal_density)**2 
            + (unc_irradiation_time/irradiation_time)**2 
            + (unc_beam_current/beam_current)**2
        )
        return cross_section, unc_cross_section

    def cross_section_subtract(self, element, isotope_parent, isotope_daughter, branching_ratio):
        # idx_energy = 0; idx_unc_left=1; idx_unc_right=2; idx_cs = 3; idx_unc_cs = 4; 
        energy, unc_left, unc_right, cross_section_d, unc_cross_section_d = self.cross_section(element, isotope_daughter, mask=False)
        energy, unc_left, unc_right, cross_section_p, unc_cross_section_p = self.cross_section(element, isotope_parent, mask=False)
        # print(daugher_cumulative)
        # energy = daugher_cumulative[idx_energy]
        # unc_left = daugher_cumulative[idx_unc_left]; unc_right = daugher_cumulative[idx_unc_right]
        daughter_independent = cross_section_d - cross_section_p*branching_ratio
        unc_daughter_independent = unc_cross_section_d - unc_cross_section_p*branching_ratio
        # print(cross_section_d)
        # print(cross_section_p)
        # print(daughter_independent)
        # daughter_independent = daugher_cumulative[idx_cs] - parent_independent[idx_cs]*branching_ratio
        # unc_daughter_independent = daugher_cumulative[idx_unc_cs] - parent_independent[idx_unc_cs]*branching_ratio
        mask = daughter_independent>0
        return energy[mask], unc_left[mask], unc_right[mask], daughter_independent[mask], unc_daughter_independent[mask]
         
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

    def subtract_isomer(self, cumulative_xs, unc_cumulative_xs, subtract_xs, unc_subtract_xs, branching_ratio):
        # independent_xs = cumulative_xs-subtract_xs*branching_ratio
        independent_xs = (cumulative_xs-subtract_xs)/branching_ratio
        unc_independent_xs = np.sqrt(
            unc_cumulative_xs**2 +
            (branching_ratio * unc_subtract_xs)**2)
        return independent_xs*1e27, unc_independent_xs*1e27
    
    def subtract_beta(self, cumulative_xs, unc_cumulative_xs, subtract_xs, unc_subtract_xs, branching_ratio):
        independent_xs = cumulative_xs - subtract_xs*branching_ratio
        unc_branching_ratio = branching_ratio*0.001
        unc_independent_xs = np.sqrt(
            unc_cumulative_xs**2 +
            (branching_ratio * unc_subtract_xs)**2 +
            (subtract_xs * unc_branching_ratio)**2)
        return independent_xs, unc_independent_xs

    def save_subtracted_beta(self, element, isotope, independent, cumulative_xs, unc_cumulative_xs, subtract_xs, unc_subtract_xs, branching_ratio):
        energy, unc_left, unc_right = weighted_average_beam_energy(element)
        independent_xs, unc_independent_xs = self.subtract_beta(cumulative_xs, unc_cumulative_xs, subtract_xs, unc_subtract_xs, branching_ratio)
        df = pd.DataFrame({
            "energy": energy,
            "unc_left": unc_left,
            "unc_right": unc_right,
            "cross_section": independent_xs,
            "unc_cross_section": unc_independent_xs  
        })
        if independent:
            str = path_manual + element + '_' + isotope + '_ind.csv'
        elif independent == False:
            str = path_manual + element + '_' + isotope + '_cum.csv'
        else:
            str = path_manual + element + '_' + isotope + '.csv'
        # print(df)
        df.to_csv(str) 


    def add(self, independent_xs, unc_independent_xs, add_xs, unc_add_xs, branching_ratio):
        cumulative_xs = independent_xs + add_xs*branching_ratio
        unc_cumulative_xs = unc_independent_xs - unc_add_xs*branching_ratio
        return cumulative_xs*1e27, unc_cumulative_xs*1e27
    
    def plot_feeding_corrected(self, xs, unc_xs, energy, unc_left, unc_right):
        mask = xs != 0
        plt.errorbar(energy[mask], xs[mask], marker='P', color='hotpink',linewidth=0.0001,
        xerr=[unc_left[mask], unc_right[mask]], yerr=unc_xs[mask], elinewidth=1.0, capthick=1.0, capsize=3.0,
        label='This work', linestyle='none')

    def plot(self, element, isotope, independent=False):
        energy, unc_left, unc_right, cross_section, unc_cross_section = self.cross_section(element, isotope, mask=True, independent=independent)
        plt.errorbar(energy, cross_section, marker='P', color='darkred',linewidth=0.0001,
        xerr=[unc_left, unc_right], yerr=unc_cross_section, elinewidth=1.0, capthick=1.0, capsize=3.0,
        label='flux weighted average cross section', linestyle='none')

    def save_manual(self, element, isotope, independent=False):
        self.cross_section_manual(element, isotope, independent=independent)

    def plot_manual(self, element, isotope, independent=None, color='dodgerblue', label=None):
        # energy, unc_left, unc_right, cross_section, unc_cross_section = self.cross_section_manual(element, isotope, mask=True, independent=independent)
        energy, unc_left, unc_right, cross_section, unc_cross_section = cross_sections(element, isotope, independent)
        print(cross_section)
        mask = cross_section>0
        if label==None:
            label='This work'
        plt.errorbar(energy[mask], cross_section[mask]*1e27, marker='P', color=color,linewidth=0.0001,
        xerr=[unc_left[mask], unc_right[mask]], yerr=unc_cross_section[mask]*1e27, elinewidth=1.0, capthick=1.0, capsize=3.0,
        label=label, linestyle='none')

    def plot_subtract(self, element, isotope_parent, isotope_daughter, branching_ratio):
        energy, unc_left, unc_right, cross_section, unc_cross_section = self.cross_section_subtract(element, isotope_parent, isotope_daughter, branching_ratio)
        plt.errorbar(energy, cross_section, marker='P', color='darkred',linewidth=0.0001,
        xerr=[unc_left, unc_right], yerr=unc_cross_section, elinewidth=1.0, capthick=1.0, capsize=3.0,
        label='flux weighted average cross section', linestyle='none')

    def latex_table_cross_section(self, element, isotope, A=0, el='', independent=True, stack=None):
        list_str = []
        energy, unc_left, unc_right, cross_section, unc_cross_section = self.cross_section_manual(element, isotope, mask=False, independent=independent)
        if stack:
            start_idx, end_idx = get_indexes_stack(stack)
            energy = energy[start_idx:end_idx]
            unc_left = unc_left[start_idx:end_idx]
            unc_right = unc_right[start_idx:end_idx]
            cross_section = cross_section[start_idx:end_idx]
            unc_cross_section = unc_cross_section[start_idx:end_idx]
        if independent:
            row_name_cs = rf'\makecell{{$^{{{A}}}${el}$_i$}} &'
        else:
            row_name_cs = rf'\makecell{{$^{{{A}}}${el}$_c$}} &'
        for i in range(len(energy)):
            cs = cross_section[i]
            unc_cs = unc_cross_section[i]
            cs = '%.2f' %cs
            unc_cs = '%.2f' %unc_cs
            # print(cs)
            if cross_section[i] != 0:
                s = rf'\makecell{{$ {cs} \pm {unc_cs} $}} & '
            else:
                s = rf'\makecell{{-}} & '
            list_str.append(s)
        row_values = ''.join(list_str)
        total_row = row_name_cs + row_values + r' \\'
        # print(total_row)
        return total_row

    def latex_table_energy(self, stack, element):
        energy, unc_left, unc_right = weighted_average_beam_energy(element)
        if stack:
            start_idx, end_idx = get_indexes_stack(stack)
            energy = energy[start_idx:end_idx]
            unc_left = unc_left[start_idx:end_idx]
            unc_right = unc_right[start_idx:end_idx]
        list_str = []
        row_name_energy = '\makecell{$E_{p}$} &'
        for i in range(len(energy)):
            E = energy[i]
            dE_l = unc_left[i]
            dE_r = unc_right[i]
            E = '%.2f' %E
            dE_l = '%.2f' %dE_l
            dE_r = '%.2f' %dE_r
            s = rf'\makecell{{$ {E}^{{+{dE_l}}}_{{-{dE_r}}} $}} & '
            list_str.append(s)
        row_values = ''.join(list_str)
        total_row = row_name_energy + row_values + r' \\'
        return total_row
