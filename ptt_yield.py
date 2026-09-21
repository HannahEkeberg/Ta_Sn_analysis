
import numpy as np
from get_variables import *
import curie as ci
import matplotlib.pyplot as plt
from nuclearanalysistools.Tendl import *

class ThickTargetYield:

    def __init__(self, element, particle):
        self.element = element
        self.particle = particle
        self.Z_beam = 1
        self.irradiation_time = 1

    def calculate_tty(self, isotope, independent, e_max=None, e_min=None):
        energy, cross_section = self.get_data(isotope, independent, e_max, e_min)
        energy = np.flip(energy); cross_section = np.flip(cross_section)
        decay_constant = ci.Isotope(isotope).decay_const('h',True)   #1/h
        smooth_xs = cross_section #* (1.0e-27)
        y = np.empty_like(energy)
        stopping_power = self.stopping_power(energy)
        element = ci.Element(self.element)
        mass = element.mass   
        for i, item in enumerate(cross_section):
            y[i] = np.trapezoid( 
                (np.concatenate(([0, 0], smooth_xs[0:i+1])) / (self.Z_beam*1.602e-19))
                / (np.concatenate((element.S([1e-3, energy[0]-0.2],particle=self.particle, density=1E-3)*1E3 , stopping_power[0:i+1]))),
                x = np.concatenate(([0.1, energy[0]-0.2], energy[0:i+1]))   ) *(6.022e23/mass) / (1E6 * 1E6)
        y = np.concatenate(([0], y))
        if e_min:
            e_range = np.concatenate(([e_min], energy))
        else:
            e_range = np.concatenate(([0], energy))
        # e_range = np.concatenate(([30], energy))
        # e_range = energy
        EOB_TTY = (1-np.exp(-decay_constant[0]*self.irradiation_time))  * y
        saturation_TTY = y
        physical_TTY = decay_constant[0] *y
        return e_range, EOB_TTY

    def calculate_tty_tendl(self, isotope, Z, A):
        energy, cross_section = self.retrieve_by_tendl(Z, A)
        decay_constant = ci.Isotope(isotope).decay_const('h',True)   #1/h
        smooth_xs = cross_section* (1.0e-27)
        # print(len(energy), len(cross_section))
        y = np.empty_like(energy)
        stopping_power = self.stopping_power(energy)
        element = ci.Element(self.element)
        mass = element.mass
        for i, item in enumerate(cross_section):
            y[i] = np.trapezoid( 
                (np.concatenate(([0, 0], smooth_xs[0:i+1])) / (self.Z_beam*1.602e-19))
                / (np.concatenate((element.S([1e-3, energy[0]-0.2],particle=self.particle, density=1E-3)*1E3 , stopping_power[0:i+1]))),
                x = np.concatenate(([0.1, energy[0]-0.2], energy[0:i+1]))   ) *(6.022e23/mass) / (1E6 * 1E6)
        
        y = np.concatenate(([0], y))
        # e_range = np.concatenate(([0], energy))
        e_range = np.concatenate(([30], energy))
        EOB_TTY = (1-np.expm1(-decay_constant[0]*self.irradiation_time))  * y
        saturation_TTY = y
        physical_TTY = decay_constant[0] *y
        return e_range, EOB_TTY

    def plot(self, tty, label, color='black', linestyle=None, linewidth=None, plot_values=False):
        e_range = tty[0]; eob_tty = tty[1]
        if plot_values:
            plt.plot(e_range, eob_tty/0.0036,'o', color=color) # dividing to go from micro Ampere hour to Coulomb
        plt.plot(e_range,eob_tty/0.0036,label=label, linestyle=linestyle, linewidth=linewidth, color=color)

    def stopping_power(self, energy):
        element = ci.Element(self.element) 
        return element.S(energy,particle=self.particle, density=1E-3)*1E3   # MeV/(g/cm^2)
    
    def get_data(self, isotope, independent, e_min=None, e_max=None):
        energy, unc_left, unc_right, xs, unc_xs = cross_sections(self.element, isotope, independent)
        
        e_min = 0; e_max=60
        mask = (energy >= e_min) & (energy <= e_max)
        energy = energy[mask]; xs = xs[mask]
        return energy, xs
    
    
    def retrieve_by_tendl(self, Z, A, Emin, Emax):
        tendl = Tendl({"Ta181": 1.0}, 'proton')
        E, Cs = tendl.tendlData(Z, A, isomerLevel=None, Emin=Emin, Elimit=Emax)
        # plt.plot(E,Cs*1e-27)
        # plt.show()
        return E, Cs
    
    # def plot_purity(self, nuclide_of_interst, other_nuclides, max_cooling_time, minimum_energy=30):
    #     cooling_times = [0, 0.5,1,2,3,4,5,6,7,8,9,10, 11, 12, 13, 14, 15]
    #     if max_cooling_time != None:
    #         extra_cooling_times = np.linspace(16, max_cooling_time, 40)
    #         cooling_times = np.concatenate((cooling_times, extra_cooling_times))

    #     purity_by_time = []
    #     for t in cooling_times:
    #         e_range, tty = self.calculate_tty(nuclide_of_interst, True)
    #         decay_constant = ci.Isotope(nuclide_of_interst).decay_const('h',False)   #1/h
    #         tty_of_interest = tty * np.exp(-decay_constant * t)
    #         tty_total = tty_of_interest.copy()
    #         for i in other_nuclides:
    #             if i == '179W':
    #                 independent=False
    #             else:
    #                 independent=True
    #             e_range_i, tty_i = self.calculate_tty(i, independent)
    #             decay_constant_i = ci.Isotope(i).decay_const('h',False)   #1/h
    #             # yie = tty_i * np.exp(-decay_constant_i * t)
    #             tty_i_cooled = tty_i * np.exp(-decay_constant_i * t)
    #             tty_total += tty_i_cooled

    #         purity = 100 * tty_of_interest / tty_total
    #         purity_by_time.append(purity)

    #     purity_by_time = np.array(purity_by_time)
    #     mask = e_range >= minimum_energy

    #     for energy, purity in zip(e_range[mask], purity_by_time[:, mask].T):
    #         plt.plot(cooling_times, purity, label=f"{energy:.1f} MeV", linewidth=0.9)

    #     plt.xlabel("Cooling time (h)")
    #     plt.ylabel(f"{nuclide_of_interst} purity (%)")
    #     plt.title(f"Radionuclidic purity of {nuclide_of_interst}")
    #     plt.legend()
    #     plt.tight_layout()
    #     plt.show()


            

    
tty = ThickTargetYield('Ta', 'p')
w176 = tty.calculate_tty('176W', True, 0, 50)
w177 = tty.calculate_tty('177W', True)
w178 = tty.calculate_tty('178W', True)
w179 = tty.calculate_tty('179W', False)
w181 = tty.calculate_tty('181W', True)
ta180 = tty.calculate_tty('180TA', True)
ta178m = tty.calculate_tty('178TAm', False)
ta178 = tty.calculate_tty('178TA', False)
ta177 = tty.calculate_tty('177TA', False)
ta176 = tty.calculate_tty('176TA', False)
ta175 = tty.calculate_tty('175TA', False)
hf180m = tty.calculate_tty('180HFm', True)
hf175 = tty.calculate_tty('175HF', False)
hf173 = tty.calculate_tty('173HF', False)
# w181 = tty.calculate_tty('180TA', True)
# w180_t = tty.calculate_tty_tendl('180W', '74', '180', 44, 55)
# w175_t = tty.calculate_tty_tendl('175W', '74', '175', 52, 55)
# w176_t = tty.calculate_tty_tendl('176W', '74', '176', 44, 55)
# w177_t = tty.calculate_tty_tendl('177W', '74', '177', 44, 55)
# w179_t = tty.calculate_tty_tendl('179W', '74', '179', 44, 55)
# w181_t = tty.calculate_tty_tendl('181W', '74', '181', 44, 55)
# w178_t = tty.calculate_tty_tendl('178W', '74', '178', 44, 55)

def plot_wolfram():
    # tty.plot(w177, r'$^{177}$W', color='#741b47', plot_values=True)
    # tty.plot(w178, r'$^{178}$W', color='#f6b26b', linewidth=1.2, linestyle='-.', plot_values=True)
    # tty.plot(w179, r'$^{179}$W', color='#cf9b65', linewidth=1.2, linestyle='-.', plot_values=True)
    # tty.plot(w176, r'$^{176}$W', color='#ac6186', linewidth=1.2, linestyle='-.', plot_values=True)
    # tty.plot(w181, r'$^{181}$W', color='#a6c499', linewidth=1.2, linestyle='-.', plot_values=True)
    
    tty.plot(ta180, r'$^{180}$Ta', color='#741b47', linewidth=1.2, linestyle='-.', plot_values=True)
    tty.plot(ta178m, r'$^{178m}$Ta', color='#d97a5f', linewidth=1.2, linestyle='-.', plot_values=True)
    tty.plot(ta178, r'$^{178}$Ta', color='#f2b66d', linewidth=1.2, linestyle='-.', plot_values=True)
    tty.plot(ta177, r'$^{177}$Ta', color='#b8875b', linewidth=1.2, linestyle='-.', plot_values=True)
    tty.plot(ta176, r'$^{176}$Ta', color='#a85f83', linewidth=1.2, linestyle='-.', plot_values=True)
    tty.plot(ta175, r'$^{175}$Ta', color='#7f6aa8', linewidth=1.2, linestyle='-.', plot_values=True)
    tty.plot(hf180m, r'$^{180m}$Hf', color='#6f9f8f', linewidth=1.2, linestyle='-.', plot_values=True)
    tty.plot(hf175, r'$^{175}$Hf',   color='#a8bf91', linewidth=1.2, linestyle='-.', plot_values=True)
    tty.plot(hf173, r'$^{173}$Hf',   color='#5f7f99', linewidth=1.2, linestyle='-.', plot_values=True)
    # tty.plot(w180_t, r'$^{180}$W (tendl)', color='#6fa8dc', linewidth=1.2, linestyle='--')
    # # tty.plot(w175_t, r'$^{175}$W (tendl)', color='#cf9b65', linewidth=1.2, linestyle='--')

    # tty.plot(w179_t, r'$^{179}$W (tendl)', color='#cf9b65', linewidth=1.2, linestyle='--')
    # tty.plot(w176_t, r'$^{176}$W (tendl)', color='#ac6186', linewidth=1.2, linestyle='--')
    # tty.plot(w177_t, r'$^{177}$W (tendl)', color='#741b47', linewidth=1.2, linestyle='--')
    # tty.plot(w181_t, r'$^{181}$W (tendl)', color='#a6c499', linewidth=1.2, linestyle='--')
    # tty.plot(w178_t, r'$^{178}$W (tendl)', color='#f6b26b', linewidth=1.2, linestyle='--')
    plt.yscale('log')
    plt.xlabel('Proton energy (MeV)')
    plt.ylabel('Physical Thick Target Yield (MBq/C)', fontsize='small')
    plt.legend(fontsize='xx-small')
    plt.savefig('generatedfiles/yields/tty_other.pdf')
    plt.show()

# def plot_purity():


plot_wolfram()

