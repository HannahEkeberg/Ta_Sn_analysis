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


class BeamCurrent():
    
    def __init__(self):
        pass

    def areal_density(self, element):
        root = os.getcwd() + '/generatedfiles/arealdensity/'
        for filename in os.listdir(root):
            if element in filename:
                 df = pd.read_csv(root + filename)
        areal_density = df['nuclei/cm2'].values
        unc_areal_density = df['unc nuclei/cm2'].values
        foils = df['foil'].values
        return foils, areal_density, unc_areal_density
            
    def eob_activity(self, foils, isotope):
        root = os.getcwd() + '/generatedfiles/activity/data/'
        eob_activity = np.zeros(len(foils)); cov_eob_activity = np.zeros(len(foils))
        for i, foil in enumerate(foils):
            for filename in os.listdir(root):
                if foil in filename and 'all_isotopes' in filename:
                    df = pd.read_csv(root + filename)
                    df_isotope = df[df['isotope'].astype(str).str.contains(isotope, case=False, na=False)]
                    if df_isotope.empty:
                        eob_activity[i] = 0.0; cov_eob_activity[i]=0.0
                    else:
                        # TODO add if test for array length > 1 
                        if len(df_isotope['fit'].values) > 1:
                            pass
                        else:
                            eob_activity[i] = df_isotope['fit'].values[0]
                            cov_eob_activity[i] = df_isotope['cov'].values[0]
        return eob_activity, cov_eob_activity
    
    def stack(self, element):
        root = os.getcwd() + '/generatedfiles/stack/'
        stack_55_fluxes = pd.read_csv(root + 'TaSn_stack_55MeV_fluxes.csv')
        stack_55 = pd.read_csv(root + 'TaSn_stack_55MeV.csv')
        stack_30_fluxes = pd.read_csv(root + 'TaSn_stack_30MeV_fluxes.csv')
        stack_30 = pd.read_csv(root + 'TaSn_stack_30MeV.csv')

        stack_55 = stack_55[stack_55['compound'] == element]
        stack_30 = stack_30[stack_30['compound'] == element]
        
        E_55 = stack_55['mu_E'].values; E_30 = stack_30['mu_E'].values
        return np.concatenate((E_55, E_30))
    
    def get_flux_energy_stack(self, foils):
        root = os.getcwd() + '/generatedfiles/stack/'
        stack_55_fluxes = pd.read_csv(root + 'TaSn_stack_55MeV_fluxes.csv')
        # generatedfiles/stack/TaSn_stack_55MeV_-25%_fluxes.csv
        # generatedfiles/stack/TaSn_stack_30MeV_-25%.csv
        stack_30_fluxes = pd.read_csv(root + 'TaSn_stack_30MeV_fluxes.csv')
        flux = []; energy = []
        total_stack = pd.concat([stack_55_fluxes, stack_30_fluxes], ignore_index=True)
        for foil in foils:
            print(foil)
            filtered_stack = total_stack[total_stack['name'] == foil]
            E = filtered_stack['energy'].values; F = filtered_stack['flux'].values
            energy.append(E); flux.append(F)
        return energy, flux

    def plot_distributions(self, energy, flux, fwhm, half_max, mean_energy, foils=None):
        # colors = ['dodgerblue', '#a4c483', 'palevioletred', 'darkorange', 'forestgreen', 'orchid', 'mediumpurple', 'navy', 'crimson', 'indianred', 'pink', 'red', 'red','blue']
        colors = Tools().colors()
        if foils == None:
            foils = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14']
        for i in range(len(energy)):
            plt.plot(energy[i], flux[i], color=colors[i], linewidth = 0.7)
            half_flux = np.max(flux[i])/2
            plt.plot(half_max[i], [half_flux, half_flux], color=colors[i], linewidth=0.8, label=foils[i]+ ' - {0:.2f}'.format(mean_energy[i])  + ' MeV (FWHM: {0:.2f})'.format(fwhm[i]))
            plt.vlines(mean_energy[i], ymin=0.0, ymax = np.max(flux[i]), linewidth=0.4, linestyle='--')#, label=r'$\mu=${}'.format(mu))
            plt.legend(fontsize='xx-small')

    def flux_weighted_average_energy(self, energy, flux, plot_distribution=False, foils=None):
        unc_energy_left = np.zeros(len(energy)); unc_energy_right = np.zeros(len(energy))
        fwhm = np.zeros(len(energy)); half_max = []; mean_energy = np.zeros(len(energy))
        # data = []
        for i in range(len(energy)):
            idx = i +1
            # print("Sn "  + str(idx))
            max_flux = np.max(flux[i])
            min_flux = np.min(flux[i])
            half_max_flux = max_flux * 0.5
        
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

        if plot_distribution:
            self.plot_distributions(energy, flux, fwhm, half_max, mean_energy)
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

        
    def flux_weighted_average_cross_section(self, element, foils, isotope):
        # will make reaction integral and uncertainty integral. 
        
        mon_energy, mon_cs, mon_unc_cs, tck, sigma_tck = self.monitor_data(element, isotope)
        energy, flux = self.get_flux_energy_stack(foils)
        flux_weighted_average_cross_section = np.zeros(len(energy))
        unc_flux_weighted_average_cross_section = np.zeros(len(energy))
        interpolated_cs_list = []
        for i in range(len(energy)):
            interpolated_cs = interpolate.splev(energy[i], tck, der=0) *1e-27 # mb--> 1e-27 cm^2. #gives interpolated cross section
            # print(np.mean(energy[i]), np.mean(flux[i]), np.mean(interpolated_cs))
            interpolated_cs_list.append(interpolated_cs)
            interpolated_unc_cs = interpolate.splev(energy[i], sigma_tck, der=0) * 1e-27
            flux_weighted_average_cross_section[i] = np.trapezoid(flux[i]*interpolated_cs, energy[i])/np.trapezoid(flux[i],energy[i])
            unc_flux_weighted_average_cross_section[i] = np.trapezoid(flux[i] * interpolated_unc_cs, energy[i])/np.trapezoid(flux[i],energy[i])    
            # print(flux_weighted_average_cross_section)
            # print(np.trapezoid(flux[i] * interpolated_unc_cs, energy[i]))
            # print(np.trapezoid(flux[i],energy[i])    )

            # plt.plot(energy[i], interpolated_cs*1e27)
        if isotope=='63ZN':
            plt.plot(mon_energy, mon_cs, label='monitor')
        return flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section

    def beamcurrent(self, element, isotope):
        foils, areal_density, unc_areal_density = self.areal_density(element)
        eob_activities, cov_eob_activities = self.eob_activity(foils, isotope)
        # eob_activities*=2500
        energy, flux = self.get_flux_energy_stack(foils)
        flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section = self.flux_weighted_average_cross_section(element, foils, isotope)
        flux_weighted_average_energy, unc_energy_left, unc_energy_right = self.flux_weighted_average_energy(energy, flux)
        
        # plt.plot(flux_weighted_average_energy, flux_weighted_average_cross_section*1e27, '.')
        # plt.show()
        # plt.plot(flux_weighted_average_energy, eob_activities/1e2, '.', label='activities')

        lamb = ci.Isotope(isotope).decay_const()
        irradiation_time = 3600  # s
        protons_per_second = eob_activities / (np.array(areal_density) *  np.exp (-lamb * irradiation_time)) * 1 / flux_weighted_average_cross_section
        beam_current = protons_per_second * elementary_charge * 1e9 # nano ampere
        
        
        def test(beam_current_expected, beam_current, eob_activities, areal_density, lamb, irradiation_time, flux_weighted_average_energy, flux_weighted_average_cross_section):
            # flux_weighted_average_energy
            beam_current_expected = np.ones(len(eob_activities))*beam_current_expected / elementary_charge
            beam_current = beam_current * 1e-9 / elementary_charge
            cross_section = eob_activities / (beam_current_expected*np.array(areal_density)*(1-np.exp(-lamb*irradiation_time)))

            activity = flux_weighted_average_cross_section / (beam_current_expected*np.array(areal_density)*(1-np.exp(-lamb*irradiation_time)))


            cross_section_assuming_wrong_activity = activity / (beam_current*np.array(areal_density)*(1-np.exp(-lamb*irradiation_time)))
            print(cross_section_assuming_wrong_activity*1e27)

            ad = eob_activities / (beam_current_expected*flux_weighted_average_cross_section)*(1-np.exp(-lamb*irradiation_time))
            # plt.plot(flux_weighted_average_energy, cross_section*1e27, '--', label='assuming beam current 100nA')
            # plt.plot(flux_weighted_average_energy, cross_section_assuming_wrong_activity*1e58, '--', label='assuming wrong activity')

            data = []
            # EITHER AD OR A0 is very wrong......
            for i in range(len(cross_section)):
                wabe = flux_weighted_average_energy[i]
                cs_expected = cross_section[i]*1e27; wacs = flux_weighted_average_cross_section[i]*1e27
                perc_diff_cs = cross_section[i]*1e27 /  (flux_weighted_average_cross_section[i]*1e27)
                bc_expected = beam_current_expected[i]*1e9*elementary_charge; bc = beam_current[i]*1e9*elementary_charge
                perc_diff = beam_current_expected[i]*1e9*elementary_charge/(beam_current[i]*1e9*elementary_charge)* 100 
                data.append([wabe, eob_activities[i], activity[i], areal_density[i], ad[i], cs_expected, wacs, perc_diff_cs, bc_expected, bc, perc_diff])
            df = pd.DataFrame(data, columns=['wabe', 'A0','activity to fit', 'ad', 'ad to fit', 'expected cs', 'measured cs', 'diff cs %', 'expected bc', 'measured bc', 'diff bc %'])
            print(df)

        # test(100*1e-9, beam_current, eob_activities, areal_density, lamb, irradiation_time, flux_weighted_average_energy, flux_weighted_average_cross_section)
        data = []
        for i in range(len(beam_current)):
            print(foils[i])
            data.append([foils[i], beam_current[i], eob_activities[i], areal_density[i], flux_weighted_average_energy[i], [flux_weighted_average_cross_section[i]]])
        
        # plt.plot(flux_weighted_average_energy, flux_weighted_average_cross_section*1e27, label=isotope)
        # plt.show()
        df = pd.DataFrame(data, columns=['foil', 'beam current (nA)', 'eob activity (Bq)', 'areal density (p/cm2)', 'flux weighted average energy (MeV)', 'flux weighted average cross section (cm^2)'])
        # print(isotope)
        print(df)

        # plt.plot(flux_weighted_average_energy, beam_current, 'o',label=isotope)
        
        return flux_weighted_average_energy, beam_current, protons_per_second


    def cross_sections(self, element, isotope, protons_second):
        foils, areal_density, unc_areal_density = self.areal_density(element)
        eob_activities, cov_eob_activities = self.eob_activity(foils, isotope)
        energy, flux = self.get_flux_energy_stack(foils)
        flux_weighted_average_energy, unc_energy_left, unc_energy_right = self.flux_weighted_average_energy(energy, flux)
        flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section = self.flux_weighted_average_cross_section(element, foils, isotope)
        beam_current = np.ones(len(energy))*100 #nA
        protons_per_second = beam_current/ elementary_charge * 1e-9 # nano ampere
        protons_per_second_calculated = self.beamcurrent(element, isotope)[2]
        lamb = ci.Isotope(isotope).decay_const()
        irradiation_time = 3600  # s
        cross_section = eob_activities / (np.array(areal_density) *  np.exp (-lamb * irradiation_time)) * 1 / protons_per_second
        cross_section_calculated = eob_activities / (np.array(areal_density) *  np.exp (-lamb * irradiation_time)) * 1 / protons_per_second_calculated
        cross_section_averaged = eob_activities / (np.array(areal_density) *  np.exp (-lamb * irradiation_time)) * 1 / protons_second
        print(cross_section*1e27)
        print(flux_weighted_average_cross_section*1e27)
        print(flux_weighted_average_energy)
        plt.plot(flux_weighted_average_energy, cross_section*1e27, '--', label='cross sections standard bc')
        plt.plot(flux_weighted_average_energy, cross_section*1e27, 'o', label='cross sections standard bc')
        plt.plot(flux_weighted_average_energy,flux_weighted_average_cross_section*1e27, 'o', label='fwacs')
        plt.plot(flux_weighted_average_energy,cross_section_calculated*1e27, 'o', label='with bc calculations')
        plt.plot(flux_weighted_average_energy, cross_section_averaged*1e27, '*', label='with av bc calculations')

# BeamCurrent().areal_density('Cu')
# flux_weighted_average_energy_63Zn, beam_current_63ZN, protons_per_second_63ZN = BeamCurrent().beamcurrent('Cu', '63ZN')





# flux_weighted_average_energy_65Zn, beam_current_65ZN = BeamCurrent().beamcurrent('Cu', '65ZN')
# flux_weighted_average_energy_62Zn, beam_current_62ZN = BeamCurrent().beamcurrent('Cu', '62ZN')

# energy, flux = BeamCurrent().get_flux_energy_stack(['Sn01','Sn02','Sn03','Sn04','Sn05','Sn06','Sn07','Sn08','Sn09','Sn10','Sn11','Sn12','Sn13','Sn014'])
# mean_energy, unc_energy_left, unc_energy_right = BeamCurrent().flux_weighted_average_energy(energy,flux)
# print(mean_energy)


"""
# BeamCurrent().beamcurrent('Cu', '58CO')
# BeamCurrent().beamcurrent('Cu', '56CO')
flux_weighted_average_energy_63Zn, beam_current_63ZN, protons_per_second_63ZN = BeamCurrent().beamcurrent('Cu', '63ZN')
flux_weighted_average_energy_62Zn, beam_current_62ZN, protons_per_second_62ZN = BeamCurrent().beamcurrent('Cu', '62ZN')
flux_weighted_average_energy_65Zn, beam_current_65ZN, protons_per_second_65ZN = BeamCurrent().beamcurrent('Cu', '65ZN')
flux_weighted_average_energy_56CO, beam_current_56CO, protons_per_second_56CO = BeamCurrent().beamcurrent('Cu', '56CO')
flux_weighted_average_energy_58CO, beam_current_58CO, protons_per_second_58CO = BeamCurrent().beamcurrent('Cu', '58CO')


avg_current = np.mean([protons_per_second_63ZN,protons_per_second_62ZN,protons_per_second_65ZN, protons_per_second_56CO, protons_per_second_58CO], axis=0)
"""

# BeamCurrent().cross_sections('Cu', '63ZN', avg_current)
# target, beamParticle
# Tendl(target = {"Cu63": 0.6915,"Cu65": 0.3085}, beamParticle='proton').plotTendl23Unique(productZ='29', productA='65', isomerLevel = None, color=None, lineStyle=None, label='TENDL23', semilog_y=False)

energy, flux = BeamCurrent().get_flux_energy_stack(['Cu01','Cu02','Cu03','Cu04','Cu05','Cu06','Cu07','Cu08','Cu09','Cu10','Cu11','Cu12','Cu13','Cu14']) 
BeamCurrent().flux_weighted_average_energy(energy, flux, plot_distribution=True)

# plt.legend()
# plt.show()


# BeamCurrent().beamcurrent('Ni', '56CO')
# BeamCurrent().monitor_data(element='Cu', isotope='63ZN')
# BeamCurrent().beamcurrent('Cu', '62CU')