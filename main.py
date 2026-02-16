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
from beamcurrent import BeamCurrent
from activity import Acitivity
from assemble import *
# from activity import Acitivity
# from areal_density import ArealDensity
sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
from nuclearanalysistools.Tendl import *
from nuclearanalysistools.Talys import *
from nuclearanalysistools.tools import *
from nuclearanalysistools.findGammas import *
from my_warnings import *


Acitivity().getA0('Ni', '60CO', plot_=True, compartment=None)



# Talys(os.getcwd() + '/talys/Sn').generateTalysFiles(element='Sn', projectile='p', mass='0', energy='60', ldmodel=None, strength=None, astro=False, potential=None, outputfile=None)
# Ni_61Cu_independent()
Ni_56Ni_cumulative()
# Cu_63Zn_monitor()
# Cu_62Zn_monitor()
# Cu_56Co_monitor()
# Cu_58Co_monitor()
# Ni_57Ni_monitor()
# Cu_54Mn_independent()
# Cu_64Cu_independent()
# Ta_177W()
# Cu_58Co_monitor()

mon_30MeV = [
    ('Cu', '63ZN'), # Looks good
    ('Cu', '65ZN'), # needs to look good
    ('Cu', '62ZN'), # looks ok
    ('Cu', '58CO'), # only 2-3 last points 
    ('Cu', '56CO'), # not valid for 30 mev
    ('Ni', '57NI') # too high
    ]

mon_55MeV = [
    ('Cu', '63ZN'), # Have taken out 63Zn
    ('Cu', '65ZN'),
    ('Cu', '62ZN'),
    ('Cu', '58CO'),
    ('Cu', '56CO'),
    ('Ni', '57NI')
    ]

"""
stack_fluxes_55 = 'TaSn_stack_55MeV_dp_1.010%_fluxes.csv'; full_stack_55 = 'TaSn_stack_55MeV_dp_1.010%.csv'
"""


flux_30_after = 'TaSn_stack_30MeV_dp_1.010%_fluxes.csv' ; full_stack_30 = 'TaSn_stack_30MeV_dp_1.010%.csv'
flux_55_after = 'TaSn_stack_55MeV_dp_1.020%_fluxes.csv' ; full_stack_55 = 'TaSn_stack_55MeV_dp_1.010%.csv'
flux_30_before = 'TaSn_stack_30MeV_dp_1.000%_fluxes.csv'; full_stack_30 = 'TaSn_stack_30MeV_dp_1.000%.csv'
flux_55_before = 'TaSn_stack_55MeV_dp_1.000%_fluxes.csv'; full_stack_55 = 'TaSn_stack_55MeV_dp_1.000%.csv'
wa_55_before = WeightedAverageFlux(stack_fluxes = flux_55_before, full_stack = full_stack_55, stack = 'stack_55_MeV')
wa_55_after = WeightedAverageFlux(stack_fluxes = flux_55_after, full_stack = full_stack_55, stack = 'stack_55_MeV')
wa_30_before = WeightedAverageFlux(stack_fluxes = flux_30_before, full_stack = full_stack_30, stack = 'stack_30_MeV')
wa_30_after = WeightedAverageFlux(stack_fluxes = flux_30_after, full_stack = full_stack_30, stack = 'stack_30_MeV')
bc_30_before = BeamCurrent(wa_30_before, stack = 'stack_30_MeV', monitor_reactions=mon_30MeV)
bc_30_after = BeamCurrent(wa_30_after, stack = 'stack_30_MeV', monitor_reactions=mon_30MeV)
bc_55_before = BeamCurrent(wa_55_before, stack = 'stack_55_MeV', monitor_reactions=mon_55MeV)
bc_55_after = BeamCurrent(wa_55_after, stack = 'stack_55_MeV', monitor_reactions=mon_55MeV)

# wa_30_after.plot_flux_distributions('Ta')
# wa_55_after.plot_flux_distributions('Ta')


# bc_30 = bc_30.get_average_and_unc()
# bc_55 = bc_55.get_average_and_unc()
# print(bc_30)



"""
save_bc = os.getcwd() + '/generatedfiles/beamcurrent/figures/' 

bc_30_before.plot_all('30 MeV before variance minimization')
plt.legend()
plt.savefig(save_bc + '30MeV_before.pdf')
plt.show()
bc_30_after.plot_all('30 MeV after variance minimization')
plt.legend()
plt.savefig(save_bc + '30MeV_after.pdf')
plt.show()

bc_55_before.plot_all('55 MeV before variance minimization')
plt.legend()
plt.savefig(save_bc + '55MeV_before.pdf')
plt.show()
bc_55_after.plot_all('55 MeV after variance minimization')
plt.legend()
plt.savefig(save_bc + '55MeV_after.pdf')
plt.show()
"""










# nuclei_cupper = ['65ZN', '63ZN', '62ZN', '60ZN', '64CU', '59CU', '62CU', '61CU', '63NI', '56NI', '57NI', '61CO', '60CO', '58CO', '57CO','56CO', '55CO','59FE', '55FE', '56MN']
# nuclei_nickel = ['64CU', '59CU', '62CU', '61CU', '63NI', '56NI', '57NI', '61CO', '60CO', '58CO', '57CO','56CO', '55CO','59FE', '55FE', '56MN']
# ag = AnalyzeGammas(nuclei_nickel)
# print(ag.findGammasSpecificIsotope('57NI'))
# ag.matchByGamma(gammaLine=379.94, gammaLineTolerance=2.5, minIntensity=None, xrays=False)
# ag.matchByGamma(gammaLine=1046.68, gammaLineTolerance=2.5, minIntensity=None, xrays=False)
# ag.matchByGamma(gammaLine=1757.55, gammaLineTolerance=2.5, minIntensity=None, xrays=False)
# ag.matchByGamma(gammaLine=1919.52, gammaLineTolerance=2.5, minIntensity=None, xrays=False)
# ag.matchByGamma(gammaLine=2804.2, gammaLineTolerance=2.5, minIntensity=None, xrays=False)






# print(ag.orderIsotopesByHalfLife())
# 669.62, 1412.08, 962.06
# varmin= VarianceMini



# varmin.plot_chi_squared('stack_30_MeV', method='p1', compartments=['12'], title_plot='30 MeV stack compartment 12')

# bc_30.plot_beam_current_isotope('Cu', '63ZN', 'blue', '63Zn',     compartments=None)#['08','09'])
# bc_30.plot_beam_current_isotope('Cu', '65ZN', 'red', '65Zn',      compartments=None)#['08','09'])
# bc_30.plot_beam_current_isotope('Cu', '62ZN', 'darkblue', '62Zn', compartments=None)#['08','09'])
# bc_30.plot_beam_current_isotope('Cu', '58CO', 'cyan', '58Co',     compartments=None)#['08','09'])
# bc_30.plot_beam_current_isotope('Cu', '56CO', 'lightpink', '56Co',compartments=None)#['08','09'])
# bc_30.plot_beam_current_isotope('Ni', '57NI', 'hotpink', '57Ni',  compartments=None)#['08','09'])
# plt.legend()
# plt.show()


# bc_30.plot_all()
# bc_55.plot_all()
# plt.show()

# i, di = bc_30.get_average_and_unc()
# print("weighted:")
# print(i)
# average_energy, average_beam_current, average_protons_per_second, average_unc_protons_per_second = bc_30.average_beam_current()
# print("non weighted:")
# print(average_beam_current)
# bc_30.plot_all()
# bc_55.plot_all()
# plt.show()



# stack_fluxes_30_after = 'TaSn_stack_30MeV_dp_1.040%_fluxes.csv'; full_stack_30_after = 'TaSn_stack_30MeV_dp_1.040%.csv'
# stack_fluxes_30_before = 'TaSn_stack_30MeV_dp_1.000%_fluxes.csv'; full_stack_30_before = 'TaSn_stack_30MeV_dp_1.000%.csv'
# wa_30_before = WeightedAverageFlux(stack_fluxes = stack_fluxes_30_before, full_stack = full_stack_30_before, stack = 'stack_30_MeV')
# wa_30_after = WeightedAverageFlux(stack_fluxes = stack_fluxes_30_after, full_stack = full_stack_30_after, stack = 'stack_30_MeV')
# bc_30_before = BeamCurrent(wa_30_before, stack = 'stack_30_MeV')
# bc_30_after = BeamCurrent(wa_30_after, stack = 'stack_30_MeV')
# bc_30_before.plot_all()
# bc_30_after.plot_all()
# plt.show()
