import sys
sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
import curie as ci
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.constants import elementary_charge
from get_variables import *
from flux_stack import WeightedAverageFlux
from beamcurrent import BeamCurrent
from activity import Acitivity
from activity_new import Activity_manual
from weighted_average_updated import *
from get_variables import *
from assemble import *
# from activity import Acitivity
# from areal_density import ArealDensity
from nuclearanalysistools.Tendl import *
from nuclearanalysistools.Talys import *
from nuclearanalysistools.tools import *
from nuclearanalysistools.findGammas import *
from my_warnings import *


# Understanding difference between curie and manual activity
# 63Zn, 58Co, 56Co...

def plot_activity_diff(element, isotope, make_manual_A0=False, independent=False):
    if make_manual_A0:
        Activity_manual().one_step_decay(element, isotope, compartment=None, plot_curve=True)
        # Activity_manual().two_step_decay(element, isotope, '177W', compartment=None, plot_curve=True)
        # Acitivity().getA0(element, isotope, plot_= True, compartment=None)
    eob_activity, std_eob_activity = eob_activity_manually(element, isotope, independent=independent)
    curie, std_curie = eob_activity_curie(element, isotope)
    energy, unc_left, unc_right = weighted_average_beam_energy(element)

    plt.errorbar(energy[7:-1],eob_activity[7:-1], yerr=std_eob_activity[7:-1], marker='o',  ls='none',color='peru',label='manual 30 MeV')
    plt.errorbar(energy[7:-1],curie[7:-1], yerr=std_curie[7:-1], marker='o',  ls='none',color='hotpink',label='curie 30 Me')
    plt.errorbar(energy[0:7],eob_activity[0:7], yerr=std_eob_activity[0:7], marker='o',  ls='none',color='dodgerblue',label='manual 55 MeV')
    plt.errorbar(energy[0:7],curie[0:7], yerr=std_curie[0:7], marker='o',  ls='none',color='forestgreen',label='curie 55 MeV')
    plt.legend()
    plt.show()


# plot_activity_diff('Ta', '177TA', independent=False)
# Acitivity().getA0(element='Sn', isotope='111SN', plot_=True, compartment=None)
Ac = Activity_manual(saveIndependent=False)
Ai = Activity_manual(saveIndependent=True)
A = Activity_manual(saveIndependent=None)
Ai_i = Activity_manual(saveIndependent=True, saveIndependentParent=True)
# A.one_step_decay('Ni', '56NI', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# A.one_step_decay('Ni', '56CO', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# Sn_112mIn()


# A.decay_using_daughter_gamma('Ta', '178W', '178TA', compartment=None, plot_curve=True, min_half_lives=24, max_half_lives=None)
# Ai.one_step_decay('Sn', '117SBg', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# Ai.one_step_decay('Ta', '178LU', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# Ac.one_step_decay('Sn', '119SBg', compartment=None, plot_curve=True, min_half_lives=50, max_half_lives=10)
# Ac.one_step_decay('Ta', '177TA', compartment=None, plot_curve=True, max_half_lives=10)
# Ai.one_step_decay('Ta', '179W', compartment=None, plot_curve=True)
# Ai.two_step_decay('Ta', '175TA', '175W', compartment=None, plot_curve=True, max_half_lives=10)









# A.one_step_decay('Cu', '58CO', compartment=None, plot_curve=True)

# Ac.one_step_decay('Ni', '45TI', compartment=None, plot_curve=True,min_half_lives=None, max_half_lives=10)

# Ac.one_step_decay('Ni', '52MN', compartment=None, plot_curve=True, filtering=True, min_half_lives=80, max_half_lives=None)
# Ai.two_step_decay('Ni', '52MN', '52FE', compartment=None, plot_curve=True, max_half_lives=None)


# Ac.one_step_decay('Cu', '57CO', compartment=None, plot_curve=True, max_half_lives=None, max_activity_uncertainty=None)
# Ac.one_step_decay('Ni', '56NI', compartment=None, plot_curve=True, max_half_lives=10, max_activity_uncertainty=None)


# Ai.two_step_decay('Ta', '56CO', '56NI', compartment='04', plot_curve=False, max_half_lives=None)


# CrossSection().cross_section_manual(element='Cu', isotope='64CU', mask=True, independent=False)
# Cu_61Cu_cumulative()
# CrossSection().latex_table('Cu', '61CU', independent=None, stack='stack_55_MeV', A = 61, el='Cu', energy=True)

# Ac.one_step_decay('Cu', '60CU', compartment=None, plot_curve=True, filtering=True, min_half_lives=None, max_half_lives=10)
# Cu_60Cu_cumulative()


# Copper reactions:
# Cu_64Cu_independent()
# Cu_60Cu_cumulative()
# Cu_57Ni_cumulative()

# OK 
# Ac.one_step_decay('Cu', '57NI', compartment=None, plot_curve=True, filtering=True, min_half_lives=None, max_half_lives=10)
# Cu_57Ni_cumulative()

# NOT MADE
# Ac.one_step_decay('Cu', '56NI', compartment=None, plot_curve=True, filtering=True, min_half_lives=None, max_half_lives=10)
# Cu_56Ni_cumulative()

# THIS ONE IS WEIRD - because using 67.412 makes more sense to existing data by Voyles 2021. Did no subtraction...
# Ac.one_step_decay('Cu', '61CO', compartment=None, plot_curve=True, filtering=True, min_half_lives=None, max_half_lives=2)
# Cu_61Co_cumulative()

#Discuss why some points are low. 
# Ac.one_step_decay('Cu', '60CO', compartment=None, plot_curve=True, filtering=True, min_half_lives=4, max_half_lives=10)
# Cu_60Co_cumulative()

# Ac.one_step_decay('Cu', '55CO', compartment=None, plot_curve=True, filtering=True, min_half_lives=None, max_half_lives=10)
# Cu_55Co_cumulative()

#OK
# Ac.one_step_decay('Cu', '59FE', compartment=None, plot_curve=True, filtering=True, min_half_lives=33, max_half_lives=10)
# Cu_59Fe_cumulative()

# OK
# Ac.one_step_decay('Cu', '56MN', compartment=None, plot_curve=True, filtering=True, min_half_lives=None, max_half_lives=10)
# Cu_56Mn_cumulative()

#OK
# Ai.one_step_decay('Cu', '54MN', compartment=None, plot_curve=True, filtering=True, min_half_lives=None, max_half_lives=10)
# Cu_54Mn_independent()



# NICKEL

#OK
# Ai.one_step_decay('Ni', '64CU', compartment=None, plot_curve=True, filtering=True, min_half_lives=None, max_half_lives=10)
# Ni_64Cu_independent()

# OK
# Ai.one_step_decay('Ni', '61CU', compartment=None, plot_curve=True, filtering=False, min_half_lives=None, max_half_lives=10)
# Ni_61Cu_independent()

# OK
# Ai.one_step_decay('Ni', '54MN', compartment=None, plot_curve=True, filtering=True, min_half_lives=None, max_half_lives=10)
# Ni_54Mn_independent()

# OK 
# Ac.one_step_decay('Ni', '52MN', compartment=None, plot_curve=True, filtering=True, min_half_lives=None, max_half_lives=10)
# Ni_52Mn_cumulative()

# Missing some points...
# Ac.one_step_decay('Ni', '51CR', compartment=None, plot_curve=True, filtering=False, min_half_lives=None, max_half_lives=None)
# Ni_51Cr_cumulative()

#Cupper reactions
# Cu_64Cu_independent()
# Cu_61Cu_cumulative()
# Cu_61Cu_monitor()
# Ni_56Ni_cumulative()
# Cu_55Co_cumulative()
# Cu_54Mn_independent()

# Cu_58Co_independent()
# Cu_58mCo_independent()

# Cu_57Co_cumulative()
# Cu_57Co_independent()

# Cu_56Co_independent() # no point since Cu_56Ni is not really produced. totalxs will just be the monitor reactions

# Ni_64Cu_independent()
# Ni_61Cu_independent()
# Ni_60Cu_independent() # FALSE
# Ni_63N # NO , Ni_57Ni No
# Ni_56Ni_cumulative()
# Ni_61Co_cumulative() # FALSE
# Ni_60Co_cumulative() # FALSE
# Ni_58Co_independent() # DID NOT WORK... either for isomer or goudnstate... 

# Ni_58Co_cumulative()

# Ni_57Co_independent()

# Ai.two_step_decay('Ni', '57CO', '57NI', compartment=None, plot_curve=True, min_half_lives=350, max_half_lives=10)

# Ni_57Co_monitor()

# Ni_56Co_cumulative() # Removed all gammas except strongest independent 1771 keV
# Ni_56Co_independent() 
# Ni_55Co_cumulative()
# 55Fe, 59Fe not found... 56Mn
# Ni_54Mn_independent() # I think we need to remove last measurement in foil 4 because too large...
# Ni_52Mn_cumulative() # cannot see anything of 52Fe or 52mMn. Only using gammas 80 hours after eob, total cumulative. 
# Ni_52mMn_cumulative() # Not possible to make any independent measurements. No data early enough to detect :(
# Ni_52Fe Not possible, not even two step decay - just led to negative eob activities.... 
# 51 Mn not observed
# Ni_51Cr_cumulative() # FERDIG


# tendl_reactions(title='ta')
# Ta_tests(Z='72', A='179', isomerLevel='46', label='isomer')
# Sn_tests(Z='50', A='109', isomerLevel=None, label='gs')

# Sn_tests(Z='51', A='116', isomerLevel='03', label='isomer')
# Sn_tests(Z='51', A='116', isomerLevel=None, label='gs', color='yellow')
# # Sn_tests(Z='49', A='113', isomerLevel='01', label='in m', color='red')
# plt.legend()
# plt.show()









# Ai.one_step_decay('Sn', '109SN',compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=50)
# Ai.two_step_decay('Sn', '109IN', '109SN',compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=None)

# str = f"$^{{nat}}$Ni(p,x)$^{{{61}}}$Cu"
# print(str)

# Acitivity().getA0(element='Ta', isotope='181W', plot_=True, compartment=None)
# Ai.one_step_decay('Ta', '181W', compartment=None, plot_curve=True, min_half_lives=40, max_half_lives=None)
# Ta_181W_independent() # Only using foils 11,12,13,14 - below threshold. Once happy with cross section data points, remove!   FERDIGFERDIG 




# testing secondary neutrons....

# Ai.one_step_decay('Ta', '182TAg', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# ta182 = eob_activity_manually(element='Ta', isotope='182TAg', independent=True, stack=None)[0]
# w181 = eob_activity_manually(element='Ta', isotope='181W', independent=True, stack=None)[0]
# print(ta182/w181)

# Ai.one_step_decay('Cu', '64CUg', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=None)


# Ai.two_step_decay('Ta', '179Wg', '179Wm', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=None)
# Ai.one_step_decay('Ta', '179Wm', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# Ac.one_step_decay('Ta', '179Wg', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# Ta_179W_cumulative() # FERDIG
# Ta_179mW_independent() # most likely false.......



# Sn_119Sb_cumulative()


# Acitivity().getA0(element='Ni', isotope='57NI', plot_=True, compartment=None)

#### MYSTERY:
# A.decay_using_daughter_gamma('Ta', '178W', '178TA', compartment=None, plot_curve=True, min_half_lives=5, max_half_lives=None)
# Ai.one_step_decay('Ta', '178W', compartment=None, plot_curve=True, min_half_lives=30, max_half_lives=None)
# Ta_178W_independent()  # Se på Shahid paper... Sjekk hvilke gammaer han har brukt! 
# Ac.one_step_decay('Ta', '178TAm', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10) # max_half_lives = 565 hours
# Ta_178Ta_cumulative()

# Ac.one_step_decay('Ta', '178TAm1', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10) # max_half_lives = 565 hours
# Ta_178mTa_cumulative()


# Ta_178Ta_cumulative()


# Ac.one_step_decay('Ta', '177TA', compartment=None, plot_curve=True, min_half_lives=80, max_half_lives=10) # max_half_lives = 565 hours
# Ai.one_step_decay('Ta', '177W', compartment=None, plot_curve=True, min_half_lives=2.6, max_half_lives=10) # max_half_lives = 22 hours
# Ta_177W_independent() # FERDIG!!!
# Ta_177Ta_cumulative() # FERDIG!!!!
# Ta_177Ta_independent() # FERDIG!!!!



# Ta_176W, 175, 174, 173W----- has no gammas or x-rays


# Ai.one_step_decay('Ta', '180TA', compartment=None, plot_curve=True, min_half_lives=55, max_half_lives=10)
# Ai.one_step_decay('Ta', '180TA', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# Acitivity().getA0(element='Ta', isotope='180TA', plot_=True, compartment=None)
#CONTINUE AFTER 175TA, 180mHF
# Ai.one_step_decay('Ta', '180TA', compartment=None, plot_curve=True, min_half_lives=12, max_half_lives=10, max_activity_uncertainty=5)
# Ta_180Ta_independent() # Not looking great, especially for foils 1-3. Curie and I disagree....

# NOT SEEN
# Ac.one_step_decay('Ta', '179TA', compartment=None, plot_curve=True, max_half_lives=None)
# Ta_179Ta_cumulative() # nothing..   If we get this try and subtract feeding as well. 


# Acitivity().getA0(element='Ta', isotope='176TA', plot_=True, compartment=None)
# Ac.one_step_decay('Ta', '176TA', compartment=None, plot_curve=True, min_half_lives=15, max_half_lives=10)
# Ai_i.two_step_decay('Ta', '176TA', '176W', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# Ta_176Ta_cumulative() # FINISH
# Ta_176Ta_independent() # FINISH!
# Ta_176W_independent() # FINISH!

# Acitivity().getA0(element='Ta', isotope='175TA', plot_=True, compartment=None)
# Ac.one_step_decay('Ta', '175TA', compartment=None, plot_curve=True, min_half_lives=5, max_half_lives=10)
# Ai_i.two_step_decay('Ta', '175TA', '175W', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# did not work to do two step decay. Only reporting 175Ta cumulative
# Ta_175Ta_cumulative() # look good! 
# Ac.one_step_decay('Ta', '175HF', compartment=None, plot_curve=True, min_half_lives=10, max_half_lives=10)
# Ta_175Hf_cumulative() # curie dont wanna plot, had to manually set decay constant. But most likely not observed....
# Ta_175Hf_independent()

# Acitivity().getA0(element='Ta', isotope='174TA', plot_=True, compartment=None)
# Ac.one_step_decay('Ta', '174TA', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# Ai_i.two_step_decay('Ta', '175TA', '175W', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# did not work to do two step decay. Only reporting 175Ta cumulative
# Ta_174Ta_cumulative() # look good! 



# Et lite notat på beta feeding og populering av ground state vs isomer. Sjekk spinn og paritetsforskjell. Sjekk om tilstanden finnes på level scheme til daugher. Hvis ikke, da populeres den ikke. 

# NOTHING TO REPORT HERE
# Ta_173Ta_cumulative()
# Ac.one_step_decay('Ta', '174TA', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=None)
# Ac.one_step_decay('Ta', '173TA', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=None)




# Only using gammas:[332.274, 500.697] INDEPENDENT 
# XS looks a little confusing so we need to have a look... 
# Ai.one_step_decay('Ta', '180HFm', compartment=None, plot_curve=True, min_half_lives=5, max_half_lives=10)
# Ta_180mHF_independent() # Looks good. However, newest version of talys does not plot. Shahid experimental data does not agree with our data. Agree with Titarenko. Codes are between both..
# Acitivity().getA0(element='Ta', isotope='180HFm', plot_=True, compartment=None)
# Ta_179Lu_independent() # Most likely false. Non-independent gammas. Weakly fed. cross sections are most likely due to production of 


# just fix figure saves, 
# Ta_172Hf_cumulative() # curie dont wanna plot, had to manually set decay constant. But most likely not observed....
# Ta_173Hf_cumulative() # curie dont wanna plot, had to manually set decay constant. But most likely not observed....


# Ta_176Hf_cumulative() # stable so why function??
#177Hfm2 should maybe be seen? need to figure out if there are gammas for the m2 isomer, lots of gammas for m1, but only seconds half life.
# 178Hf m2 need to check spectrum. should be gammas...

# Ai.one_step_decay('Ta', '180HFm', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=5)
# Ai.one_step_decay('Ta', '178HFm2', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=None) # NONE


# Acitivity().getA0(element='Ta', isotope='179LU', plot_=True, compartment=None)
# Ai.one_step_decay('Ta', '179LU', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# Ta_179Lu_independent()

# Ac.one_step_decay('Ta', '177LU', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# Ta_177Lu_cumulative()
# 179Lu # Elisabeth har..... 




# TEST SN
# Acitivity().getA0(element='Sn', isotope='115SB', plot_=True, compartment=None)
# Sn_115Sb_cumulative()
# Sn_111Sn_cumulative()
# Ac.one_step_decay('Sn', '113SN',compartment=None, plot_curve=True, min_half_lives=40, max_half_lives=None)
# Sn_113Sn_cumulative()


# Ac.one_step_decay('Sn', '110INm',compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=None)
# Ac.one_step_decay('Sn', '110IN',compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=None)
# Ai.two_step_decay('Sn', '110IN', '110INm', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=None)
# Sn_110In_cumulative()
# Sn_110In_independent()
# Sn_110mIn_cumulative()


# FOR 175W --> 175TA --> 175HF, please double check that we do not have a point for 175Ta in Ta01. 
# Ac.one_step_decay('Ta', '175HF', compartment=None, plot_curve=True, min_half_lives=10, max_half_lives=10)
# Acitivity().getA0(element='Ta', isotope='175TA', plot_=True, compartment=None)
# Only using independent gammas. Especially in Ta01, remove all spectra bfore 100 hours as feeding from 175Ta might be affecting and feeding. 
# [[343.4, 433.0]]. Adding in more in Foils Ta05-08 because of no feeding. 
# Ai.two_step_decay('Ta', '175HF', '175TA', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=None)
# [343.4, 433.0]
# Ta_175Hf_cumulative()
# Ta_175Hf_independent()


# Ac.one_step_decay('Ta', '178HFm2', compartment=None, plot_curve=True, min_half_lives=60, max_half_lives=10)
# Ac.one_step_decay('Ta', '179HFm2', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# Ta_178m2Hf_cumulative()

# Ac.one_step_decay('Ta', '173HF', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=None)
# Ta_173Hf_cumulative() # cumulative, but realistically no feeding... 

 
 # NOPE
# Ac.one_step_decay('Ta', '172HF', compartment=None, plot_curve=True, min_half_lives=12, max_half_lives=None)
# Ta_172Hf_cumulative()
# Ac.one_step_decay('Ta', '179LU', compartment=None, plot_curve=True, min_half_lives=5, max_half_lives=10)
# Ta_179Lu_independent()

# Ta_Lu()



# CrossSection().latex_table('Cu', '65ZN', independent=True, stack='stack_30_MeV', A = 65, el='Zn')
# CrossSection().latex_table('Cu', '63ZN', independent=True, stack='stack_55_MeV', A = 63, el='Zn')
# CrossSection().latex_table('Cu', '62ZN', independent=True, stack='stack_55_MeV', A = 62, el='Co')
# CrossSection().latex_table('Cu', '56CO', independent=False, stack='stack_55_MeV', A = 56, el='Co')
# CrossSection().latex_table('Cu', '58CO', independent=False, stack='stack_55_MeV', A = 58, el='Co')
# s = CrossSection().latex_table('Cu', '60CO', independent=False, stack='stack_55_MeV', A = 60, el='Co', energy=False)


# s = CrossSection().latex_table_cross_section('Cu', '61CO', independent=False, stack='stack_55_MeV', A = 61, el='Co')
# print(s)








def print_latex_cu(stack, element):
    cs = CrossSection()
    s = ''
    s += rf'\hline'
    s += cs.latex_table_energy(stack, element)
    # s += rf'\hline'
    s += "\n" + r"\noalign{\vskip 0.3cm}"
    s += "\n" + r"\multicolumn{8}{c}{\underline{Monitor reactions}} \\[0.2cm]" + "\n"
    s += cs.latex_table_cross_section(element, '65ZN', independent=None, stack=stack, A = 65, el='Zn', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '63ZN', independent=None, stack=stack, A = 63, el='Zn', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '62ZN', independent=None, stack=stack, A = 62, el='Zn', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '58CO', independent=None, stack=stack, A = 58, el='Co', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '56CO', independent=None, stack=stack, A = 56, el='Co', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '61CU', independent=None, stack=stack, A = 61, el='Cu', subscript_independent=False)
    s += "\n" + r"\noalign{\vskip 0.1cm}"
    s += rf'\hline'
    s += "\n" + r"\noalign{\vskip 0.1cm}"
    s += cs.latex_table_cross_section(element, '64CU', independent=True, stack=stack, A = 64, el='Cu' , subscript_independent=True)
    s += cs.latex_table_cross_section(element, '60CU', independent=False, stack=stack, A = 60, el='Cu', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '57NI', independent=False, stack=stack, A = 57, el='Ni', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '61CO', independent=False, stack=stack, A = 61, el='Co', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '60CO', independent=False, stack=stack, A = 60, el='Co', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '55CO', independent=False, stack=stack, A = 55, el='Co', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '59FE', independent=False, stack=stack, A = 59, el='Fe', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '56MN', independent=False, stack=stack, A = 56, el='Mn', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '54MN', independent=True, stack=stack, A = 54, el='Mn', subscript_independent=True)
    return s

def print_latex_ni(stack, element):
    cs = CrossSection()
    s = ''
    s += rf'\hline'
    s += cs.latex_table_energy(stack, element)
    s += "\n" + r"\noalign{\vskip 0.3cm}"
    s += "\n" + r"\multicolumn{8}{c}{\underline{Monitor reactions}} \\[0.2cm]" + "\n"

    s += cs.latex_table_cross_section(element, '57NI', independent=None, stack=stack, A = 57, el='Ni', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '56NI', independent=None, stack=stack, A = 56, el='Ni', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '58CO', independent=None, stack=stack, A = 58, el='Co', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '57CO', independent=None, stack=stack, A = 57, el='Co', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '56CO', independent=None, stack=stack, A = 56, el='Co', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '55CO', independent=None, stack=stack, A = 55, el='Co', subscript_independent=False)
    s += "\n" + r"\noalign{\vskip 0.1cm}"
    s += rf'\hline'
    s += "\n" + r"\noalign{\vskip 0.1cm}"
    s += cs.latex_table_cross_section(element, '64CU', independent=True, stack=stack, A = 64, el='Cu', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '61CU', independent=True, stack=stack, A = 61, el='Cu', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '57CO', independent=True, stack=stack, A = 57, el='Co', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '56CO', independent=True, stack=stack, A = 56, el='Co', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '54MN', independent=True, stack=stack, A = 54, el='Mn', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '52MN', independent=False, stack=stack, A = 52, el='Mn', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '51CR', independent=False, stack=stack, A = 51, el='Cr', subscript_independent=False)
    return s


def print_latex_ta(stack, element):
    cs = CrossSection()
    s = ''
    s += rf'\hline \\'
    s += cs.latex_table_energy(stack, element)
    s += rf'\hline \\'
    s += cs.latex_table_cross_section(element, '181W', independent=True, stack=stack, A = 181, el='W', isomer_state='g', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '179W', independent=False, stack=stack, A = 179, el='W', isomer_state='m+g', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '178W', independent=True, stack=stack, A = 178, el='W', isomer_state='g', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '177W', independent=True, stack=stack, A = 178, el='W', isomer_state='g', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '176W', independent=True, stack=stack, A = 176, el='W', isomer_state='g', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '180Ta', independent=True, stack=stack, A = 180, el='Ta', isomer_state='g', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '178TAm1', independent=False, stack=stack, A = 178, el='Ta', isomer_state='m1', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '178TA', independent=False, stack=stack, A = 178, el='Ta', isomer_state='g', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '177TA', independent=False, stack=stack, A = 177, el='Ta', isomer_state='g', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '177TA', independent=True, stack=stack, A = 177, el='Ta', isomer_state='g', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '176TA', independent=False, stack=stack, A = 176, el='Ta', isomer_state='g', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '176TA', independent=True, stack=stack, A = 176, el='Ta', isomer_state='g', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '175TA', independent=False, stack=stack, A = 175, el='Ta', isomer_state='g', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '180HFm', independent=True, stack=stack, A = 180, el='Hf', isomer_state='m', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '175Hf', independent=False, stack=stack, A = 175, el='Hf', isomer_state='g', subscript_independent=False)
    s += cs.latex_table_cross_section(element, '175Hf', independent=True, stack=stack, A = 175, el='Hf', isomer_state='g', subscript_independent=True)
    s += cs.latex_table_cross_section(element, '173Hf', independent=False, stack=stack, A = 173, el='Hf', isomer_state='g', subscript_independent=False)
    return s



# s = print_latex_cu('stack_55_MeV', 'Cu')
# s = print_latex_cu('stack_30_MeV', 'Cu')
# s = print_latex_ni('stack_55_MeV', 'Ni')
# s = print_latex_ni('stack_30_MeV', 'Ni')
# s = print_latex_ta('stack_55_MeV', 'Ta')
# for line in s.split(r'\\'):
#     print(line + r'\\')




# Ni_56Co_independent()
# Ni_56Co_cumulative()
# Ni_56Ni_cumulative()
# Ni_56Co_independent()
# Ni_58Co_cumulative()
# Ni_58Co_independent()
# Ni_58mCo_independent()
# Ai.two_step_decay('Ta', '176TA', '176W', compartment=None, plot_curve=True, max_half_lives=10, max_activity_uncertainty=0.50)
# Ac.one_step_decay('Ta', '177HFm2', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# Ta_177m2Hf_cumulative()
# Ac.one_step_decay('Ta', '179HFm2', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# Ta_179m2Hf_cumulative()
# Ta_175Ta_cumulative()
# Ac.one_step_decay('Ta', '176TA', compartment=None, plot_curve=True, min_half_lives=1, max_half_lives=10)
# Ta_176Ta_cumulative()
# Ta_176Ta_independent()


# Ni_61Co_cumulative()
# Ni_60Co_cumulative()
#.two_step_decay('Ta', '178TA', '178W', compartment=None, plot_curve=True)
# Ni_61Cu_independent()


# Talys(os.getcwd() + '/talys/Sn').generateTalysFiles(element='Sn', projectile='p', mass='0', energy='60', ldmodel=None, strength=None, astro=False, potential=None, outputfile=None)
# Ni_61Cu_independent()
# Ni_56Ni_cumulative()
# Ta_180Ta_independent()
# Ta_177Ta_independent()
# Ta_177W_independent()
# Ta_177Ta_cumulative()
# Cu_63Zn_monitor()
# Cu_65Zn_monitor()
# Cu_62Zn_monitor()
# Cu_56Co_monitor()
# Cu_58Co_monitor()
# Cu_58Co_independent()
# Cu_58mCo_independent()



# Cu_56Mn_cumulative()
# Cu_59Fe_cumulative()

# Cu_54Mn_independent()
# Cu_64Cu_independent()
# Ta_177W()
# Cu_58Co_monitor()








################ MONITOR REACTIONS ###################

# A.one_step_decay('Cu', '62ZN', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# A.one_step_decay('Cu', '63ZN', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# A.one_step_decay('Cu', '65ZN', compartment=None, plot_curve=True, min_half_lives=6, max_half_lives=10)
# A.one_step_decay('Cu', '56CO', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# A.one_step_decay('Cu', '58CO', compartment=None, plot_curve=True, min_half_lives=80, max_half_lives=10)
# Ac.one_step_decay('Cu', '61CU', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# A.one_step_decay('Ni', '55CO', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# A.one_step_decay('Ni', '56NI', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# A.one_step_decay('Ni', '56CO', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10) # impossible to pull out 56ni
# A.one_step_decay('Ni', '57NI', compartment=None, plot_curve=True, min_half_lives=None, max_half_lives=10)
# A.one_step_decay('Ni', '57CO', compartment=None, plot_curve=True, min_half_lives=350, max_half_lives=10)
# A.one_step_decay('Ni', '58CO', compartment=None, plot_curve=True, min_half_lives=80, max_half_lives=10)

Cu_62Zn_monitor()
Cu_63Zn_monitor()
Cu_65Zn_monitor()
Cu_61Cu_monitor()
Cu_56Co_monitor()
Cu_58Co_monitor()
Ni_55Co_monitor()
Ni_56Ni_monitor()
Ni_56Co_monitor()
Ni_57Ni_monitor()
Ni_57Co_monitor()
Ni_58Co_monitor()



