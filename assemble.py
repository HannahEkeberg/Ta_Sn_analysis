from crosssection import CrossSection
import os

import sys
sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
from nuclearanalysistools.Tendl import *
from nuclearanalysistools.Exfor import *
from nuclearanalysistools.tools import *
from nuclearanalysistools.Talys import *
import matplotlib.pyplot as plt

from warnings import *

natCu = {"Cu63": 0.6915, "Cu65": 0.3085}
natNi = {"Ni58": 0.680769, "Ni60": 0.262231, "Ni61": 0.011399, "Ni62": 0.036345, "Ni64": 0.009256}
natTa = {"Ta181": 1.0}

exfor_path = os.getcwd() + '/exfor/'
talys_path = os.getcwd() + '/talys/'
tendl_Ta = Tendl(natTa, 'proton')
exfor = Exfor(exfor_path)
talys = Talys(talys_path)
tendl_Cu = Tendl(natCu, 'proton')
tendl_Ni = Tendl(natNi, 'proton')

def Cu_63Zn_monitor():
    element = 'Cu'; isotope = '63ZN'
    CrossSection().plot(element, isotope)
    CrossSection().monitor_cross_section(element,isotope)
    Exfor(exfor_path).plotExforDataFromFilename(filename='Cu_63Zn.txt', setLegend=False)
    tendl_Cu.plotTendl23Unique(productZ='30', productA='63', isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    plt.xlim(0,60)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Cu_62Zn_monitor():
    element = 'Cu'; isotope = '62ZN'
    CrossSection().monitor_cross_section(element,isotope)
    Exfor(exfor_path).plotExforDataFromFilename(filename='Cu_62Zn.txt', setLegend=False, maxE=80)
    tendl_Cu.plotTendl23Unique(productZ='30', productA='62', isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    plt.xlim(0,80)
    plt.ylim(0,200)
    CrossSection().plot(element, isotope)
    legend(fontsize='7')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Cu_65Zn_monitor():
    element = 'Cu'; isotope = '65ZN'
    Exfor(exfor_path).plotExforDataFromFilename(filename='Cu_65Zn.txt', setLegend=False)
    tendl_Cu.plotTendl23Unique(productZ='30', productA='65', isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().monitor_cross_section(element,isotope)
    CrossSection().plot(element, isotope)
    plt.xlim(0,80)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Cu_58Co_monitor():
    element = 'Cu'; isotope = '58CO'
    Exfor(exfor_path).plotExforDataFromFilename(filename='Cu_58Co.txt', setLegend=False)
    CrossSection().monitor_cross_section(element,isotope)
    CrossSection().plot(element, isotope)
    plt.xlim(0,80)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Cu_56Co_monitor():
    element = 'Cu'; isotope = '56CO'
    CrossSection().monitor_cross_section(element,isotope)
    CrossSection().plot(element, isotope)
    Exfor(exfor_path).plotExforDataFromFilename(filename='Cu_56Co.txt', setLegend=False, maxE = 80)
    plt.xlim(0,80)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Ni_57Ni_monitor():
    element = 'Ni'; isotope = '57NI'
    Exfor(exfor_path).plotExforDataFromFilename(filename='Ni_57Ni.txt', setLegend=False, maxE = 80)
    CrossSection().monitor_cross_section(element,isotope)
    CrossSection().plot(element, isotope)
    plt.xlim(0,80)
    plt.ylim(0,300)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Ta_175Ta():
    element = 'Ta'; isotope = '175TA'
    # Exfor(exfor_path).plotExforDataFromFilename(filename='Ni_57Ni.txt', setLegend=False, maxE = 80)
    CrossSection().plot(element, isotope)
    plt.xlim(0,80)
    # plt.ylim(0,300)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Ta_177Ta():
    element = 'Ta'; isotope = '177TA'; max_energy = 80; max_cs=850
    Z = '73'; A = '177'; betaFeeding='beta+'; branchingratio = 1.0; 
    betaPlusDecayChain = {isotope: ['74', branchingratio, None]}
    exfor.plotExforDataFromFilename(filename='Ta_177Ta_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Ta_180Ta():
    element = 'Ta'; isotope = '180TA'; max_energy = 80; max_cs=350
    Z = '73'; A = '180'; betaFeeding=None; branchingratio = None; 
    betaPlusDecayChain = None
    exfor.plotExforDataFromFilename(filename='Ta_180Ta_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def legend(fontsize='small', loc='best'):
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(),fontsize=fontsize, loc=loc)

def set_titles(title=None, element=None, isotope=None, independent=False):
    plt.xlabel('Proton energy (MeV)')
    plt.ylabel('Proton energy (MeV)')
    if title:
        plt.title(title)
    elif element:
        i = 0
        while isotope[i].isdigit():
            i += 1
            A = isotope[:i]; Z = isotope[i:].capitalize()
            str = f"$^{{nat}}${element}(p,x)$^{{{A}}}${Z}"
            if independent:
                str += ' - independent'
            else: 
                str += ' - cumulative'
            plt.title(str)
    else:
        plt.title('')

Ni_57Ni_monitor()
Cu_63Zn_monitor()
Cu_62Zn_monitor()
Cu_65Zn_monitor()  # THIS LOOKS WEIRD
Cu_58Co_monitor()
Cu_56Co_monitor()
# Ta_177Ta()
# Ta_180Ta()
# Sn_119Sb()