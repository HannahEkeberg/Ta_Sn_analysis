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

save_cs_figures_to = os.getcwd() + '/generatedfiles/crossections/figures/'

natCu = {"Cu63": 0.6915, "Cu65": 0.3085}
natNi = {"Ni58": 0.680769, "Ni60": 0.262231, "Ni61": 0.011399, "Ni62": 0.036345, "Ni64": 0.009256}
natTa = {"Ta181": 1.0}
natSn = {"Sn112": 0.0097,"Sn114": 0.0066,"Sn115": 0.0034,"Sn116": 0.1454,"Sn117": 0.0768,"Sn118": 0.2422,"Sn119": 0.0859,"Sn120": 0.3258,"Sn122": 0.0463,"Sn124": 0.0579
}

exfor_path = os.getcwd() + '/exfor/'
talys_path = os.getcwd() + '/talys/'
tendl_Ta = Tendl(natTa, 'proton')
exfor = Exfor(exfor_path)
talys = Talys(talys_path)
tendl_Cu = Tendl(natCu, 'proton')
tendl_Ni = Tendl(natNi, 'proton')
tendl_Sn = Tendl(natSn, 'proton')

def Cu_63Zn_monitor():
    element = 'Cu'; isotope = '63ZN'; max_energy = 80; max_cs= 475
    CrossSection().plot(element, isotope)
    CrossSection().monitor_cross_section(element,isotope)
    Exfor(exfor_path).plotExforDataFromFilename(filename='Cu_63Zn.txt', setLegend=False)
    tendl_Cu.plotTendl23Unique(productZ='30', productA='63', isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    # legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    legend(fontsize='7', loc='best')
    plt.savefig(save_cs_figures_to + 'Cu_63Zn_monitor.pdf')
    plt.show()

def Cu_62Zn_monitor():
    element = 'Cu'; isotope = '62ZN'; max_energy = 80; max_cs= 150
    CrossSection().monitor_cross_section(element,isotope)
    Exfor(exfor_path).plotExforDataFromFilename(filename='Cu_62Zn.txt', setLegend=False, maxE=80)
    tendl_Cu.plotTendl23Unique(productZ='30', productA='62', isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    CrossSection().plot(element, isotope)
    set_titles(element=element, isotope=isotope)
    legend(fontsize='7', loc='best')
    plt.savefig(save_cs_figures_to + 'Cu_62Zn_monitor.pdf')
    plt.show()

def Cu_65Zn_monitor():
    element = 'Cu'; isotope = '65ZN'; max_energy = 60; max_cs= 300
    Exfor(exfor_path).plotExforDataFromFilename(filename='Cu_65Zn.txt', setLegend=False)
    tendl_Cu.plotTendl23Unique(productZ='30', productA='65', isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().monitor_cross_section(element,isotope)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    set_titles(element=element, isotope=isotope)
    legend(fontsize='6', loc='best')
    plt.savefig(save_cs_figures_to + 'Cu_65Zn_monitor.pdf')
    plt.show()

def Cu_58Co_monitor():
    element = 'Cu'; isotope = '58CO'; max_energy = 80; max_cs= 130
    Exfor(exfor_path).plotExforDataFromFilename(filename='Cu_58Co.txt', setLegend=False)
    CrossSection().monitor_cross_section(element,isotope)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='8', loc='best')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Cu_58Co_monitor.pdf')
    plt.show()

def Cu_56Co_monitor():
    element = 'Cu'; isotope = '56CO'; max_energy = 80; max_cs = 16
    CrossSection().monitor_cross_section(element,isotope)
    CrossSection().plot(element, isotope)
    Exfor(exfor_path).plotExforDataFromFilename(filename='Cu_56Co.txt', setLegend=False, maxE = 80)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small', loc='best')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Cu_56Co_monitor.pdf')
    plt.show()

def Ni_57Ni_monitor():
    element = 'Ni'; isotope = '57NI'; max_energy = 80; max_cs = 300
    Exfor(exfor_path).plotExforDataFromFilename(filename='Ni_57Ni.txt', setLegend=False, maxE = 80)
    CrossSection().monitor_cross_section(element,isotope)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='6', loc='best')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Ni_57Ni_monitor.pdf')
    plt.show()

def Cu_61Cu_cumulative():
    element = 'Cu'; isotope = '61CU'; max_energy = 80; max_cs= 250
    Z = '29'; A = '61'; betaFeeding='beta+'; branchingratio = 1.0;
    betaPlusDecayChain = {isotope: ['30', branchingratio, None]}
    exfor.plotExforDataFromFilename(filename='Cu_61Cu_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Cu_64Cu_independent():
    element = 'Cu'; isotope = '64CU'; max_energy = 80; max_cs= 250
    Z = '29'; A = '64'; betaFeeding = None; branchingratio = None;
    betaPlusDecayChain = None
    exfor.plotExforDataFromFilename(filename='Cu_64Cu_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Cu_57Ni_cumulative():
    element = 'Cu'; isotope = '57NI'; max_energy = 80; max_cs= 4
    Z = '28'; A = '57'; betaFeeding = None; branchingratio = None;
    betaPlusDecayChain = None
    exfor.plotExforDataFromFilename(filename='Cu_57Ni_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Cu_60Co_cumulative():
    element = 'Cu'; isotope = '60CO'; max_energy = 80; max_cs= 14
    Z = '27'; A = '60'; betaFeeding = None; branchingratio = None;
    betaPlusDecayChain = None
    exfor.plotExforDataFromFilename(filename='Cu_60Co_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Cu_57Co_cumulative():
    element = 'Cu'; isotope = '57CO'; max_energy = 80; max_cs= 90
    Z = '27'; A = '57'; betaFeeding='beta+'; branchingratio = 1.0;
    betaPlusDecayChain = {isotope: ['28', branchingratio, None]}
    # betaPlusDecayChain = None
    exfor.plotExforDataFromFilename(filename='Cu_57Co_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Cu_57Co_independent():
    pass

def Cu_58Co_independent():
    pass

def Cu_58mCo_independent():
    pass

def Cu_59Fe_cumulative():
    element = 'Cu'; isotope = '59FE'; max_energy = 80; max_cs= 1
    Z = '26'; A = '59'; betaFeeding=None; branchingratio = None;
    betaPlusDecayChain = None#{isotope: ['30', branchingratio, None]}
    exfor.plotExforDataFromFilename(filename='Cu_59Fe_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Cu_56Mn_cumulative():
    element = 'Cu'; isotope = '56MN'; max_energy = 80; max_cs= 1
    Z = '25'; A = '56'; betaFeeding=None; branchingratio = None;
    betaPlusDecayChain = None#{isotope: ['30', branchingratio, None]}
    exfor.plotExforDataFromFilename(filename='Cu_56Mn_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Cu_54Mn_independent():
    element = 'Cu'; isotope = '54MN'; max_energy = 80; max_cs= 10.5
    Z = '25'; A = '54'; betaFeeding=None; branchingratio = None;
    betaPlusDecayChain = None#{isotope: ['30', branchingratio, None]}
    exfor.plotExforDataFromFilename(filename='Cu_54Mn_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Ni_64Cu_independent():
    element = 'Ni'; isotope = '64CU'; max_energy = 80; max_cs= 10.5
    Z = '29'; A = '64'; betaFeeding=None; branchingratio = None;
    betaPlusDecayChain = None 
    # exfor.plotExforDataFromFilename(filename='Cu_54Mn_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Ni_64Cu_independent.pdf')
    plt.show()

def Ni_61Cu_independent():
    element = 'Ni'; isotope = '61CU'; max_energy = 80; max_cs= 14
    Z = '29'; A = '61'; betaFeeding=None; branchingratio = None; decaychain=None
    exfor.plotExforDataFromFilename(filename='Ni_61Cu_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl(element, Z, A)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Ni_61Cu_independent.pdf')
    plt.show()

def Ni_56Ni_cumulative():
    element = 'Ni'; isotope = '56NI'; max_energy = 80; max_cs= 18
    Z = '28'; A = '56'; betaFeeding='beta+'; branchingratio = 1.0; 
    decaychain = {isotope: ['29', branchingratio, None]}
    exfor.plotExforDataFromFilename(filename='Ni_56Ni_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl(element, Z, A, betaPlusDecayChain=decaychain)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Ni_56Ni_cumulative.pdf')
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

def Ta_177W():
    element = 'Ta'; isotope = '177W'; max_energy = 80; max_cs=850
    Z = '74'; A = '177'; betaFeeding=None; branchingratio = None; 
    betaPlusDecayChain = None
    exfor.plotExforDataFromFilename(filename='Ta_177W.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=True)
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

def tendl(element, Z, A, isomerLevel=None, betaPlusDecayChain = None, betaMinusDecayChain = None, isomer_decay_chain = None):
    # betaPlusDecayChain = None, betaMinusDecayChain = None, isomerDecayChain = None
    # {isotope: [productZ, branchingRatio isomerLevel]} #beta+/beta-
    # {isotope: [branchingRatio isomerLevel]} #isomer
    if element == 'Cu':
        tendl = tendl_Cu
    elif element == 'Ni':
        tendl = tendl_Ni
    elif element == 'Ta':
        tendl = tendl_Ta
    elif element == 'Sn':
        tendl = tendl_Sn
    else:
        raise Exception('Invalid element: ' + element)
    if betaPlusDecayChain:
        tendl.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = isomerLevel,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    elif betaMinusDecayChain:
        tendl.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = isomerLevel,  
                                  betaMinusDecayChain = betaMinusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    elif isomer_decay_chain:
        tendl.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = isomerLevel,  
                                  isomer_decay_chain = None,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    else:
        print(tendl)
        tendl.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = isomerLevel, color=None, lineStyle=None, label=None, semilog_y=False)



# Ni_57Ni_monitor()
# Cu_63Zn_monitor()
# Cu_62Zn_monitor()
# Cu_65Zn_monitor()  # THIS LOOKS WEIRD
# Cu_58Co_monitor()
# Cu_56Co_monitor()
# Ta_177Ta()
# Ta_180Ta()
# Sn_119Sb()