from crosssection import CrossSection
import os

import sys
sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
from nuclearanalysistools.Tendl import *
from nuclearanalysistools.Exfor import *
from nuclearanalysistools.tools import *
from nuclearanalysistools.Talys import *
import matplotlib.pyplot as plt
from get_variables import *

from warnings import *

save_cs_figures_to = os.getcwd() + '/generatedfiles/crossections/figures/'

natCu = {"Cu63": 0.6915, "Cu65": 0.3085}
natNi = {"Ni58": 0.680769, "Ni60": 0.262231, "Ni61": 0.011399, "Ni62": 0.036345, "Ni64": 0.009256}
natTa = {"Ta181": 0.9998799, "Ta180": 0.0001201}
natSn = {"Sn112": 0.0097,"Sn114": 0.0066,"Sn115": 0.0034,"Sn116": 0.1454,"Sn117": 0.0768,"Sn118": 0.2422,"Sn119": 0.0859,"Sn120": 0.3258,"Sn122": 0.0463,"Sn124": 0.0579
}

exfor_path = os.getcwd() + '/exfor/'
talys_path_old = os.getcwd() + '/talys/'
talys_path = os.getcwd() + '/talys_analysis/'
tendl_Ta = Tendl(natTa, 'proton')
exfor = Exfor(exfor_path)
talys_old = Talys(talys_path_old)
talys = Talys(talys_path)
tendl_Cu = Tendl(natCu, 'proton')
tendl_Ni = Tendl(natNi, 'proton')
tendl_Sn = Tendl(natSn, 'proton')

def Cu_63Zn_monitor():
    element = 'Cu'; isotope = '63ZN'; max_energy = 80; max_cs= 475
    CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope)
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
    CrossSection().plot_manual(element, isotope)
    set_titles(element=element, isotope=isotope)
    legend(fontsize='7', loc='best')
    plt.savefig(save_cs_figures_to + 'Cu_62Zn_monitor.pdf')
    plt.show()

def Cu_65Zn_monitor():
    element = 'Cu'; isotope = '65ZN'; max_energy = 60; max_cs= 300
    Exfor(exfor_path).plotExforDataFromFilename(filename='Cu_65Zn.txt', setLegend=False)
    tendl_Cu.plotTendl23Unique(productZ='30', productA='65', isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().monitor_cross_section(element,isotope)
    CrossSection().plot_manual(element, isotope)
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
    CrossSection().plot_manual(element, isotope, independent=False)
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
    CrossSection().plot_manual(element, isotope)
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
    CrossSection().plot_manual(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='6', loc='best')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Ni_57Ni_monitor.pdf')
    plt.show()

def Cu_61Cu_cumulative():
    element = 'Cu'; isotope = '61CU'; max_energy = 80; max_cs= 300
    Z = '29'; A = '61'; betaFeeding='beta+'; branchingratio = 1.0;
    betaPlusDecayChain = {isotope: ['30', branchingratio, None]}
    exfor.plotExforDataFromFilename(filename='Cu_61Cu_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Cu_61Cu_cumulative.pdf')
    plt.show()

def Cu_60Cu_cumulative():
    element = 'Cu'; isotope = '60CU'; max_energy = 80; max_cs= 50
    Z = '29'; A = '60'; betaFeeding='beta+'; branchingratio = 1.0;
    independent=False
    betaPlusDecayChain = {isotope: ['30', branchingratio, None]}
    exfor.plotExforDataFromFilename(filename='Cu_60Cu_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Cu_60Cu_cumulative.pdf')
    plt.show()

def Cu_62Cu_independent():
    element = 'Cu'; isotope = '62CU'; max_energy = 80; max_cs= 1400
    Z = '29'; A = '62'; betaFeeding='beta+'; branchingratio = 1.0;
    independent=True
    betaPlusDecayChain = {isotope: ['30', branchingratio, None]}
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope, independent=independent)
    # CrossSection().plot_manual(element, isotope='62ZN', independent=independent, color='hotpink')
    # CrossSection().plot_manual(element, isotope='62CU', independent=False, color='darkred')
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Cu_62Cu_independent.pdf')
    plt.show()

def Cu_64Cu_independent():
    element = 'Cu'; isotope = '64CU'; max_energy = 80; max_cs= 250
    Z = '29'; A = '64'; betaFeeding = None; branchingratio = None;
    betaPlusDecayChain = None; independent=True
    exfor.plotExforDataFromFilename(filename='Cu_64Cu_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope,independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Cu_64Cu_independent.pdf')
    plt.show()

def Cu_57Ni_cumulative():
    element = 'Cu'; isotope = '57NI'; max_energy = 80; max_cs= 4
    Z = '28'; A = '57'; betaFeeding = None; branchingratio = None;
    betaPlusDecayChain = None; independent=False
    exfor.plotExforDataFromFilename(filename='Cu_57Ni_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Cu_57Ni_cumulative.pdf')
    plt.show()

def Cu_56Ni_cumulative():
    element = 'Cu'; isotope = '56NI'; max_energy = 80; max_cs= 2
    Z = '28'; A = '56'; betaFeeding = None; branchingratio = None;
    betaPlusDecayChain = None; independent=False
    exfor.plotExforDataFromFilename(filename='Cu_56Ni_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Cu_56Ni_cumulative.pdf')
    plt.show()

def Cu_61Co_cumulative():
    # only using gamma 917
    element = 'Cu'; isotope = '61CO'; max_energy = 80; max_cs=15
    Z = '27'; A = '61'; betaFeeding = 'beta-'; branchingratio = 1.0;
    betaMinusDecayChain = {isotope: ['26', branchingratio, None]}; independent=False
    exfor.plotExforDataFromFilename(filename='Cu_61Co_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaMinusDecayChain = betaMinusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Cu_61Co_cumulative.pdf')
    plt.show()

def Cu_60Co_cumulative():
    # Everything looks good except Cu01. 
    element = 'Cu'; isotope = '60CO'; max_energy = 80; max_cs= 20
    Z = '27'; A = '60'; betaFeeding = None; branchingratio = None;
    betaPlusDecayChain = None
    exfor.plotExforDataFromFilename(filename='Cu_60Co_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope, independent=False)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=False)
    plt.savefig(save_cs_figures_to + 'Cu_60Co_cumulative.pdf')
    plt.show()

def Cu_55Co_cumulative():
    # Everything looks good except Cu01. 
    element = 'Cu'; isotope = '55CO'; max_energy = 80; max_cs= 2
    Z = '27'; A = '55'; betaFeeding = None; branchingratio = None;
    betaPlusDecayChain = None
    exfor.plotExforDataFromFilename(filename='Cu_55Co_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope, independent=False)
    plt.xlim(0,max_energy)
    plt.ylim(-0.1,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=False)
    plt.savefig(save_cs_figures_to + 'Cu_55Co_cumulative.pdf')
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
    CrossSection().plot_manual(element, isotope, independent=False)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=False)
    plt.savefig(save_cs_figures_to + 'Cu_57Co_cumulative.pdf')
    plt.show()

def Cu_57Co_independent():
    # CrossSection().plot(element, isotope)
    element = 'Cu'; isotope = '57CO'; max_energy = 80; max_cs= 90
    Z = '27'; A = '57'; betaFeeding=None; branchingratio = None;
    exfor.plotExforDataFromFilename(filename='Cu_57Co_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)

    CrossSection().plot_manual(element, isotope, independent=True)
    # CrossSection().plot_manual(element, isotope, independent=False, color='hotpink')
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=True)
    plt.savefig(save_cs_figures_to + 'Cu_57Co_independent.pdf')
    plt.show()

def Cu_58Co_independent():
    element = 'Cu'; isotope='58CO'; max_energy=80; max_cs= 80
    Z = '27'; A = '58'; betaFeeding = None; branchingratio = None;
    exfor.plotExforDataFromFilename(filename='Cu_58Co_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = None,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    
    CrossSection().plot_manual(element, isotope, independent=True, color='hotpink')
    CrossSection().plot_manual(element, isotope, independent=False, color='dodgerblue')


    # total, unc_total = cross_sections(element, isotope, independent=False)
    # gs, unc_gs = cross_sections(element, isotope, independent=True)
    # CrossSection().plot_manual(element, isotope, independent=True, color='hotpink')
    # CrossSection().plot_manual(element, isotope, independent=False, color='dodgerblue')
    # br = 0.999988
    # xs, unc_xs = CrossSection().subtract(total, unc_total, gs, unc_gs, br)
    # energy, unc_left, unc_right = weighted_average_beam_energy(element)
    # CrossSection().plot_feeding_corrected(xs, unc_xs, energy, unc_left, unc_right)


    set_titles(element=element, isotope=isotope, independent=True)
    plt.xlim(0,max_energy)
    plt.ylim(-0.1,max_cs)
    legend()
    plt.savefig(save_cs_figures_to + 'Cu_58Co_independent.pdf')
    plt.show()

def Cu_58mCo_independent():
    br = 0.999988
    element = 'Cu'; isotope='58CO'; max_energy=80; max_cs = 80
    Z = '27'; A = '58'; betaFeeding = None; branchingratio = None; isomerLevel='01';
    exfor.plotExforDataFromFilename(filename='Cu_58mCo_independent.txt', setLegend=False, maxE = max_energy)
    energy, unc_left, unc_right = weighted_average_beam_energy(element)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = isomerLevel, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = isomerLevel,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    total, unc_total = cross_sections(element, isotope, independent=False)
    gs, unc_gs = cross_sections(element, isotope, independent=True)
    CrossSection().plot_manual(element, isotope, independent=True, color='hotpink')
    CrossSection().plot_manual(element, isotope, independent=False, color='dodgerblue')
    xs, unc_xs = CrossSection().subtract_isomer(total, unc_total, gs, unc_gs, br)

    CrossSection().plot_feeding_corrected(xs, unc_xs, energy, unc_left, unc_right)
    set_titles(element=element, isotope='58mCO', independent=True)
    # xs2 = CrossSection().add(gs, xs, br)
    # print(xs2*1e27)

    # plt.plot(energy, xs*1e27, '*')
    # plt.plot(energy, xs2*1e27, '*', color='cyan')
    # isomer = total-gs*br
    # unc_isomer = unc_total-unc_gs*br
    # plt.errorbar(energy, isomer, yerr=unc_isomer, xerr=[unc_left, unc_right], marker='o')
    # plt.plot(energy, isomer, 'o')
    # plt.errorbar(energy[0:7], total[0:7], yerr=unc_total[0:7], xerr=[unc_left[0:7], unc_right[0:7]], marker='o', ls='none')
    # plt.errorbar(energy[0:7], gs[0:7], yerr=unc_gs[0:7], xerr=[unc_left[0:7], unc_right[0:7]], marker='o', ls='none')
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend()
    plt.savefig(save_cs_figures_to + 'Cu_58mCo_independent.pdf')
    plt.show()

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
    CrossSection().plot_manual(element, isotope, independent=False)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Cu_59Fe_cumulative.pdf')
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
    CrossSection().plot_manual(element, isotope)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Cu_56Mn_cumulative.pdf')
    plt.show()

def Cu_54Mn_independent():
    element = 'Cu'; isotope = '54MN'; max_energy = 80; max_cs= 10.5
    Z = '25'; A = '54'; betaFeeding=None; branchingratio = None;
    independent=True
    betaPlusDecayChain = None#{isotope: ['30', branchingratio, None]}
    exfor.plotExforDataFromFilename(filename='Cu_54Mn_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Cu.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Cu_54Mn_independent.pdf')
    plt.show()

def Ni_64Cu_independent():
    # cross sections in foils 13,14 looks good. In 1-2, probably false?
    element = 'Ni'; isotope = '64CU'; max_energy = 80; max_cs= 10.5
    Z = '29'; A = '64'; betaFeeding=None; branchingratio = None;
    betaPlusDecayChain = None; independent=True
    # exfor.plotExforDataFromFilename(filename='Cu_54Mn_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ni.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope,independent=independent)
    plt.savefig(save_cs_figures_to + 'Ni_64Cu_independent.pdf')
    plt.show()

def Ni_61Cu_independent():
    element = 'Ni'; isotope = '61CU'; max_energy = 80; max_cs= 14
    Z = '29'; A = '61'; betaFeeding=None; branchingratio = None; decaychain=None
    independent=True
    exfor.plotExforDataFromFilename(filename='Ni_61Cu_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl(element, Z, A)
    # CrossSection().plot(element, isotope, )
    CrossSection().plot_manual(element, isotope, independent=True)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ni_61Cu_independent.pdf')
    plt.show()

def Ni_60Cu_independent():
    # FALSE. Did not see......
    element = 'Ni'; isotope = '60CU'; max_energy = 80; max_cs= 104
    Z = '29'; A = '60'; betaFeeding=None; branchingratio = None; decaychain=None
    independent=True
    exfor.plotExforDataFromFilename(filename='Ni_60Cu_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl(element, Z, A)
    # CrossSection().plot(element, isotope, )
    CrossSection().plot_manual(element, isotope, independent=True)
    # CrossSection().plot_manual(element='Ni', isotope='56CO', independent=False, color='darkred')
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    # plt.savefig(save_cs_figures_to + 'Ni_61Cu_independent.pdf')
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

def Ni_61Co_cumulative():
    # Looks false
    element = 'Ni'; isotope = '61CO'; max_energy = 80; max_cs= 2
    Z = '27'; A = '61'; betaFeeding='beta-'; branchingratio = 1.0; 
    decaychain = {isotope: ['61', branchingratio, None]}
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    CrossSection().plot_manual(element, isotope, independent=False)
    # CrossSection().plot(element, isotope, independent=False)
    tendl(element, Z, A, betaPlusDecayChain=decaychain)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Ni_61Co_cumulative.pdf')
    plt.show()

def Ni_60Co_cumulative():
    element = 'Ni'; isotope = '60CO'; max_energy = 80; max_cs= 150
    Z = '27'; A = '60'; betaFeeding=None; branchingratio = None; 
    decaychain = {isotope: ['60', branchingratio, None]}
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    CrossSection().plot_manual(element, isotope, independent=False)
    # CrossSection().plot(element, isotope, independent=False)
    tendl(element, Z, A, betaPlusDecayChain=decaychain)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Ni_60Co_cumulative.pdf')
    plt.show()

def Ni_57Co_cumulative():
    element = 'Ni'; isotope = '57CO'; max_energy = 80; max_cs= 670
    Z = '27'; A = '57'; betaFeeding='beta+'; branchingratio = 1.0; 
    betaPlusDecayChain = {isotope: ['28', branchingratio, None]}
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ni.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot_manual(element, isotope, independent=False)
    # CrossSection().plot(element, isotope)
    exfor.plotExforDataFromFilename(filename='Ni_57Co_cumulative.txt', setLegend=False, maxE = max_energy)
    # CrossSection().plot(element, isotope, independent=False)



    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Ni_57Co_cumulative.pdf')
    plt.show()

def Ni_57Co_independent():
    element = 'Ni'; isotope = '57CO'; max_energy = 80; max_cs= 580
    Z = '27'; A = '57'; betaFeeding=None; branchingratio = 1.0; 
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ni.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = None,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    exfor.plotExforDataFromFilename(filename='Ni_57Co_independent.txt', setLegend=False, maxE = max_energy)
    # CrossSection().plot_manual(element, isotope, independent=True)
    # CrossSection().plot_manual(element, isotope, independent=False, color='hotpink')
    # CrossSection().plot_manual(element, isotope='57NI', independent=True, color='hotpink')
    # CrossSection().plot(element, isotope, independent=False)

    total, unc_total = cross_sections(element, isotope, independent=False)
    feeding, unc_feeding = cross_sections(element, isotope='57NI', independent=False)
    xs, unc_xs = CrossSection().subtract_beta(total, unc_total, feeding, unc_feeding, branchingratio)
    energy, unc_left, unc_right = weighted_average_beam_energy(element)
    # CrossSection().plot_manual(element, isotope, independent=True, color='hotpink')
    # CrossSection().plot_manual(element, isotope, independent=False, color='dodgerblue')
    CrossSection().plot_feeding_corrected(xs, unc_xs, energy, unc_left, unc_right)

    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=True)
    plt.savefig(save_cs_figures_to + 'Ni_57Co_independent.pdf')
    plt.show()


def Ni_56Co_cumulative():
    element = 'Ni'; isotope = '56CO'; max_energy = 80; max_cs= 580
    Z = '27'; A = '56'; betaFeeding='beta+'; branchingratio = 1.0; 
    betaPlusDecayChain = {isotope: ['28', branchingratio, None]}
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ni.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot_manual(element, isotope, independent=False)
    CrossSection().plot_manual(element, isotope='56NI', independent=False, color='hotpink')
    exfor.plotExforDataFromFilename(filename='Ni_56Co_cumulative.txt', setLegend=False, maxE = max_energy)
    # CrossSection().plot(element, isotope, independent=False)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Ni_56Co_cumulative.pdf')
    plt.show()
    
def Ni_56Co_independent():
    element = 'Ni'; isotope = '56CO'; max_energy = 80; max_cs= 580
    Z = '27'; A = '56'; betaFeeding=None; branchingratio = None; 
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ni.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = None,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot_manual(element, isotope, independent=True)
    CrossSection().plot_manual(element, isotope, independent=False, color='forestgreen')
    CrossSection().plot_manual(element, isotope='56NI', independent=True, color='hotpink')
    # CrossSection().plot(element, isotope, independent=False)


    total, unc_total = cross_sections(element, isotope, independent=False)
    feeding, unc_feeding = cross_sections(element, isotope='56NI', independent=False)
    xs, unc_xs = CrossSection().subtract_beta(total, unc_total, feeding, unc_feeding, branching_ratio=1.0)
    energy, unc_left, unc_right = weighted_average_beam_energy(element)
    # CrossSection().plot_manual(element, isotope, independent=True, color='hotpink')
    # CrossSection().plot_manual(element, isotope, independent=False, color='dodgerblue')
    CrossSection().plot_feeding_corrected(xs, unc_xs, energy, unc_left, unc_right)

    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=True)
    plt.savefig(save_cs_figures_to + 'Ni_56Co_independent.pdf')
    plt.show()

def Ni_55Co_cumulative():
    element = 'Ni'; isotope = '55CO'; max_energy = 80; max_cs= 55
    Z = '27'; A = '55'; betaFeeding='beta+'; branchingratio = 1.0; 
    betaPlusDecayChain = {isotope: ['28', branchingratio, None]}
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ni.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot_manual(element, isotope, independent=False)
    # CrossSection().plot_manual(element, isotope='56NI', independent=False, color='hotpink')
    exfor.plotExforDataFromFilename(filename='Ni_55Co_cumulative.txt', setLegend=False, maxE = max_energy)
    # CrossSection().plot(element, isotope, independent=False)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.savefig(save_cs_figures_to + 'Ni_55Co_cumulative.pdf')
    plt.show()

def Ni_58Co_cumulative():
    element = 'Ni'; isotope='58CO'; max_energy=80; max_cs=100
    Z = '28'; A = '56'; betaFeeding='beta+'; branchingratio = 0.099988; 
    decaychain = {isotope: ['29', branchingratio]}
    # {isotope: [branchingRatio isomerLevel]} #isomer
    independent=False
    tendl(element, Z, A, isomerLevel=None, betaPlusDecayChain = None, betaMinusDecayChain = None, isomer_decay_chain = decaychain)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = None, branchingRatio =branchingratio, parentIsomerLevel = '01')
    exfor.plotExforDataFromFilename(filename='Ni_58Co_cumulative.txt', setLegend=False, maxE = max_energy)
    CrossSection().plot_manual(element, isotope, independent=independent)
    # CrossSection().plot_manual(element, isotope, independent=False, color='brown')
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ni_58Co_cumulative.pdf')
    plt.show()

def Ni_54Mn_independent():
    element = 'Ni'; isotope = '54MN'; max_energy = 80; max_cs= 14
    Z = '25'; A = '54'; betaFeeding=None; branchingratio = None; decaychain=None
    independent=True
    exfor.plotExforDataFromFilename(filename='Ni_54Mn_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl(element, Z, A)
    # CrossSection().plot(element, isotope, )
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ni_54Mn_independent.pdf')
    plt.show()


def Ni_52Mn_cumulative():
    element = 'Ni'; isotope='52MN'; max_energy=80; max_cs=31
    Z = '25'; A = '52'; betaFeeding='beta+'; branchingratio = 1.0; 
    decaychain = {isotope: ['26', branchingratio]}
    # {isotope: [branchingRatio isomerLevel]} #isomer
    independent=False
    tendl(element, Z, A, isomerLevel=None, betaPlusDecayChain = None, betaMinusDecayChain = None, isomer_decay_chain = decaychain)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    # talys.plotTalys(productZ='26', productA=A, targetFoil=element, isomerLevel = None, betaFeeding = None)
    exfor.plotExforDataFromFilename(filename='Ni_52Mn_cumulative.txt', setLegend=False, maxE = max_energy)
    CrossSection().plot_manual(element, isotope, independent=independent)
    CrossSection().plot(element, isotope)
    # CrossSection().plot_manual(element, isotope, independent=False, color='brown')
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ni_52Mn_cumulative.pdf')
    plt.show()

def Ni_52mMn_cumulative():
    #FALSE
    element = 'Ni'; isotope='52MNm'; max_energy=80; max_cs=31
    Z = '25'; A = '52'; betaFeeding='beta+'; branchingratio = 1.0; 
    decaychain = {isotope: ['26', branchingratio]}
    # {isotope: [branchingRatio isomerLevel]} #isomer
    independent=False
    tendl(element, Z, A, isomerLevel='01', betaPlusDecayChain = None, betaMinusDecayChain = None, isomer_decay_chain = decaychain)
    # talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = '01', betaFeeding = betaFeeding, branchingRatio =branchingratio)
    # talys.plotTalys(productZ='26', productA=A, targetFoil=element, isomerLevel = None, betaFeeding = None)
    # exfor.plotExforDataFromFilename(filename='Ni_52Mn_cumulative.txt', setLegend=False, maxE = max_energy)
    CrossSection().plot_manual(element, isotope, independent=independent)
    # CrossSection().plot_manual(element, isotope, independent=False, color='brown')
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    # plt.savefig(save_cs_figures_to + 'Ni_52Mn_cumulative.pdf')
    plt.show()

def Ni_51Cr_cumulative():
    element = 'Ni'; isotope='51CR'; max_energy=80; max_cs=31
    Z = '24'; A = '51'; betaFeeding='beta+'; branchingratio = 1.0; 
    decaychain = {isotope: ['25', branchingratio]}
    # {isotope: [branchingRatio isomerLevel]} #isomer
    independent=False
    tendl(element, Z, A, isomerLevel=None, betaPlusDecayChain = None, betaMinusDecayChain = None, isomer_decay_chain = decaychain)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    # talys.plotTalys(productZ='26', productA=A, targetFoil=element, isomerLevel = None, betaFeeding = None)
    exfor.plotExforDataFromFilename(filename='Ni_51Cr_cumulative.txt', setLegend=False, maxE = max_energy)
    CrossSection().plot_manual(element, isotope, independent=independent)
    # CrossSection().plot_manual(element, isotope, independent=False, color='brown')
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ni_51Cr_cumulative.pdf')
    plt.show()

# def Ni_58Co_independent():
#     element = 'Ni'; isotope = '58CO'; max_energy = 80; max_cs= 580
#     Z = '27'; A = '58'; betaFeeding=None; branchingratio = None; 
#     talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
#     CrossSection().plot_manual(element, isotope, independent=True)
#     CrossSection().plot_manual(element, isotope, independent=False, color='darkred')
#     CrossSection().plot_manual(element, isotope='58COm', independent=True, color='hotpink')
#     # CrossSection().plot(element, isotope, independent=False)
#     plt.xlim(0,max_energy)
#     plt.ylim(0,max_cs)
#     legend(fontsize='xx-small')
#     set_titles(element=element, isotope=isotope, independent=True)
#     # plt.savefig(save_cs_figures_to + 'Ni_56Co_independent.pdf')
#     # element = 'Ni'; isotope='58CO'; max_energy=80
#     # exfor.plotExforDataFromFilename(filename='Ni_58Co_independent.txt', setLegend=False, maxE = max_energy)
#     # CrossSection().plot_manual(element, isotope, independent=True)
#     # # CrossSection().plot_manual(element, isotope, independent=False, color='brown')
#     legend()
#     plt.show()

# def Ni_58mCo_independent():
#     br = 0.999988
#     element = 'Ni'; isotope='58CO'; max_energy=80
#     # exfor.plotExforDataFromFilename(filename='Cu_58mCo_independent.txt', setLegend=False, maxE = max_energy)
#     energy, unc_left, unc_right = weighted_average_beam_energy(element)
#     total, unc_total = cross_sections(element, isotope, independent=False)
#     gs, unc_gs = cross_sections(element, isotope, independent=True)
#     print(gs)
#     isomer = total-gs*br
#     print(isomer)
#     unc_isomer = unc_total-unc_gs*br
#     # print(yerr)
#     # print(unc_isomer)
#     plt.errorbar(energy, isomer, yerr=np.abs(unc_isomer), xerr=[unc_left, unc_right], marker='o', ls='none', label='isomer')
# #     plt.plot(energy, isomer, 'o')
#     plt.errorbar(energy, total, yerr=unc_total, xerr=[unc_left, unc_right], marker='o', ls='none', label='total')
#     plt.errorbar(energy,    gs, yerr=np.abs(unc_gs), xerr=[unc_left,    unc_right], marker='o', ls='none', label='groundstate')
#     legend()
#     plt.show()

def Ta_175Ta():
    element = 'Ta'; isotope = '175TA'
    # Exfor(exfor_path).plotExforDataFromFilename(filename='Ni_57Ni.txt', setLegend=False, maxE = 80)
    CrossSection().plot(element, isotope)
    plt.xlim(0,80)
    # plt.ylim(0,300)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope)
    plt.show()

def Ta_177Ta_cumulative():
    element = 'Ta'; isotope = '177TA'; max_energy = 80; max_cs=850
    Z = '73'; A = '177'; betaFeeding='beta+'; branchingratio = 1.0;
    betaPlusDecayChain = {isotope: ['74', branchingratio, None]}
    exfor.plotExforDataFromFilename(filename='Ta_177Ta_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    # CrossSection().save_manual(element, isotope, independent=False)
    CrossSection().plot_manual(element, isotope, independent=False)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=False)
    plt.savefig(save_cs_figures_to + 'Ta_177Ta_cumulative.pdf')
    plt.show()

def Ta_177W_independent():
    element = 'Ta'; isotope = '177W'; max_energy = 80; max_cs=1000
    Z = '74'; A = '177'; betaFeeding=None; branchingratio = None; 
    betaPlusDecayChain = None
    independent=True
    exfor.plotExforDataFromFilename(filename='Ta_177W.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    CrossSection().save_manual(element, isotope, independent=independent)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ta_177W_independent.pdf')
    plt.show()

def Ta_179W_cumulative():
    element = 'Ta'; isotope = '179W'; max_energy = 80; max_cs=1267
    Z = '74'; A = '179'; betaFeeding=None; branchingratio = None; 
    betaPlusDecayChain = None
    independent=False
    # exfor.plotExforDataFromFilename(filename='Ta_177W.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    # CrossSection().save_manual(element, isotope, independent=independent)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ta_179W_cumulative.pdf')
    plt.show()

def Ta_179mW_independent():
    element = 'Ta'; isotope = '179Wm'; max_energy = 80; max_cs=1267
    Z = '74'; A = '179'; betaFeeding=None; branchingratio = None; 
    betaPlusDecayChain = None; isomerLevel='02'
    independent=True
    # exfor.plotExforDataFromFilename(filename='Ta_177W.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = isomerLevel, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = isomerLevel,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = '00',  
    #                               betaPlusDecayChain = betaPlusDecayChain,
    #                               color='orange', lineStyle=None, label=None, semilog_y=False)
    # tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
    #                               betaPlusDecayChain = betaPlusDecayChain,
    #                               color='yellow', lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    # CrossSection().save_manual(element, isotope, independent=independent)
    # CrossSection().plot_manual(element, isotope='179W', independent=False, color='red')
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ta_179mW_independent.pdf')
    plt.show()

def Ta_179Ta_cumulative():
    element = 'Ta'; isotope = '179TA'; max_energy = 80; max_cs=1200
    Z = '73'; A = '179'; betaFeeding=None; branchingratio = None; 
    betaPlusDecayChain = None; independent=False
    # exfor.plotExforDataFromFilename(filename='Ta_177W.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    CrossSection().save_manual(element, isotope, independent)
    CrossSection().plot_manual(element, isotope, independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=False)
    plt.savefig(save_cs_figures_to + 'Ta_179Ta_cumulative.pdf')
    plt.show()

    # ADD SUBTRACTION 

def Ta_178W_independent():
    element = 'Ta'; isotope = '178W'; max_energy = 80; max_cs=1200
    Z = '74'; A = '178'; betaFeeding=None; branchingratio = None; 
    betaPlusDecayChain = None
    exfor.plotExforDataFromFilename(filename='Ta_178W_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope, independent=True)
    # CrossSection().save_manual(element, isotope, independent=True)
    CrossSection().plot_manual(element, isotope, independent=True)
    # CrossSection().save_manual(element, isotope, independent=None)
    CrossSection().plot_manual(element, isotope, independent=None, color='darkred')
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=True)
    plt.savefig(save_cs_figures_to + 'Ta_178W_independent.pdf')
    plt.show()

def Ta_178Ta_cumulative():
    element = 'Ta'; isotope = '178TA'; max_energy = 80; max_cs=250
    Z = '73'; A = '178'; betaFeeding=None; branchingratio = None; 
    betaPlusDecayChain = None; independent=False
    # exfor.plotExforDataFromFilename(filename='Ta_178W_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope, independent=True)
    # CrossSection().save_manual(element, isotope, independent=independent)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ta_178Ta_cumulative.pdf')
    plt.show()

def Ta_176Ta_independent():
    element = 'Ta'; isotope = '176TA'; max_energy = 80; max_cs=70
    Z = '73'; A = '176'; betaFeeding=None; branchingratio = None; 
    betaPlusDecayChain = None
    # exfor.plotExforDataFromFilename(filename='Ta_177W.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope, independent=True)
    CrossSection().plot_manual(element, isotope, independent=True)
    CrossSection().plot_manual(element, isotope='176W', independent=True, color='blue')
    CrossSection().plot_manual(element, isotope, independent=False, color='red')
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=True)
    plt.savefig(save_cs_figures_to + 'Ta_176Ta_independent.pdf')
    plt.show()

def Ta_176Ta_cumulative():
    element = 'Ta'; isotope = '176TA'; max_energy = 80; max_cs=670; independent=False
    Z = '73'; A = '176'; betaFeeding='beta+'; branchingratio = 1.0; 
    betaPlusDecayChain = {isotope: ['74', branchingratio, None]}
    exfor.plotExforDataFromFilename(filename='Ta_176Ta_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    # CrossSection().save_manual(element, isotope, independent=False)
    CrossSection().plot_manual(element, isotope, independent=False)
    plt.xlim(0,max_energy)
    plt.ylim(-5,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ta_176Ta_cumulative.pdf')
    plt.show()

def Ta_176Ta_independent():
    element = 'Ta'; isotope = '176TA'; max_energy = 80; max_cs=660; independent=False
    Z = '73'; A = '176'; betaFeeding=None; branchingratio = None; 
    betaPlusDecayChain = None
    # exfor.plotExforDataFromFilename(filename='Ta_176Ta_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    # CrossSection().save_manual(element, isotope, independent=True)
    CrossSection().plot_manual(element, isotope, independent=True, color='yellow', label='176Ta independent')
    CrossSection().plot_manual(element, isotope, independent=False, color='darkred', label='176Ta cumulative')
    # CrossSection().save_manual(element, isotope='176W', independent=True)
    CrossSection().plot_manual(element, isotope='176W', color='hotpink', independent=True, label='176W independent')
    # CrossSection().plot_manual(element, isotope='176W', independent=True)
    plt.xlim(0,max_energy)
    plt.ylim(-5,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=True)
    plt.savefig(save_cs_figures_to + 'Ta_176Ta_independent.pdf')
    plt.show()

def Ta_176W_independent():
    element = 'Ta'; isotope = '176W'; max_energy = 80; max_cs=670; independent=False
    Z = '74'; A = '176'; betaFeeding=None; branchingratio = None; 
    betaPlusDecayChain = None
    # exfor.plotExforDataFromFilename(filename='Ta_176Ta_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    # CrossSection().save_manual(element, isotope, independent=True)
    # CrossSection().plot_manual(element, isotope, independent=True, color='yellow', label='176Ta independent')
    # CrossSection().plot_manual(element, isotope, independent=False, color='darkred', label='176Ta cumulative')
    # CrossSection().save_manual(element, isotope='176W', independent=True)
    CrossSection().plot_manual(element, isotope='176W', independent=True)# color='hotpink', independent=True, label='176W independent')
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=True)
    plt.savefig(save_cs_figures_to + 'Ta_176W_independent.pdf')
    plt.show()

def Ta_175Ta_cumulative():
    element = 'Ta'; isotope = '175TA'; max_energy = 80; max_cs=245; independent=False
    Z = '73'; A = '175'; betaFeeding='beta+'; branchingratio = 1.0; 
    betaPlusDecayChain = {isotope: ['74', branchingratio, None]}
    exfor.plotExforDataFromFilename(filename='Ta_175Ta_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    CrossSection().save_manual(element, isotope, independent=False)
    CrossSection().plot_manual(element, isotope, independent=False)
    plt.xlim(0,max_energy)
    plt.ylim(-1,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ta_175Ta_cumulative.pdf')
    plt.show()


def Ta_173Ta_cumulative():
    element = 'Ta'; isotope = '173TA'; max_energy = 80; max_cs=245; independent=False
    Z = '73'; A = '173'; betaFeeding='beta+'; branchingratio = 1.0; 
    betaPlusDecayChain = {isotope: ['74', branchingratio, None]}
    exfor.plotExforDataFromFilename(filename='Ta_175Ta_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    CrossSection().save_manual(element, isotope, independent=independent)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(-1,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ta_173Ta_cumulative.pdf')
    plt.show()

def Ta_180mHF_independent():
    element = 'Ta'; isotope = '180HFm'; max_energy = 80; max_cs=4.2; independent=True
    Z = '72'; A = '180'; betaFeeding=None; branchingratio = None;
    betaPlusDecayChain = None
    exfor.plotExforDataFromFilename(filename='Ta_180mHf_independent.txt', setLegend=False, maxE = max_energy)
    # talys_old.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = '07', betaFeeding = betaFeeding, branchingRatio =branchingratio)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = '06', betaFeeding = betaFeeding, branchingRatio =branchingratio)
    # talys_old.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = '06', betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = '06', color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    # CrossSection().save_manual(element, isotope, independent=independent)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')

    plt.title(f"$^{{nat}}${element}(p,x)$^{{180m}}$Hf")
    
    plt.savefig(save_cs_figures_to + 'Ta_180mHf_independent.pdf')
    plt.show()

def Ta_179Lu_independent():
    element = 'Ta'; isotope = '179LU'; max_energy = 80; max_cs=None; independent=True
    Z = '71'; A = '179'; betaFeeding=None; branchingratio = None;
    betaPlusDecayChain = None
    # exfor.plotExforDataFromFilename(filename='Ta_180mHf_independent.txt', setLegend=False, maxE = max_energy)
    # talys_old.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = '07', betaFeeding = betaFeeding, branchingRatio =branchingratio)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(-5,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    # plt.savefig(save_cs_figures_to + 'Ta_179Lu_independent.pdf')
    plt.show()

def Ta_172Hf_cumulative():
    element = 'Ta'; isotope = '172HF'; max_energy = 80; max_cs=11; independent=False
    Z = '72'; A = '172'; betaFeeding=None; branchingratio = None;
    betaPlusDecayChain = None
    exfor.plotExforDataFromFilename(filename='Ta_172Hf_cumulative.txt', setLegend=False, maxE = max_energy)
    # talys_old.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = '07', betaFeeding = betaFeeding, branchingRatio =branchingratio)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    CrossSection().save_manual(element, isotope, independent=independent)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ta_172Hf_cumulative.pdf')
    plt.show()

def Ta_173Hf_cumulative():
    element = 'Ta'; isotope = '173HF'; max_energy = 80; max_cs=21; independent=False
    Z = '72'; A = '173'; betaFeeding=None; branchingratio = None;
    betaPlusDecayChain = None
    exfor.plotExforDataFromFilename(filename='Ta_173Hf_cumulative.txt', setLegend=False, maxE = max_energy)
    # talys_old.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = '07', betaFeeding = betaFeeding, branchingRatio =branchingratio)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    CrossSection().save_manual(element, isotope, independent=independent)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ta_173Hf_cumulative.pdf')


    # tendl_Ta.get_q_value('73', '173', None)
    # tendl_Ta.get_q_value('74', '173', None)
    plt.show()


def Ta_175Hf_cumulative(): # FIX feeding
    element = 'Ta'; isotope = '175HF'; max_energy = 80; max_cs=21; independent=False
    Z = '72'; A = '175'; betaFeeding='beta+'; branchingratio = 1.0;
    betaPlusDecayChain = {isotope: ['73', branchingratio, None]}
    exfor.plotExforDataFromFilename(filename='Ta_175Hf_cumulative.txt', setLegend=False, maxE = max_energy)
    # talys_old.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, betaPlusDecayChain=betaPlusDecayChain, isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    # CrossSection().save_manual(element, isotope, independent=independent)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ta_175Hf_cumulative.pdf')
    plt.show()

def Ta_175Hf_independent():
    element = 'Ta'; isotope = '175HF'; max_energy = 80; max_cs=21; independent=True
    Z = '72'; A = '175'; betaFeeding=None; branchingratio = None;
    betaPlusDecayChain = None
    exfor.plotExforDataFromFilename(filename='Ta_175Hf_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    
    energy, unc_left, unc_right, subtract_xs, unc_subtract_xs = cross_sections(element, isotope='175TA', independent=False)
    energy, unc_left, unc_right, cumulative_xs, unc_cumulative_xs = cross_sections(element, isotope='175HF', independent=False)
    CrossSection().plot_manual(element, isotope='175TA', independent=False, color='peru', label='175TA c')
    CrossSection().plot_manual(element, isotope='175HF', independent=False, color='forestgreen', label='175HF c')
    CrossSection().save_subtracted_beta(element, isotope, independent, cumulative_xs, unc_cumulative_xs, subtract_xs, unc_subtract_xs, branching_ratio=1.0)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ta_175Hf_independent.pdf')
    plt.show()

def Ta_176Hf_cumulative():
    element = 'Ta'; isotope = '176HF'; max_energy = 80; max_cs=21; independent=False
    Z = '72'; A = '176'; betaFeeding=None; branchingratio = None;
    betaPlusDecayChain = None
    # exfor.plotExforDataFromFilename(filename='Ta_173Hf_cumulative.txt', setLegend=False, maxE = max_energy)
    # talys_old.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = '07', betaFeeding = betaFeeding, branchingRatio =branchingratio)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    # plt.savefig(save_cs_figures_to + 'Ta_179Lu_independent.pdf')
    plt.show()

def Ta_179Lu_independent():
    element = 'Ta'; isotope = '179LU'; max_energy = 80; max_cs=None; independent=True
    Z = '71'; A = '172'; betaFeeding=None; branchingratio = None;
    betaPlusDecayChain = None
    # exfor.plotExforDataFromFilename(filename='Ta_172Hf_cumulative.txt', setLegend=False, maxE = max_energy)
    # talys_old.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = '07', betaFeeding = betaFeeding, branchingRatio =branchingratio)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    CrossSection().save_manual(element, isotope, independent=independent)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ta_179Lu_independent.pdf')
    plt.show()


def Ta_Lu():
    tendl_Ta.plotTendl23Unique(productZ='71', productA='179', isomerLevel = None, color=None, lineStyle=None, label='179Lu', semilog_y=False)
    plt.xlim(0,80)
    plt.show()

def Sn_test():
    element = 'Sn'; isotope = '117SBg'; max_energy = 80; max_cs=None; independent=True
    Z = '51'; A = '117'; betaFeeding=None; branchingratio = None;
    betaPlusDecayChain = None
    # exfor.plotExforDataFromFilename(filename='Ta_180mHf_independent.txt', setLegend=False, maxE = max_energy)
    # talys_old.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = '07', betaFeeding = betaFeeding, branchingRatio =branchingratio)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Sn.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(-5,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    # plt.savefig(save_cs_figures_to + 'Ta_179Lu_independent.pdf')
    plt.show()

def Ta_177Ta_independent():
    element = 'Ta'; isotope = '177TA'; max_energy = 80; max_cs=1000
    Z = '73'; A = '177'
    # CrossSection().plot_subtract(element, isotope_parent='177W', isotope_daughter=isotope, branching_ratio=1.0)
    tendl(element, Z, A)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element) #, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    independent=True
    # total, unc_total = cross_sections(element, isotope, independent=False)
    # feeding, unc_feeding = cross_sections(element, isotope='177W', independent=True)
    energy, unc_left, unc_right = weighted_average_beam_energy(element)
    energy, unc_left, unc_right, subtract_xs, unc_subtract_xs = cross_sections(element, isotope='177W', independent=True)
    energy, unc_left, unc_right, cumulative_xs, unc_cumulative_xs = cross_sections(element, isotope='177TA', independent=False)
    CrossSection().plot_manual(element, isotope='177W', independent=True, color='peru', label='177W')
    CrossSection().plot_manual(element, isotope='177TA', independent=False, color='dodgerblue', label='177TA c')
    xs, unc_xs = CrossSection().subtract_beta(cumulative_xs, unc_cumulative_xs, subtract_xs, unc_subtract_xs, branching_ratio=1.0)
    CrossSection().plot_feeding_corrected(xs, unc_xs, energy, unc_left, unc_right)
    
    # CrossSection().plot_manual(element, isotope, independent=False)
    # CrossSection().save_manual(element, isotope, independent=True)
    # CrossSection().plot_manual(element, isotope, independent=True)

    # CrossSection().plot_manual(element, isotope, independent=False)
    # CrossSection().plot_manual(element, isotope='177W', independent=True, color='darkred')

    plt.xlim(0,max_energy)
    plt.ylim(-230,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=True)
    plt.savefig(save_cs_figures_to + 'Ta_177Ta_independent.pdf')
    plt.show()

def Ta_180Ta_independent():
    element = 'Ta'; isotope = '180TA'; max_energy = 80; max_cs=350
    # element = 'Ta'; isotope = '180TA'; max_energy = 80; max_cs=350
    Z = '73'; A = '180'; betaFeeding=None; branchingratio = None; 
    # betaPlusDecayChain = {isotope: ['74', branchingratio, None]}
    betaPlusDecayChain = None; independent=True
    exfor.plotExforDataFromFilename(filename='Ta_180Ta_cumulative.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Ta.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
    tendl_Ta.plotTendl23Unique(productZ='75', productA=A, isomerLevel = None, color='red', lineStyle=None, label=None, semilog_y=False)
    CrossSection().plot(element, isotope)
    # CrossSection().save_manual(element, isotope, independent=independent)
    CrossSection().plot_manual(element, isotope, independent=independent)
    # CrossSection().plot_manual(element, isotope, independent=True)
    # CrossSection().plot_manual(element, isotope='175TA', independent=True, color='hotpink')
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=True)
    plt.savefig(save_cs_figures_to + 'Ta_180Ta_independent.pdf')
    plt.show()

# def Ta_178Ta_independent():   # Must use x ray detector... 
#     element = 'Ta'; isotope = '180TA'; max_energy = 80; max_cs=350
#     Z = '73'; A = '180'; betaFeeding=None; branchingratio = None; 
#     # betaPlusDecayChain = {isotope: ['74', branchingratio, None]}
#     betaPlusDecayChain = None
#     exfor.plotExforDataFromFilename(filename='Ta_180Ta_cumulative.txt', setLegend=False, maxE = max_energy)
#     talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
#     tendl_Ta.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)
#     # CrossSection().plot(element, isotope, independent=True)
#     CrossSection().plot_manual(element, isotope, independent=True)
#     plt.xlim(0,max_energy)
#     plt.ylim(0,max_cs)
#     legend(fontsize='xx-small')
#     set_titles(element=element, isotope=isotope, independent=True)
#     plt.savefig(save_cs_figures_to + 'Ta_180Ta_independent.pdf')
#     plt.show()

# def Ta_178Ta_cumulative():
#     element = 'Ta'; isotope = '178TAm'; max_energy = 80; max_cs=850
#     # Z = '73'; A = '178'; betaFeeding='beta+'; branchingratio = 1.0;
#     Z = '73'; A = '178'; betaFeeding=None; branchingratio = None;
#     independent=False
#     betaPlusDecayChain = None#{isotope: ['74', branchingratio, None]}
#     exfor.plotExforDataFromFilename(filename='Ta_178Ta_cumulative.txt', setLegend=False, maxE = max_energy)
#     talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
#     tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
#                                   betaPlusDecayChain = betaPlusDecayChain,
#                                   color=None, lineStyle=None, label=None, semilog_y=False)
#     # CrossSection().plot(element, isotope)
#     CrossSection().save_manual(element, isotope, independent=independent)
#     CrossSection().plot_manual(element, isotope, independent=independent)
#     plt.xlim(0,max_energy)
#     plt.ylim(0,max_cs)
#     legend(fontsize='xx-small')
#     set_titles(element=element, isotope=isotope, independent=independent)
#     # plt.savefig(save_cs_figures_to + 'Ta_178Ta_cumulative.pdf')
#     plt.show()

# def Ta_178Ta_cumulative():
#     element = 'Ta'; isotope = '178TAm'; max_energy = 80; max_cs=850
#     # Z = '73'; A = '178'; betaFeeding='beta+'; branchingratio = 1.0;
#     Z = '73'; A = '178'; betaFeeding=None; branchingratio = None;
#     betaPlusDecayChain = None#{isotope: ['74', branchingratio, None]}
#     exfor.plotExforDataFromFilename(filename='Ta_178Ta_cumulative.txt', setLegend=False, maxE = max_energy)
#     talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
#     tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
#                                   betaPlusDecayChain = betaPlusDecayChain,
#                                   color=None, lineStyle=None, label=None, semilog_y=False)
#     # CrossSection().plot(element, isotope)
#     CrossSection().save_manual(element, isotope, independent=True)
#     CrossSection().plot_manual(element, isotope, independent=True)
#     # CrossSection().plot_manual(element, isotope, independent=False)
#     plt.xlim(0,max_energy)
#     plt.ylim(0,max_cs)
#     legend(fontsize='xx-small')
#     set_titles(element=element, isotope=isotope, independent=False)
#     plt.savefig(save_cs_figures_to + 'Ta_178Ta_cumulative.pdf')
#     plt.show()

# def Ta_179W_cumulative():
#     element = 'Ta'; isotope = '179W'; max_energy = 80; max_cs= 200
#     Z = '74'; A = '179'; betaFeeding=None; branchingratio = None; decaychain=None
#     independent=False
#     # exfor.plotExforDataFromFilename(filename='Ni_61Cu_independent.txt', setLegend=False, maxE = max_energy)
#     talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
#     tendl(element, Z, A)
#     # CrossSection().plot(element, isotope)
#     CrossSection().plot_manual(element, isotope, independent=independent)
#     plt.xlim(0,max_energy)
#     plt.ylim(0,max_cs)
#     legend(fontsize='xx-small')
#     set_titles(element=element, isotope=isotope, independent=independent)
#     # plt.savefig(save_cs_figures_to + 'Ni_61Cu_independent.pdf')
#     plt.show()


def Ta_181W_independent():
    element = 'Ta'; isotope = '181W'; max_energy = 80; max_cs= 100
    Z = '74'; A = '181'; betaFeeding=None; branchingratio = None; decaychain=None
    independent=True
    # exfor.plotExforDataFromFilename(filename='Ni_61Cu_independent.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = None, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl(element, Z, A)
    # tendl(element, Z='7', A='172')
    # CrossSection().plot(element, isotope)
    

    # CrossSection().save_manual(element, isotope, independent=True)
    CrossSection().plot_manual(element, isotope, independent=True)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    plt.savefig(save_cs_figures_to + 'Ta_181W_independent.pdf')
    plt.show()


def Sn_119Sb_cumulative():
    element = 'Sn'; isotope = '119SBg'; max_energy = 80; max_cs=1267
    Z = '51'; A = '119'; betaFeeding=None; branchingratio = None; 
    betaPlusDecayChain = None; isomerLevel=None
    independent=False
    # exfor.plotExforDataFromFilename(filename='Ta_177W.txt', setLegend=False, maxE = max_energy)
    talys.plotTalys(productZ=Z, productA=A, targetFoil=element, isomerLevel = isomerLevel, betaFeeding = betaFeeding, branchingRatio =branchingratio)
    tendl_Sn.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = isomerLevel,  
                                  betaPlusDecayChain = betaPlusDecayChain,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    # tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = '00',  
    #                               betaPlusDecayChain = betaPlusDecayChain,
    #                               color='orange', lineStyle=None, label=None, semilog_y=False)
    # tendl_Ta.plotTendl23Unique_feeding(productZ=Z, productA=A, isomerLevel = None,  
    #                               betaPlusDecayChain = betaPlusDecayChain,
    #                               color='yellow', lineStyle=None, label=None, semilog_y=False)
    # CrossSection().plot(element, isotope)
    CrossSection().save_manual(element, isotope, independent=independent)
    # CrossSection().plot_manual(element, isotope='179W', independent=False, color='red')
    CrossSection().plot_manual(element, isotope, independent=independent)
    plt.xlim(0,max_energy)
    plt.ylim(0,max_cs)
    legend(fontsize='xx-small')
    set_titles(element=element, isotope=isotope, independent=independent)
    # plt.savefig(save_cs_figures_to + 'Ta_179mW_independent.pdf')
    plt.show()

def legend(fontsize='small', loc='best'):
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(),fontsize=fontsize, loc=loc)

def set_titles(title=None, element=None, isotope=None, independent=False):
    plt.xlabel('Proton energy (MeV)')
    plt.ylabel('Cross section (mb)')
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
                                  isomerDecayChain = None,
                                  color=None, lineStyle=None, label=None, semilog_y=False)
    else:
        # print(tendl)
        tendl.plotTendl23Unique(productZ=Z, productA=A, isomerLevel = isomerLevel, color=None, lineStyle=None, label=None, semilog_y=False)


def tendl_reactions(ylim=None):
    tendl = tendl_Ta
    colors = [ 'crimson', 'dodgerblue', 'forestgreen', 'palevioletred', 'sienna', 'indianred','darkgoldenrod', 'firebrick', 'orchid',  'darkorange','dodgerblue', 'lime','mediumpurple', 'hotpink', 'pink', 'darkred', 'olive', 'cyan', 'indigo', 'springgreen', 'plum', 'darkslategrey', 'teal', 'rosybrown', 'yellow', 'burlywood']
    tendl.plotTendl23Unique(productZ='74', productA='181', isomerLevel = None, color=colors[0], lineStyle=None, label='181W (10 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='74', productA='179', isomerLevel = None, color=colors[2], lineStyle=None, label='179W (16 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='74', productA='178', isomerLevel = None, color=colors[3], lineStyle=None, label='178W (23 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='74', productA='177', isomerLevel = None, color=colors[4], lineStyle=None, label='177W (32 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='74', productA='176', isomerLevel = None, color=colors[5], lineStyle=None, label='176W (39 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='74', productA='175', isomerLevel = None, color=colors[6], lineStyle=None, label='175W ()', semilog_y=False)
    tendl.plotTendl23Unique(productZ='73', productA='180', isomerLevel = None, color=colors[7], lineStyle=None, label='180Ta (5 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='73', productA='179', isomerLevel = None, color=colors[8], lineStyle=None, label='179Ta (5 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='73', productA='178', isomerLevel = None, color=colors[9], lineStyle=None, label='178Ta (13 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='73', productA='177', isomerLevel = None, color=colors[10], lineStyle=None, label='177Ta (20 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='73', productA='176', isomerLevel = None, color=colors[11], lineStyle=None, label='176Ta (29 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='73', productA='175', isomerLevel = None, color=colors[12], lineStyle=None, label='175Ta (36 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='73', productA='174', isomerLevel = None, color=colors[13], lineStyle=None, label='174Ta (45 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='72', productA='180', isomerLevel = '06', color=colors[14], lineStyle=None, label='180mHf (6 MeV)',  semilog_y=False)
    tendl.plotTendl23Unique(productZ='72', productA='179', isomerLevel = '46', color=colors[15], lineStyle=None, label='179m2Hf (11 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='72', productA='178', isomerLevel = '99', color=colors[16], lineStyle=None, label='178m2Hf (0 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='72', productA='177', isomerLevel = '99', color=colors[17], lineStyle=None, label='177m2Hf (0 MeV)', semilog_y=False)
    tendl.plotTendl23Unique(productZ='72', productA='175', isomerLevel = None, color=colors[18], lineStyle=None, label='175Hf (13 MeV)',   semilog_y=False)
    tendl.plotTendl23Unique(productZ='72', productA='173', isomerLevel = None, color=colors[19], lineStyle=None, label='173Hf (28 MeV)',   semilog_y=False)
    tendl.plotTendl23Unique(productZ='72', productA='173', isomerLevel = None, color=colors[20], lineStyle=None, label='172Hf (47 MeV)',   semilog_y=False)


    # tendl.plotTendl23Unique(productZ='71', productA='179', isomerLevel = None, color=colors[21], lineStyle=None, label='179Lu (14 MeV)',   semilog_y=False)
    # tendl.plotTendl23Unique(productZ='71', productA='178', isomerLevel = None, color=colors[22], lineStyle=None, label='178Lu (13 MeV)',   semilog_y=False)
    # tendl.plotTendl23Unique(productZ='71', productA='177', isomerLevel = None, color=colors[23], lineStyle=None, label='177Lu (0 MeV)',   semilog_y=False)
    # tendl.plotTendl23Unique(productZ='71', productA='176', isomerLevel = '01', color=colors[24], lineStyle=None, label='176mLu (3 MeV)',   semilog_y=False)
    tendl.plotTendl23Unique(productZ='71', productA='175', isomerLevel = None, color=colors[25], lineStyle=None, label='175Lu (3 MeV)',   semilog_y=False)
    tendl.plotTendl23Unique(productZ='71', productA='174', isomerLevel = None, color=colors[25], lineStyle=None, label='174Lu (11 MeV)',   semilog_y=False)
    tendl.plotTendl23Unique(productZ='71', productA='173', isomerLevel = None, color=colors[25], lineStyle=None, label='173Lu (17 MeV)',   semilog_y=False)
    tendl.plotTendl23Unique(productZ='71', productA='172', isomerLevel = None, color=colors[25], lineStyle=None, label='172Lu (37 MeV)',   semilog_y=False)
    tendl.plotTendl23Unique(productZ='71', productA='171', isomerLevel = None, color=colors[25], lineStyle=None, label='171Lu (33 MeV)',   semilog_y=False)
    tendl.plotTendl23Unique(productZ='71', productA='170', isomerLevel = None, color=colors[25], lineStyle=None, label='170Lu (53 MeV)',   semilog_y=False)
    energies = weighted_average_beam_energy('Ta')[0]
    for e in energies:
        plt.axvline(e, linewidth='0.3')
    plt.xlim(0,60)
    if ylim:
        plt.ylim(0,ylim)
    plt.legend(fontsize='small')
    plt.show()
