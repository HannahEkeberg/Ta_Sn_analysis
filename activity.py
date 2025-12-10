import os 
import curie as ci
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
from scipy import interpolate
from scipy.optimize import curve_fit
from analyze_spectra import *


import sys
sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
from nuclearanalysistools.Tendl import *
from nuclearanalysistools.tools import *

pathToPeakFiles = os.getcwd() + '/generatedfiles/peakdata/data/'
pathToActivityFiles = os.getcwd() + '/generatedfiles/activity/data/'

class Acitivity:
          
    def __init__(self):
        self.eob_stack30 = '09/23/2025 18:35:00'
        self.eob_stack55 = '09/24/2025 14:43:00'

    def getA0_single_isotope(self, isotope, foil, listOfPeakDataSummaries=None, guess=3.7e5, units = 'h', fitByA=True, plot=True, overwriteData=True, saveDecayChain=False):
        if listOfPeakDataSummaries==None:
            listOfPeakDataSummaries = self.listOfPeakSummaries(foil)
        peak_data = self.concat_peakData(listOfPeakDataSummaries)
        eob = self.getEob(foil)
        filter_peak_data = peak_data[peak_data['isotope'].astype(str).str.contains(isotope, case=False, na=False)]

        self.decayrate_to_activity(filter_peak_data, eob)

        if fitByA:
            isotopes, fit, cov = self.fitByA(eob, peak_data, isotope, guess, units, plot)
        else:
            isotopes, fit, cov = self.fitByR(eob, peak_data, isotope, None, units, plot)

        df = self.eob_activity_dataframe(foil, isotopes, fit, cov, saveDecayChain)
        if overwriteData:
            df.to_csv(pathToActivityFiles + foil + '_' + isotope + '_' + '.csv')
        return df
    
    def getA0_mulitple_isotopes(self, listOfIsotopes, foil, listOfPeakDataSummaries=None, guess=3.7e5, units = 'h', fitByA=True, plot=True, overwriteData=True, saveDecayChain=False):
        if listOfPeakDataSummaries==None:
            listOfPeakDataSummaries = self.listOfPeakSummaries(foil)
        peak_data = self.concat_peakData(listOfPeakDataSummaries)
        if not listOfIsotopes:
            # listOfIsotopes = getListOfIsotopesPerFoil(foil)
            listOfIsotopes = peak_data['isotope'].unique().tolist()
        eob = self.getEob(foil)

        data = []
        dfs = []
        for idx, isotope in enumerate(listOfIsotopes):
            try:
                if fitByA:
                    isotopes, fit, cov = self.fitByA(eob, peak_data, isotope, guess, units, plot)
                else:
                    isotopes, fit, cov = self.fitByR(eob, peak_data, isotope, None, units, plot)
            except:
                print("Did not get any activitites for isotope: " + isotope)

            for i in range(len(isotopes)):
                data.append([foil, isotopes[i], fit[i], cov[i]])

            df = self.eob_activity_dataframe(foil, isotopes, fit, cov, saveDecayChain)
            dfs.append(df)

        df_concat = pd.concat(dfs, axis=0)
        if overwriteData:
            df_concat.to_csv(pathToActivityFiles + foil + '_all_isotopes.csv')
        return df

    def fitByR(self, eob, peak_data, isotope, R_guess=[[1e4, 1]], units='h', plot=False):
        dc = ci.DecayChain(parent_isotope=isotope, R=R_guess, units=units)
        dc.get_counts(spectra='', EoB=eob, peak_data=peak_data)
        isotopes, fit, cov = dc.fit_R()
        A0 = dc.activity(isotope, time=0)
        if plot:
            dc.plot()
        return isotopes, fit, cov
    
    def fitByA(self, eob, peak_data, isotope, A0_guess=1e4, units='h', plot=False):
        dc = ci.DecayChain(parent_isotope=isotope, A0=A0_guess, units=units)
        dc.get_counts(spectra='', EoB=eob, peak_data=peak_data)
        isotopes, fit, cov = dc.fit_A0()
        A0 = dc.activity(isotope, time=0)
        if plot:
            dc.plot()
        return isotopes, fit, cov

    def eob_activity_dataframe(self, foil, isotopes, fit, cov, saveDecayChain=False):
        data = []
        if len(fit)>1 and saveDecayChain:
            for i in range(len(isotopes)):
                data.append([foil, isotopes[i], fit[i], cov[i]])
        else:
            data.append([foil, isotopes[0], fit[0], cov[0][0]])
        df = pd.DataFrame(data, columns=['foil', 'isotope', 'fit', 'cov'] )
        return df
        # if overwriteData:
            # df.to_csv(pathToActivityFiles + foil + '_' + isotope + '_' + '.csv')

    def extractActivityManually(self, foil):
        peakDataSummaries = self.listOfPeakSummaries(foil)
        peak_data = self.concat_peakData(peakDataSummaries)
        data = []
        for index, row in peak_data.iterrows():
            isotope = row['isotope']
            fname = row['filename']
            foilName = foil
            E = row['energy']
            Nc = row['counts']
            d_Nc = row['unc_counts']
            I_gamma = row['intensity']
            d_I_gamma = row['unc_intensity']
            eps = row['efficiency']
            d_eps = row['unc_efficiency']
            specStartTime = row['start_time']
            dt1 = datetime.strptime(self.getEob(foil), '%m/%d/%Y %H:%M:%S')
            dt2 = datetime.strptime(specStartTime, '%m/%d/%Y %H:%M:%S')
            delayTime = (dt2 - dt1).total_seconds()
            countTime = row['real_time']
            A, dA = self.activity(Nc, d_Nc, eps, d_eps, I_gamma,d_I_gamma, countTime, delayTime, isotope)
            data.append([isotope, foilName, E, A, dA, Nc, I_gamma, eps, delayTime, countTime ])
        
        new_df = pd.DataFrame(data, columns = ['isotope', 'foil', 'E gamma', 'A', 'dA', 'Nc', 'I gamma', 'efficiency', 'delay time', 'count time'])
        return new_df
    
    def plotActivityManually(self, isotope, foil, data=None):
        self.lamb = ci.Isotope(isotope).decay_const()
        if data == None:
            data = self.extractActivityManually(foil)
        data_isotope = data[data['isotope'].str.contains(isotope)]
        print(data_isotope)
        A = data_isotope['A'].values; dA = data_isotope['dA'].values; delay_time = data_isotope['delay time'].values
        popt, pcov = curve_fit(self.singleDecayCurve, delay_time, A, p0=1e6, sigma=dA, absolute_sigma=True)
        # print(popt)
        time = np.max(delay_time)
        xplot = np.linspace(0,time,1000)
        
        A0_estimated = self.singleDecayCurve(popt, 0)
        print(A0_estimated)
        sigma_activity_estimated = np.sqrt(np.diagonal(pcov))   #Uncertainty in the fitting parameters# print(A_est)

        plt.plot(xplot,self.singleDecayCurve(xplot*3600,*popt),'r-', color='red')
        # print(A, delay_time/3600)
        plt.plot(delay_time, A, '.')
        plt.errorbar(delay_time, A, color='green', linewidth=0.001,yerr=dA, elinewidth=0.5, ecolor='k', capthick=0.5)   # cap thickness for error bar color='blue')
        # plt.show()

    def activity(self, Nc, dNc, eps, deps, I_gamma, dI_gamma, t_count, t_delay, isotope):
        lamb = np.log(2) / ci.Isotope(isotope).half_life()
        corr = (1 - np.exp(-lamb * t_count)) * np.exp(-lamb * t_delay)
        activity = Nc / (eps * I_gamma * corr)
        sigma_activity = activity * np.sqrt((dNc/Nc)**2  + (deps/eps)**2 + (dI_gamma/I_gamma)**2 )
        return activity, sigma_activity
    
    def singleDecayCurve(self, A0_guess, t):
        A_est = A0_guess * np.exp(-self.lamb*t)
        return A_est

    def decayrate_to_activity(self, peak_data, eob):
        print(peak_data)
        # Decay rate = counts / live_time * 1 / (efficiency_gamma * I_gamma)  
        decay_rate = peak_data['decay_rate'].values
        counts = peak_data['counts'].values
        unc_decay_rate = peak_data['unc_decay_rate'].values
        real_time = peak_data['real_time'].values
        live_time = peak_data['live_time'].values
        start_time = peak_data['start_time'].values
        efficiency = peak_data['efficiency'].values        
        I = peak_data['intensity'].values*100
        E = peak_data['energy'].values
        isotope = peak_data['isotope'].values
        lamb = ci.Isotope(isotope[0]).decay_const()
        A = decay_rate*lamb
        data = []
        for i in range(len(peak_data)):
            time_since_eob = Tools().date_diff(eob, start_time[i], units=None)
            data.append([isotope[i], time_since_eob, A[i], decay_rate[i], E[i], I[i]])
        df = pd.DataFrame(data, columns=['isotope','time since eob', 'activity', 'decay rate', 'E gamma', 'I gamma (%)'])
        print(df)

    def concat_peakData(self, peakDataSummaries):
        dataframes = []
        for i in peakDataSummaries:
            dataframes.append(pd.read_csv(pathToPeakFiles + i))
        df_concat = pd.concat(dataframes, axis=0)
        return df_concat
    
    def listOfPeakSummaries(self, foil):
        root = os.getcwd() + '/generatedfiles/peakdata/data/'
        filenames = []
        for filename in os.listdir(root):
            if foil in filename:
                filenames.append(filename)
        return filenames

    def getEob(self, foil):
        element = foil[0:2]
        stack1_numbs = ['01', '02', '03', '04', '05', '06', '07']
        stack2_numbs = ['08', '09','10', '11', '12', '13', '14']

        if foil in [element + number for number in stack1_numbs]:
            return self.eob_stack55
        elif foil in [element + number for number in stack2_numbs]:
            return self.eob_stack30
        else:
            raise ValueError("No valid foil: " + foil)
        

    def printActivitiesInFoils(self, element, isotope, units='h'):
        stack_numbs = ['01', '02', '03', '04', '05', '06', '07', '08', '09','10', '11', '12', '13', '14']
        foils = [element + number for number in stack_numbs]
        for f in foils:
            df =  self.getA0(isotope, f, units=units, plot=False, guess=1e6, fitByA=True)
            print(isotope + ' -->End of beam activity: ' + f)
            print(df)
            print("****")


# listOfIsotopes = getListOfIsotopesPerFoil('Cu01')

Acitivity().getA0_mulitple_isotopes(None, 'Cu01', units='h', plot=False, guess=1e6, fitByA=True)
Acitivity().getA0_mulitple_isotopes(None, 'Cu02', units='h', plot=False, guess=1e6, fitByA=True)
Acitivity().getA0_mulitple_isotopes(None, 'Cu03', units='h', plot=False, guess=1e6, fitByA=True)
Acitivity().getA0_mulitple_isotopes(None, 'Cu04', units='h', plot=False, guess=1e6, fitByA=True)
Acitivity().getA0_mulitple_isotopes(None, 'Cu05', units='h', plot=False, guess=1e6, fitByA=True)
Acitivity().getA0_mulitple_isotopes(None, 'Cu06', units='h', plot=False, guess=1e6, fitByA=True)
Acitivity().getA0_mulitple_isotopes(None, 'Cu07', units='h', plot=False, guess=1e6, fitByA=True)
Acitivity().getA0_mulitple_isotopes(None, 'Cu08', units='h', plot=False, guess=1e6, fitByA=True)
Acitivity().getA0_mulitple_isotopes(None, 'Cu09', units='h', plot=False, guess=1e6, fitByA=True)
Acitivity().getA0_mulitple_isotopes(None, 'Cu10', units='h', plot=False, guess=1e6, fitByA=True)
Acitivity().getA0_mulitple_isotopes(None, 'Cu11', units='h', plot=False, guess=1e6, fitByA=True)
Acitivity().getA0_mulitple_isotopes(None, 'Cu12', units='h', plot=False, guess=1e6, fitByA=True)
Acitivity().getA0_mulitple_isotopes(None, 'Cu13', units='h', plot=False, guess=1e6, fitByA=True)
Acitivity().getA0_mulitple_isotopes(None, 'Cu14', units='h', plot=False, guess=1e6, fitByA=True)

# listOfPeakDataSummaries = ['HA10242025_Det2_Cu01_10cm_job_peak_data.csv','BR09242025_Cu01_52cm_IDM_peak_data.csv', 'BY09242025_Cu01_52cm_IDM_peak_data.csv', 'CM09242025_Cu01_40cm_IDM_peak_data.csv', 'DA09252025_Cu01_30cm_IDM_peak_data.csv', 'EL09262025_Cu01_10cm_IDM_peak_data.csv', 'FH09282025_Cu01_10cm_IDM_peak_data.csv']
# Acitivity().getA0('65ZN', listOfPeakDataSummaries)
# Acitivity().getA0('62ZN', listOfPeakDataSummaries)
# Acitivity().getA0('63ZN', listOfPeakDataSummaries)
# Acitivity().getA0('61CU', 'Cu01', plot=False)
# Acitivity().getA0('58COm', 'Cu01', units='d')
# Acitivity().getA0('58CO', 'Cu01', units='d', guess=1e2)


# Acitivity().plotActivityManually('62ZN', 'Cu09')
# Acitivity().getA0('65ZN', 'Cu01', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0('65ZN', 'Cu02', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0('65ZN', 'Cu03', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0('65ZN', 'Cu05', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0('65ZN', 'Cu05', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0('65ZN', 'Cu06', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0('65ZN', 'Cu07', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0('65ZN', 'Cu08', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0('65ZN', 'Cu09', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0('65ZN', 'Cu10', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0('65ZN', 'Cu11', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0('65ZN', 'Cu12', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0('65ZN', 'Cu13', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0('65ZN', 'Cu14', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0('63ZN', 'Cu06', units='h', plot=True, guess=1e6, fitByR=False)
# Acitivity().getA0_single_isotope('62CU', 'Cu02', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_single_isotope('62ZN', 'Cu02', units='h', plot=False, guess=1e6, fitByA=True, saveDecayChain=False)




# Acitivity().printActivitiesInFoils('Cu', '63ZN')


# Tendl({"Cu63": 0.6915, "Cu65": 0.3085}, 'proton').plotTendl23Unique(productZ='29', productA='63', isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False)



# Tendl({"Cu63": 0.6915, "Cu65": 0.3085}, 'proton').plotTendl23Unique(productZ='27', productA='58', isomerLevel = '00', color=None, lineStyle=None, label=None, semilog_y=False)
# Tendl({"Cu63": 0.6915, "Cu65": 0.3085}, 'proton').plotTendl23Unique(productZ='27', productA='58', isomerLevel = '01', color='r', lineStyle=None, label=None, semilog_y=False)
# Tendl({"Cu63": 0.6915, "Cu65": 0.3085}, 'proton').plotTendl23Unique(productZ='27', productA='58', isomerLevel = None, color='green', lineStyle=None, label=None, semilog_y=False)
# plt.show()
