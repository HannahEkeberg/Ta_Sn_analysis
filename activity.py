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
pathToPeakFiles_isotope = os.getcwd() + '/generatedfiles/peakdata/data_isotope/'
pathToActivityFiles = os.getcwd() + '/generatedfiles/activity/data/'
pathToActivityFiles_isotope = os.getcwd() + '/generatedfiles/activity/data_isotope/'

class Acitivity:
          
    def __init__(self):
        # self.eob_stack30 = '09/23/2025 18:35:00'
        self.eob_stack30 = '09/23/2025 18:40:00'
        # self.eob_stack55 = '09/24/2025 14:43:00'
        self.eob_stack55 = '09/24/2025 15:45:00'
        self.R_stack30 = [1.0, 1]
        self.R_stack55 = [1.0, 1]

    def getA0_single_isotope(self, isotope, foil, listOfPeakDataSummaries=None, guess=3.7e5, units = 'h', fitByA=True, plot=True, overwriteData=False, saveDecayChain=False):
        if listOfPeakDataSummaries==None:
            listOfPeakDataSummaries = self.listOfPeakSummaries(foil)
        peak_data = self.concat_peakData(listOfPeakDataSummaries)
        eob = self.getEob(foil)
        filter_peak_data = peak_data[peak_data['isotope'].astype(str).str.contains(isotope, case=False, na=False)]
        # print(filter_peak_data)
        self.decayrate_to_activity(filter_peak_data, eob)

        if fitByA:
            # print(peak_data)
            isotopes, fit, cov = self.fitByA(eob, peak_data, isotope, guess, units, plot)
        else:
            isotopes, fit, cov = self.fitByR(eob, peak_data, isotope, None, units, plot)

        df = self.eob_activity_dataframe(foil, isotopes, fit, cov, saveDecayChain)
        if overwriteData:
            df.to_csv(pathToActivityFiles + foil + '_' + isotope + '_' + '.csv')
        return df

    def getA0(self, element, isotope, plot_=False, compartment=None):
        foils = self.foils(element)
        eob_activity_list = []
        # generatedfiles/peakdata/data_isotope/Cu01_63ZNg_gammas.csv
        # generatedfiles/peakdata/data_isotope/Cu05_63ZN_gammas.csv
        for foil in foils:
            str = foil + '_' + isotope + '_gammas.csv'
            peak_data = pathToPeakFiles_isotope + str
            eob = self.getEob(foil)
            if not os.path.exists(peak_data):
                # if run == 0:
                    # isotope=isotope + 'g'
                # run +=1
                str = foil + '_' + isotope + 'g_gammas.csv'
                peak_data = pathToPeakFiles_isotope + str
            try:
                df = pd.read_csv(peak_data, comment="#")
                if compartment and compartment in foil and not plot_:
                    e = df['energy'].values; counts = df['counts'].values; unc_counts = df['unc_counts'].values; file = df['filename'].values
                    for i in range(len(e)):
                        print(e[i], unc_counts[i]/counts[i]*100, file[i])
                    plot=True
                else:
                    plot=plot_
                isotopes, fit, cov = self.fitByA(eob, df, isotope, plot=plot)
                eob_activity = fit[0]; unc_eob_activity = cov[0][0]
                if not np.isfinite(eob_activity) or not np.isfinite(unc_eob_activity):
                    eob_activity = 0; unc_eob_activity = 0
                data = [foil, isotope, eob_activity, unc_eob_activity]
                # df_activity = self.eob_activity_dataframe(foil, isotopes, fit, cov, saveDecayChain=saveDecayChain)
                eob_activity_list.append(data)
            except:
                print('no peak data found in foil ' + foil + ' for isotope ' + isotope)
                data = [foil, isotope, 0, 0]
                eob_activity_list.append(data)
        df_eob = pd.DataFrame(eob_activity_list, columns=['foil', 'isotope', 'fit', 'cov'])
        # if (df_eob["fit"] != 0.0).any():
        df_eob.to_csv(pathToActivityFiles_isotope +  element + '_' + isotope + '.csv') 

    def foils(self, element):
        n=14
        return [f"{element}{i:02d}" for i in range(1, n+1)]

    def getA0_mulitple_isotopes(self, listOfIsotopes, foil, listOfPeakDataSummaries=None, guess=3.7e5, units = 'h', fitByA=True, plot=True, overwriteData=True, saveDecayChain=False):
        if listOfPeakDataSummaries==None:
            listOfPeakDataSummaries = self.listOfPeakSummaries(foil)
        if listOfPeakDataSummaries:
            peak_data = self.concat_peakData(listOfPeakDataSummaries)
        else:
            peak_data = pd.DataFrame(columns=['filename','isotope','energy','counts','unc_counts','intensity','unc_intensity','efficiency','unc_efficiency','decays','unc_decays','decay_rate','unc_decay_rate','chi2','start_time','live_time','real_time'])
        if not listOfIsotopes:
            # listOfIsotopes = getListOfIsotopesPerFoil(foil)
            listOfIsotopes = peak_data['isotope'].unique().tolist()
        eob = self.getEob(foil)
        dfs = []
        for idx, isotope in enumerate(listOfIsotopes):
            try:
                if fitByA:
                    isotopes, fit, cov = self.fitByA(eob, peak_data, isotope, guess, units, plot)
                else:
                    isotopes, fit, cov = self.fitByR(eob, peak_data, isotope, None, units, plot)
            except:
                print("Did not get any activitites for isotope: " + isotope)
                isotopes = None
            if isotopes:
                df = self.eob_activity_dataframe(foil, isotopes, fit, cov, saveDecayChain)
                dfs.append(df)

        if dfs:
            df_concat = pd.concat(dfs, axis=0)
            if overwriteData:
                df_concat.to_csv(pathToActivityFiles + foil + '_all_isotopes.csv')
        # return df

    def fitByR(self, eob, peak_data, isotope, R_guess=[[1e4, 1]], units='h', plot=False):
        dc = ci.DecayChain(parent_isotope=isotope, R=R_guess, units=units)
        dc.get_counts(spectra='', EoB=eob, peak_data=peak_data)
        isotopes, fit, cov = dc.fit_R()
        A0 = dc.activity(isotope, time=0)
        if plot:
            dc.plot()
        return isotopes, fit, cov
    
    def fitByA(self, eob, peak_data, isotope, production_rate=1e6, units='h', plot=False):
        dc = ci.DecayChain(parent_isotope=isotope, A0=production_rate, units=units)
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
                data.append([foil, isotopes[i], fit[i], cov[i][i]])
        else:
            data.append([foil, isotopes[0], fit[0], cov[0][0]])
        df = pd.DataFrame(data, columns=['foil', 'isotope', 'fit', 'cov'] )
        return df
        # if overwriteData:
            # df.to_csv(pathToActivityFiles + foil + '_' + isotope + '_' + '.csv')

    def extractActivityManually(self, foil):
        peakDataSummaries = self.listOfPeakSummaries(foil)
        peak_data = self.concat_peakData(peakDataSummaries)
        # print(peak_data)
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
        
        new_df = pd.DataFrame(data, columns = ['isotope', 'foil', 'E gamma (keV)', 'A (Bq)', 'dA (Bq)', 'Nc', 'I gamma', 'efficiency', 'delay time (s)', 'count time (s)'])
        # print(new_df[new_df['isotope'] == '63ZNg'])
        return new_df
    
    def plotActivityManually(self, isotope, foil, data=None):
        self.decay_const = ci.Isotope(isotope).decay_const()
        if data == None:
            data = self.extractActivityManually(foil)
        data_isotope = data[data['isotope'].str.contains(isotope)]
        # print(data_isotope)
        # A = data_isotope['A (Bq)'].values; dA = data_isotope['dA (Bq)'].values; delay_time = data_isotope['delay time (s)'].values
        # popt, pcov = curve_fit(self.singleDecayCurve, delay_time, A, p0=1e6, sigma=dA, absolute_sigma=True)
        # time = np.max(delay_time)/3600 # hours
        # xplot = np.linspace(0,time,1000)
        # A0_estimated = self.singleDecayCurve(0, popt)
        # sigma_activity_estimated = np.sqrt(np.diagonal(pcov))   #Uncertainty in the fitting parameters# print(A_est)
        # plt.plot(xplot,self.singleDecayCurve(xplot*3600,*popt), color='tan', linewidth=0.9, label='fit')
        # plt.errorbar(delay_time/3600, A, color='darkolivegreen', linewidth=0.001,yerr=dA, elinewidth=0.5, ecolor='k', capthick=0.5,marker='*', label='activity')   # cap thickness for error bar color='blue')
        # plt.errorbar(0, A0_estimated, color='darkblue', linewidth=0.001,yerr=sigma_activity_estimated, elinewidth=0.5, ecolor='k', capthick=0.5,marker='+', label='eob activity: %.2f MBq' %(A0_estimated*1e-6 ))   # cap thickness for error bar color='blue')
        # plt.xlabel('Time since eob (h)')
        # plt.ylabel('Activity (Bq)')
        # plt.legend()
        # plt.title(foil + ' - ' + isotope)
        # plt.savefig('test.png')
        plt.show()

    def activity(self, Nc, dNc, eps, deps, I_gamma, dI_gamma, t_count, t_delay, isotope):
        decay_const = ci.Isotope(isotope).decay_const()
        activity = (Nc*decay_const) / (eps * I_gamma * (1-np.exp(-decay_const*t_count)) * np.exp(-decay_const * t_delay)    )
        sigma_activity = activity * np.sqrt((dNc/Nc)**2  + (deps/eps)**2 + (dI_gamma/I_gamma)**2 )
        return activity, sigma_activity
    
    def singleDecayCurve(self, t, A0_guess):
        A_est = A0_guess * np.exp(-self.decay_const*t)
        return A_est

    def decayrate_to_activity(self, peak_data, eob):
        # print(peak_data)
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
        # print(df)

    def concat_peakData(self, peakDataSummaries):
        dataframes = []
        for i in peakDataSummaries:
            dataframes.append(pd.read_csv(pathToPeakFiles + i))
        return pd.concat(dataframes, axis=0)
        # return df_concat
    
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

    def eob_activity_from_files(self, foils, isotope):
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
                            if eob_activity[i] == 0:
                                cov_eob_activity[i]  = 0
                            else:
                                cov_eob_activity[i] = df_isotope['cov'].values[0]
        return eob_activity, cov_eob_activity

# listOfIsotopes = getListOfIsotopesPerFoil('Cu01')


# Acitivity().getA0_single_isotope(isotope='65ZN', foil='Cu12', plot=True)
# Acitivity().getA0_single_isotope(isotope='65ZN', foil='Cu13', plot=True)
# Acitivity().getA0_single_isotope(isotope='65ZN', foil='Cu14', plot=True)
# Acitivity().plotActivityManually('57NI', 'Ni01')
# Acitivity().getA0('Cu', '63ZN')



# Acitivity().getA0_mulitple_isotopes(None, 'Ni01', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ni02', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ni03', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ni04', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ni05', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ni06', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ni07', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ni08', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ni09', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ni10', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ni11', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ni12', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ni13', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ni14', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ta01', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ta02', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ta03', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ta04', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ta05', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ta06', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ta07', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ta08', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ta09', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ta10', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ta11', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ta12', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ta13', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Ta14', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Cu01', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Cu02', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Cu03', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Cu04', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Cu05', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Cu06', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Cu07', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Cu08', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Cu09', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Cu10', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Cu11', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Cu12', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Cu13', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Cu14', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Sn01', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Sn02', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Sn03', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Sn04', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Sn05', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Sn06', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Sn07', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Sn08', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Sn09', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Sn10', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Sn11', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Sn12', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Sn13', units='h', plot=False, guess=1e6, fitByA=True)
# Acitivity().getA0_mulitple_isotopes(None, 'Sn14', units='h', plot=False, guess=1e6, fitByA=True)

# listOfPeakDataSummaries = ['HA10242025_Det2_Cu01_10cm_job_peak_data.csv','BR09242025_Cu01_52cm_IDM_peak_data.csv', 'BY09242025_Cu01_52cm_IDM_peak_data.csv', 'CM09242025_Cu01_40cm_IDM_peak_data.csv', 'DA09252025_Cu01_30cm_IDM_peak_data.csv', 'EL09262025_Cu01_10cm_IDM_peak_data.csv', 'FH09282025_Cu01_10cm_IDM_peak_data.csv']
# Acitivity().getA0('65ZN', listOfPeakDataSummaries)
# Acitivity().getA0('62ZN', listOfPeakDataSummaries)
# Acitivity().getA0('63ZN', listOfPeakDataSummaries)
# Acitivity().getA0('61CU', 'Cu01', plot=False)
# Acitivity().getA0('58COm', 'Cu01', units='d')
# Acitivity().getA0('58CO', 'Cu01', units='d', guess=1e2)



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
