import os
import pandas as pd
import numpy as np
import curie as ci
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import traceback


pathToPeakFiles_isotope = os.getcwd() + '/generatedfiles/peakdata/data_isotope/'
pathToActivityFiles_isotope = os.getcwd() + '/generatedfiles/activity/data_isotope2/'
pathToFigures = os.getcwd() + '/generatedfiles/activity/figures/'
test_list = []
class Activity_manual:

    def __init__(self, saveIndependent=None, saveIndependentParent=None):
        self.eob_stack30 = '09/23/2025 18:40:05'
        self.eob_stack55 = '09/24/2025 15:45:32'
        self.saveIndependent = saveIndependent
        self.saveIndependentParent = saveIndependentParent


    def decay_using_daughter_gamma(self, element, isotope_parent, isotope_daughter, compartment, plot_curve, min_half_lives=None, filtering=True, max_half_lives=None, max_activity_uncertainty=None):
        iso = ci.Isotope(isotope_parent)
        decay_constant = iso.decay_const()
        half_life = iso.half_life()
        foils = self.foils(element)
        eob_activity_list = []
        for foil in foils:
            peak_data = self.retrieve_peak_data(foil, isotope_daughter)
            print(peak_data)
            try:
                df = pd.read_csv(peak_data, comment="#")
                # print(df)
                delay_time, activity, unc_activity = self.calculate_activities(df, foil, decay_constant)
                if filtering:
                    delay_time, activity, unc_activity = self.filter_activities_on_delay_time(delay_time, activity, unc_activity, half_life,
                                        min_half_lives=min_half_lives, max_half_lives=max_half_lives, max_activity_uncertainty=max_activity_uncertainty)
                if compartment and compartment in foil and plot_curve==False:
                    plot=True
                else:
                    plot=plot_curve
                A0, unc_A0 = self.fit_activity_first_order(delay_time, activity, unc_activity, decay_constant, plot=plot)
                data = [foil, isotope_parent, A0, unc_A0]
                eob_activity_list.append(data)
                if plot:
                    self.savefig(element, isotope_parent, foil)
            except:
                print('no peak data found in foil ' + foil + ' for isotope ' + isotope_parent)
                data = [foil, isotope_parent, 0, 0]
                eob_activity_list.append(data)
        df = pd.DataFrame(eob_activity_list, columns=['foil', 'isotope', 'fit', 'cov'])
        self.savecsv(df, element, isotope_parent)
        pass

    def one_step_decay(self, element, isotope, compartment, plot_curve, min_half_lives=None, filtering=True, max_half_lives=None, max_activity_uncertainty=None):
        iso = ci.Isotope(isotope)
        decay_constant = iso.decay_const()
        half_life = iso.half_life()

        foils = self.foils(element)
        eob_activity_list = []
        for foil in foils:
            peak_data = self.retrieve_peak_data(foil, isotope)
            try:
                df = pd.read_csv(peak_data, comment="#")
                # print(df)
                delay_time, activity, unc_activity = self.calculate_activities(df, foil, decay_constant)
                if filtering:
                    df = self.filter_counts(df)
                    delay_time, activity, unc_activity = self.filter_activities_on_delay_time(delay_time, activity, unc_activity, half_life,
                                        min_half_lives=min_half_lives, max_half_lives=max_half_lives, max_activity_uncertainty=max_activity_uncertainty)
                if compartment and compartment in foil and plot_curve==False:
                    plot=True
                else:
                    plot=plot_curve
                A0, unc_A0 = self.fit_activity_first_order(delay_time, activity, unc_activity, decay_constant, plot=plot)
                data = [foil, isotope, A0, unc_A0]
                eob_activity_list.append(data)
                if plot:
                    self.savefig(element, isotope, foil)
            except Exception as e:
                if not isinstance(e, FileNotFoundError):
                    print(f"{foil}: + Error type: {type(e).__name__} – {e}")
                    # raise e
                # else:
                print('no peak data found in foil ' + foil + ' for isotope ' + isotope)
                data = [foil, isotope, 0, 0]
                eob_activity_list.append(data)
        df = pd.DataFrame(eob_activity_list, columns=['foil', 'isotope', 'fit', 'cov'])
        self.savecsv(df, element, isotope)


    def two_step_decay(self, element, isotope, isotope_parent, compartment, plot_curve, min_half_lives=None, max_half_lives=None, max_activity_uncertainty=None, independent_parent=None,):
        # one-step to a feeding curve, want total cumulative, so remove gammas before with halflife_p *10
        i = ci.Isotope(isotope); ip = ci.Isotope(isotope_parent)
        decay_constant = i.decay_const(); decay_constant_parent = ip.decay_const()
        half_life = i.half_life(); half_life_parent = ip.half_life()
        foils = self.foils(element)
        eob_activity_list = []
        eob_activity_parent_list = []
        eob_activity_parent = self.activity_file_parent(element, isotope_parent)
        # print(eob_activity_parent)
        if os.path.exists(eob_activity_parent):
            df = pd.read_csv(eob_activity_parent, comment="#")
            eob_activity_parent = df['fit']; unc_eob_activity_parent = df['cov']
        else:
            eob_activity_parent = None; unc_eob_activity_parent = None
        for i, foil in enumerate(foils):
            peak_data = self.retrieve_peak_data(foil, isotope)
            try:
                df = pd.read_csv(peak_data, comment="#")
                df = self.filter_counts(df)
                delay_time, activity, unc_activity = self.calculate_activities(df, foil, decay_constant)
                delay_time,  activity, unc_activity = self.filter_activities_on_delay_time(delay_time, activity, unc_activity, half_life, min_half_lives, max_half_lives, max_activity_uncertainty)
                if compartment and compartment in foil and plot==False:
                    plot = True
                else:
                    plot = plot_curve
                if eob_activity_parent is not None and eob_activity_parent[i]>0:
                    print("01")
                    # guess = [1e6, 1e6]
                    guess = [eob_activity_parent[i], 1e6]
                    A0, unc_A0, A0_parent, unc_A0_parent = self.fit_activities_second_order(delay_time, activity, unc_activity, decay_constant, decay_constant_parent=decay_constant_parent, guess=guess, plot=plot)
                elif eob_activity_parent is not None and eob_activity_parent[i] == 0:
                    print("02")
                    guess = [1e6]
                    A0_parent = 0; unc_A0_parent=0
                    A0, unc_A0 = self.fit_activity_first_order(delay_time, activity, unc_activity, decay_constant, plot=plot)
                elif eob_activity_parent == None:
                    print("03")
                    guess = [1e2, 1e4]
                    try:
                        print("trying to fit two parameters")
                        A0, unc_A0, A0_parent, unc_A0_parent = self.fit_activities_second_order(delay_time, activity, unc_activity, decay_constant, decay_constant_parent=decay_constant_parent, guess=guess,plot=plot)
                    except:
                        print("Failed fitting two parameters, will try first order")
                        A0_parent = 0; unc_A0_parent=0
                        A0, unc_A0 = self.fit_activity_first_order(delay_time, activity, unc_activity, decay_constant, plot=plot)
                else:
                    raise print('Invalid input')
                data = [foil, isotope, A0, unc_A0]
                eob_activity_list.append(data)
                # print(A0_parent)
                # if A0 > 0 and A0_parent is not None and A0_parent > 0:
                if A0_parent is not None and A0_parent > 0:
                    data_parent = [foil, isotope_parent, A0_parent, unc_A0_parent]
                    eob_activity_parent_list.append(data_parent)
                else:
                    data_parent = [foil, isotope_parent, 0, 0]
                    eob_activity_parent_list.append(data_parent)
                if plot:
                    self.savefig(element, isotope, foil)
            except Exception as e:
                # print("An exception occurred:", e) 
                # To get the exception name and message
                if not isinstance(e, FileNotFoundError):
                    print(f"{foil}: + Error type: {type(e).__name__} – {e}")
                    tb = traceback.extract_tb(e.__traceback__)
                    last = tb[-1]  # der feilen faktisk skjedde
                    print(f"{foil}: {type(e).__name__} – {e}")
                    print(f"File: {last.filename}, line: {last.lineno}, function: {last.name}")

                else:
                    print('no peak data found in foil ' + foil + ' for isotope ' + isotope)
                    data = [foil, isotope, 0, 0]
                    data_parent = [foil, isotope_parent, 0, 0]
                    eob_activity_list.append(data)
                    eob_activity_parent_list.append(data_parent)
        df = pd.DataFrame(eob_activity_list, columns=['foil', 'isotope', 'fit', 'cov'])
        df_parent = pd.DataFrame(eob_activity_parent_list, columns=['foil', 'isotope', 'fit', 'cov'])
        self.savecsv(df, element, isotope)
        self.savecsv_parent(df_parent, element, isotope_parent)
        # df_parent.to_csv()

    def first_order_bateman(self, t, p0, decay_constant):
        return p0 * np.exp(-decay_constant*t)
    
    def second_order_bateman(self, t, p0, decay_constant_parent, decay_constant_daughter):
        p0_parent = p0[0]; p0_daughter = p0[1]
        first_order_bateman = p0_daughter * np.exp(-decay_constant_daughter*t)
        feeding = p0_parent * decay_constant_daughter / (decay_constant_parent-decay_constant_daughter) * (np.exp(-decay_constant_daughter*t) - np.exp(-decay_constant_parent *t))
        return feeding + first_order_bateman

    def foils(self, element):
        n=14
        return [f"{element}{i:02d}" for i in range(1, n+1)]
    
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
    
    def calculate_activities(self, peak_data, foil, decay_constant):
        peak_data = peak_data.sort_values(by='start_time', ascending=True).reset_index(drop=True)
        # print(peak_data)
        E = peak_data['energy']
        counts = peak_data['counts']
        unc_counts = peak_data['unc_counts']
        intensity = peak_data['intensity']
        unc_intensity = peak_data['unc_intensity']
        efficiency = peak_data['efficiency']
        unc_efficiency = peak_data['unc_efficiency']
        # count_time = peak_data['real_time']
        count_time = peak_data['live_time']
        eob = pd.to_datetime(self.getEob(foil))
        peak_data['start_time'] = pd.to_datetime(peak_data['start_time'])
        delay_time = (peak_data['start_time'] - eob).dt.total_seconds()

        activity, unc_activity = self.activity(counts, unc_counts, efficiency, unc_efficiency, intensity, unc_intensity, count_time, decay_constant)
        activity, unc_activity = self.activity(counts, unc_counts, efficiency, unc_efficiency, intensity, unc_intensity, count_time, decay_constant)
        A = pd.Series(activity, name='activity'); unc_A = pd.Series(unc_activity, name='unc activity')
        filename = pd.Series(peak_data['filename'], name='filename')
        E = pd.Series(peak_data['energy'], name='energy')
        df = pd.concat([filename, E, delay_time/3600, A*1e-6, unc_A*1e-6], axis=1)
        print(df)

        return delay_time, activity, unc_activity

    def savefig(self, element, isotope, foil):
        # test_list.append('RUN')
        if self.saveIndependent == None:
            fig_filename = pathToFigures +  element + '/' + element + '_' + isotope + '_' + foil + '.pdf'
            title = element + ' ' + isotope + ' ' + foil
        elif self.saveIndependent == True:
            fig_filename = pathToFigures +  element + '/' + element + '_' + isotope + '_' + foil + '_ind.pdf'
            title = element + ' ' + isotope + ' ' + foil + ' - independent' 
        elif self.saveIndependent == False: 
            fig_filename = pathToFigures +  element + '/' + element + '_' + isotope + '_' + foil + '_cum.pdf'
            title = element + ' ' + isotope + ' ' + foil + ' - cumulative' 
        plt.xlabel('Time since eob (hours)')
        plt.ylabel('Activity (Bq)')
        plt.title(title)
        plt.legend()
        # print(fig_filename)
        plt.savefig(fig_filename)
        plt.show()
        # plt.close()
    
    def savecsv(self, df, element, isotope):
        # If none, should be a indirect  measurement for instance like 178Ta/178W.
        if self.saveIndependent == None:
            filename = pathToActivityFiles_isotope +  element + '_' + isotope + '.csv'
        elif self.saveIndependent == True:
            filename = pathToActivityFiles_isotope +  element + '_' + isotope + '_ind.csv'
        elif self.saveIndependent == False:
            filename = pathToActivityFiles_isotope +  element + '_' + isotope + '_cum.csv'
        # print("Saving to")
        # print(filename)
        df.to_csv(filename)


    def activity_file_parent(self, element, isotope_parent):
        if self.saveIndependentParent == None:
            filename = pathToActivityFiles_isotope +  element + '_' + isotope_parent + '.csv'
        elif self.saveIndependentParent == True:
            filename = pathToActivityFiles_isotope +  element + '_' + isotope_parent + '_ind.csv'
        elif self.saveIndependentParent == False:
            filename = pathToActivityFiles_isotope +  element + '_' + isotope_parent + '_cum.csv'
        return filename
    
    def savecsv_parent(self, df, element, isotope_parent):
        filename = self.activity_file_parent(element, isotope_parent)
        df.to_csv(filename)
        
    def activity(self, counts, unc_counts, efficiency, unc_efficiency, intensity, unc_intensity, count_time, decay_constant):
        activity = (counts*decay_constant) / (efficiency * intensity * (1-np.exp(-decay_constant*count_time)) )
        unc_activity = activity * np.sqrt((unc_counts/counts)**2  + (unc_efficiency/efficiency)**2 + (unc_intensity/intensity)**2 )
        return activity, unc_activity

    def filter_activities_on_delay_time(self, delay_time, activity, unc_activity, half_life,
                                    min_half_lives=None, max_half_lives=None, max_activity_uncertainty=None):
        delay_time = np.asarray(delay_time)
        activity = np.asarray(activity)
        unc_activity = np.asarray(unc_activity)
        mask = np.ones_like(delay_time, dtype=bool)
        if min_half_lives is not None:
            min_half_lives = min_half_lives*3600 # hours
            mask &= delay_time >= (min_half_lives)
        if max_half_lives is not None:
            mask &= delay_time <= (half_life * max_half_lives)
        if max_activity_uncertainty:
            mask &= (unc_activity / activity) < max_activity_uncertainty
        mask &= np.isfinite(activity) & np.isfinite(unc_activity) & (unc_activity > 0)

        return delay_time[mask], activity[mask], unc_activity[mask]

    def filter_counts(self, peakdata):
        max_error=0.4; min_counts=1
        mask = ((peakdata['counts'] > min_counts) & (peakdata['unc_counts'] / peakdata['counts'] < max_error))
        return peakdata[mask]

    def fit_activities_second_order(self, delay_time, activity, unc_activity, decay_constant, decay_constant_parent=None, guess=[1e5, 1e5], plot=False):
        delay_time_hours = delay_time/3600
        time = np.linspace(0, np.max(delay_time), 1000)
        popt, pcov = curve_fit(lambda t, p0, p1: self.second_order_bateman(t, (p0, p1), decay_constant_parent, decay_constant),delay_time,activity,p0=guess,sigma=unc_activity,absolute_sigma=True)
        p0_fit, p1_fit = popt  # parent and daughter initial activities (or whatever you defined)
        A0_parent = p0_fit; unc_A0_parent = np.sqrt(pcov[0,0])
        A0 = p1_fit; unc_popt = np.sqrt(pcov[1,1])
        time_hours = time/3600
        if plot:
            A_fit_feeding = self.first_order_bateman(time, popt[0], decay_constant_parent)
            A_fit = self.second_order_bateman(time, (popt[0], popt[1]), decay_constant_parent, decay_constant)
            plt.errorbar(delay_time_hours, activity, color='darkolivegreen', linewidth=0.001,yerr=unc_activity, elinewidth=0.5, ecolor='k', capthick=0.5,marker='.', label='activity')   # cap thickness for error bar color='blue')
            plt.plot(time_hours, A_fit, color='tan', linewidth=0.9, label='fit')
            plt.plot(time_hours, A_fit_feeding, color='dodgerblue', linewidth=0.9, label='feeding')
            plt.errorbar(0, A0, marker='*', yerr=unc_popt, label=r'eob: %.f ($\pm %.f$) Bq' %(A0,unc_popt), color='maroon')
        return A0, unc_popt, A0_parent, unc_A0_parent
    
    def fit_activity_first_order(self, delay_time, activity, unc_activity, decay_constant, guess=[1e5], plot=False):
        delay_time_hours = delay_time/3600
        time = np.linspace(0, np.max(delay_time), 1000) 
        popt, pcov = curve_fit(lambda t, p0: self.first_order_bateman(t, p0, decay_constant),delay_time,activity,p0=guess,sigma=unc_activity,absolute_sigma=True)
        unc_popt = np.sqrt(np.diagonal(pcov))
        A0 = self.first_order_bateman(0, *popt, decay_constant)
        time_hours = time/3600
        print(A0, type(A0), np.shape(A0))
        print(unc_popt, type(unc_popt), np.shape(unc_popt))
        if isinstance(unc_popt, np.ndarray):
            unc = unc_popt[0]
        else:
            unc = unc_popt
        if plot:
            A_fit = self.first_order_bateman(time, *popt, decay_constant)
            unc_fit_plus = self.first_order_bateman(time,*(popt+unc_popt), decay_constant); unc_fit_minus=self.first_order_bateman(time,*(popt-unc_popt), decay_constant)
            plt.errorbar(delay_time_hours, activity, color='darkolivegreen', linewidth=0.001,yerr=unc_activity, elinewidth=0.5, ecolor='k', capthick=0.5,marker='.', label='activity')   # cap thickness for error bar color='blue')
            plt.plot(time_hours, A_fit, color='tan', linewidth=0.9, label='fit')
            plt.fill_between(time_hours, unc_fit_plus, unc_fit_minus, color='tan', alpha=0.1)
            plt.errorbar(0, A0, marker='*', yerr=unc_popt, label=r'eob: %.f ($\pm %.f$) Bq' %(A0,unc), color='maroon')
        return A0, unc_popt[0]

    def uncertainty_A0(self, pcov, decay_constant_daughter, decay_constant_parent,  time):
        # TODO
        try:
            deriv_Ad_Ap = decay_constant_daughter/(decay_constant_parent-decay_constant_daughter)*(np.exp(-decay_constant_daughter*time)-np.exp(-decay_constant_parent*time))
            deriv_Ad_Ad = np.exp(-decay_constant_daughter*time)
            J = np.array((deriv_Ad_Ap, deriv_Ad_Ad)) #Jacobian
            dA0_d = np.sqrt( np.dot(np.dot(J,pcov),J.T ) )
        except:
            print('Missing decay constant parent?')
        return dA0_d

    def retrieve_peak_data(self, foil, isotope):
        str = foil + '_' + isotope + '_gammas.csv'
        peak_data = pathToPeakFiles_isotope + str
        if not os.path.exists(peak_data):
            str = foil + '_' + isotope + 'g_gammas.csv'
            peak_data = pathToPeakFiles_isotope + str
        return peak_data
