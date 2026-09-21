import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge
from get_variables import *
import curie as ci

class WeightedAverageBeamCurrent:

    def __init__(self, number_of_compartments, number_of_monitor_foils, monitor_reactions_per_foil, wa, stack, list_of_isotopes):
        self.number_of_compartments = number_of_compartments # 7, 
        self.number_of_monitor_foils = number_of_monitor_foils # 2
        self.monitor_reactions_per_foil = monitor_reactions_per_foil # np.array([6, 6])  for cupper and nickel, respectively
        self.wa = wa # WeightedAverageFlux file
        self.stack = stack
        # self.list_of_isotopes = list_of_isotopes # ['Cu_63Zn', 'Cu_62Zn', ...]
        self.list_of_isotopes = list_of_isotopes   #[params_Cu_65Zn, params_Cu_63Zn, params_Cu_62Zn, params_Cu_58CO, params_Cu_56CO, params_Ni_57NI]
        self.irradiation_time, self.unc_irradiation_time = get_irradiation_time(stack) 

    def weighted_average_beam_current(self, csv_filename=None): # ALL ANDREW WORK WHAT IS GOING ON??? haha
        matrix_A0, matrix_sigma_A0, matrix_lambda_, matrix_mass_density, matrix_sigma_mass_density, matrix_reaction_integral, matrix_uncertainty_integral, matrix_irr_time, matrix_sigma_irr_time = self.reshaping_parameters()
        number_of_monitor_reactions = 0
        print(self.list_of_isotopes)
        submatrix_lower_indices = np.zeros(self.number_of_monitor_foils)
        submatrix_upper_indices = np.zeros(self.number_of_monitor_foils)


        for i in range(0, self.number_of_monitor_foils):
            # print(i)
            submatrix_lower_indices[i] = number_of_monitor_reactions
            number_of_monitor_reactions += self.monitor_reactions_per_foil[i]
            submatrix_upper_indices[i] = number_of_monitor_reactions
            # print(submatrix_lower_indices, submatrix_upper_indices)

        submatrix_lower_indices=submatrix_lower_indices.astype(int)
        submatrix_upper_indices=submatrix_upper_indices.astype(int)
        

        # Set up correlation matrices
	    # Lambda is completely uncorrelated, except on the diagonal
        corr_lambda = np.zeros((number_of_monitor_reactions,number_of_monitor_reactions))
	    # Areal density is completely uncorrelated, except within one foil's submatrix
        corr_areal_density = np.zeros((number_of_monitor_reactions,number_of_monitor_reactions))
	    # Reaction integral is completely uncorrelated, except within one foil's submatrix
        corr_reaction_integral = np.zeros((number_of_monitor_reactions,number_of_monitor_reactions))
	    # Irradiation length is completely correlated (same for all foils)
        corr_t_irradiation = np.ones((number_of_monitor_reactions,number_of_monitor_reactions))
	    # EoB activitoes are partially uncorrelated (similar subset of efficiencies)
        corr_EoB_activities = 0.3 * np.ones((number_of_monitor_reactions,number_of_monitor_reactions))    #just set to 0.3 since we do not have MC simulations
    
        # Set up lists to hold output data
        foil_index = []
        weighted_average = []
        unc_weighted_average = []
        percent_unc_weighted_average = []

        # Get correlation submtarix for each monitor foil - n x n, where n= # of reactions per foil
        for i in range(self.number_of_monitor_foils): # Monitor reactions per foil [3,3,1]
            submatrix = np.ones((self.monitor_reactions_per_foil[i], self.monitor_reactions_per_foil[i]))
            corr_areal_density[submatrix_lower_indices[i]:submatrix_upper_indices[i], submatrix_lower_indices[i]:submatrix_upper_indices[i]] = submatrix
            corr_reaction_integral[submatrix_lower_indices[i]:submatrix_upper_indices[i], submatrix_lower_indices[i]:submatrix_upper_indices[i]] = 0.3*submatrix
        # Ensure all diagonal elements are still ones
        np.fill_diagonal(corr_lambda,1)
        np.fill_diagonal(corr_areal_density,1)
        np.fill_diagonal(corr_reaction_integral,1)
        np.fill_diagonal(corr_t_irradiation,1)
        np.fill_diagonal(corr_EoB_activities,1)

        # Loop over all beam positions
        number_of_energies = len(matrix_mass_density)

        ## Hold curents as we go along...
        currents = np.zeros((number_of_energies,number_of_monitor_reactions))
        unc_currents = np.zeros((number_of_energies,number_of_monitor_reactions))
        function_dictionary = {'0':self.dIdA0, '1':self.dIdRhoDr, '2':self.dIdLambda, '3':self.dIdTIrradiation, '4':self.dIdIntegral}

        # unc_rxn_integral = uncertainty_integral * reaction_integral#rxn_int * 	percent_rn_uncertainties      #rxn = reactions
        unc_rxn_integral = matrix_uncertainty_integral * matrix_reaction_integral   #rxn_int * 	percent_rn_uncertainties      #rxn = reactions
        for i_energy in range(0, number_of_energies):
            # print("energy: " + str(i_energy))
        #    # Get nonzero entries in A0:
            # nonzero_indices = np.nonzero(matrix_A0[i_energy,:])
            nonzero_indices = np.where((matrix_A0[i_energy,:] > 0) & (matrix_reaction_integral[i_energy,:] > 0))
            ad = matrix_mass_density[i_energy,:]
            unc_ad = matrix_sigma_mass_density[i_energy,:]
            A0 = matrix_A0[i_energy,:]
            unc_A0 = matrix_sigma_A0[i_energy,:]
            rxn_int = matrix_reaction_integral[i_energy,:]
            # print("rxn int: "+ str(rxn_int))
            delta_t = np.ones(number_of_monitor_reactions) * matrix_irr_time[i_energy]
            unc_delta_t = np.ones(number_of_monitor_reactions) * matrix_sigma_irr_time[i_energy]
            unc_rxn_int = unc_rxn_integral[i_energy,:]
            loop_lambdas = matrix_lambda_
            uncertainty_lambdas = loop_lambdas * 0.001

            # print(len(np.transpose(nonzero_indices)), number_of_monitor_reactions)
            if len(np.transpose(nonzero_indices)) == number_of_monitor_reactions:
                # No nonzero indices!!!
                # Keep normal correlation matrices
                loop_corr_lambda = corr_lambda
                loop_corr_areal_density = corr_areal_density
                loop_corr_reaction_integral = corr_reaction_integral
                loop_corr_t_irradiation = corr_t_irradiation
                loop_corr_EoB_activities = corr_EoB_activities

            else:
                # Some nonzero indices
                # Find which indices are missing!
                temp3 = np.asarray(nonzero_indices[0])
                temp4 = np.array(range(0, number_of_monitor_reactions))
                disjoint_indices = np.setdiff1d(temp4,temp3,assume_unique=False).tolist()
                #print('disjoint indices: ', disjoint_indices)
                # Delete rows and columns in correlation matries of disjoint indices
                # gen = (x for x in xyz if x not in a)
                if len(disjoint_indices) != 1:
                    loop_corr_lambda = np.delete(corr_lambda,np.array(disjoint_indices),0)
                    loop_corr_lambda = np.delete(loop_corr_lambda,np.array(disjoint_indices),1)
                    loop_corr_areal_density = np.delete(corr_areal_density,np.array(disjoint_indices),0)
                    loop_corr_areal_density = np.delete(loop_corr_areal_density,np.array(disjoint_indices),1)
                    loop_corr_reaction_integral = np.delete(corr_reaction_integral,np.array(disjoint_indices),0)
                    loop_corr_reaction_integral = np.delete(loop_corr_reaction_integral,np.array(disjoint_indices),1)
                    loop_corr_t_irradiation = np.delete(corr_t_irradiation,np.array(disjoint_indices),0)
                    loop_corr_t_irradiation = np.delete(loop_corr_t_irradiation,np.array(disjoint_indices),1)
                    loop_corr_EoB_activities = np.delete(corr_EoB_activities,np.array(disjoint_indices),0)
                    loop_corr_EoB_activities = np.delete(loop_corr_EoB_activities,np.array(disjoint_indices),1)
                else:
                    for disjoint_index in disjoint_indices:
                        loop_corr_lambda = np.delete(corr_lambda,disjoint_index,0)
                        loop_corr_lambda = np.delete(loop_corr_lambda,disjoint_index,1)
                        loop_corr_areal_density = np.delete(corr_areal_density,disjoint_index,0)
                        loop_corr_areal_density = np.delete(loop_corr_areal_density,disjoint_index,1)
                        loop_corr_reaction_integral = np.delete(corr_reaction_integral,disjoint_index,0)
                        loop_corr_reaction_integral = np.delete(loop_corr_reaction_integral,disjoint_index,1)
                        loop_corr_t_irradiation = np.delete(corr_t_irradiation,disjoint_index,0)
                        loop_corr_t_irradiation = np.delete(loop_corr_t_irradiation,disjoint_index,1)
                        loop_corr_EoB_activities = np.delete(corr_EoB_activities,disjoint_index,0)
                        loop_corr_EoB_activities = np.delete(loop_corr_EoB_activities,disjoint_index,1)
            # print(loop_corr_areal_density.shape)
            #print('beam_current inputs: ', A0, ad, loop_lambdas[i_energy,:], delta_t, rxn_int)
            ## A0, areal_density, decay_constant, irr_time, reaction_integral):
            temp_currents =  self.beam_current(A0, ad, loop_lambdas[i_energy,:], delta_t, rxn_int)
            #print(temp_currents)
            #currents[i_energy, :] =  temp_currents
            #print('temp_currents: ', temp_currents)
            #print('unc_beam_current inputs: ', A0, ad, loop_lambdas[i_energy,:], delta_t, rxn_int, unc_A0, unc_ad, uncertainty_lambdas[i_energy,:], unc_delta_t, unc_rxn_int)
            unc_temp_currents = self.unc_beam_current_approx(A0, ad, loop_lambdas[i_energy,:], delta_t, rxn_int, unc_A0, unc_ad, uncertainty_lambdas[i_energy,:], unc_delta_t, unc_rxn_int)
            unc_currents[i_energy,:] = unc_temp_currents

            value_array = np.array([A0, ad, loop_lambdas[0], delta_t, rxn_int])
            uncertainty_array = np.array([unc_A0, unc_ad, uncertainty_lambdas[0], unc_delta_t, unc_rxn_int])
            # ENDRET:
            # value_array = np.array([A0, ad, loop_lambdas[i_energy], delta_t, rxn_int])
            # uncertainty_array = np.array([unc_A0, unc_ad, uncertainty_lambdas[i_energy], unc_delta_t, unc_rxn_int])
            
            correlation_array = np.array([loop_corr_EoB_activities, loop_corr_areal_density, loop_corr_lambda, loop_corr_t_irradiation, loop_corr_reaction_integral])


            # Set up covariance matrix for current energy position
            cov = np.zeros((len(nonzero_indices[0]),len(nonzero_indices[0])))

            # NaN handling - replace range(0,number_of_monitor_reactions) with indices of nonzero elements of A0?
            # Fill correlation matrices
            for i_index,i_element in enumerate(nonzero_indices[0]):
                # if i_energy == 5:
                #     # print(f"\nFunction key: {dict_key}")

                #     print("INPUTS i:")
                #     print("A0:", A0[i_element])
                #     print("ad:", ad[i_element])
                #     print("lambda:", loop_lambdas[0, i_element])
                #     print("delta_t:", delta_t[i_element])
                #     print("rxn_int:", rxn_int[i_element])
                for j_index,j_element in enumerate(nonzero_indices[0]):
                    for dict_index,dict_key in enumerate(function_dictionary):
                        if rxn_int[i_element]==0:    # FIX FOR .... 62ZN
                            rxn_int[i_element] = 1e-30
                        if rxn_int[j_element]==0:
                            rxn_int[j_element] = 1e-30

                        dIdxi = function_dictionary[dict_key](A0[i_element], ad[i_element], loop_lambdas[0,i_element], delta_t[i_element], rxn_int[i_element])
                        # if np.isnan(dIdxi):
                        #     print(rxn_int[i_element])
                        # if np.isnan(dIdxi):
                        #     print(function_dictionary[dict_key])
                        #     print(dIdxi)
                        dIdxj = function_dictionary[dict_key](A0[j_element], ad[j_element], loop_lambdas[0,j_element], delta_t[j_element], rxn_int[j_element])
                        # if np.isnan(dIdxj):
                        #     print(i_index, j_index)
                        #     print(function_dictionary[dict_key])
                        #     print(A0[j_element], ad[j_element], loop_lambdas[0,j_element], delta_t[j_element], rxn_int[j_element])
                        #     print("*****")
                        cov[i_index,j_index] += dIdxi * uncertainty_array[int(dict_key),i_element] *  correlation_array[int(dict_key),i_index,j_index] *  uncertainty_array[int(dict_key),j_element] * dIdxj
                        # if np.isnan(dIdxi) or np.isnan(dIdxj):
                        # if i_energy==5:
                        #     print("***")
                        #     print(i_index, j_index)
                        #     print(cov[i_index,j_index])
                        #     print(dIdxi * uncertainty_array[int(dict_key),i_element] *  correlation_array[int(dict_key),i_index,j_index] *  uncertainty_array[int(dict_key),j_element] * dIdxj)
                        #     if np.isnan(dIdxj):
                        #         print("INPUTS j:")
                        #         print("A0:", A0[j_element])
                        #         print("ad:", ad[j_element])
                        #         print("lambda:", loop_lambdas[0, j_element])
                        #         print("delta_t:", delta_t[j_element])
                        #         print("rxn_int:", rxn_int[j_element])
            
            
            # print("Final covariance matrix energy: " + str(i_energy))
            # print( cov)
            inverted_covariance = np.linalg.inv(cov)
            # print(inverted_covariance)
            numerator = 0.0
            denominator = 0.0

            sigma = np.sqrt(np.diag(cov))
            corr = cov / np.outer(sigma, sigma)

            #print(corr)
            one_vector = np.ones(len(cov))

            weights = inverted_covariance @ corr
            weights = weights / (corr @ inverted_covariance @ corr)
            print("weights")
            print(weights)
            print("****")

            for i_index,i_element in enumerate(nonzero_indices[0]):
                for j_index,j_element in enumerate(nonzero_indices[0]):
                    numerator += temp_currents[j_element] * inverted_covariance[i_index,j_index]
                    denominator +=  inverted_covariance[i_index,j_index]

            weighted_average_current = numerator/denominator
            uncertainty_weighted_average_current = np.sqrt(1.0/denominator)


           ###print("weighted_average_current: ",weighted_average_current, " +/- ",uncertainty_weighted_average_current, " nA     (", 100*uncertainty_weighted_average_current/weighted_average_current ," %)")

           # Append values for current energy
            foil_index.append(i_energy)
            weighted_average.append(weighted_average_current)
            unc_weighted_average.append(uncertainty_weighted_average_current)
            percent_unc_weighted_average.append(100*uncertainty_weighted_average_current/weighted_average_current)

        outfile = np.stack((np.transpose(foil_index),np.transpose(weighted_average),np.transpose(unc_weighted_average),np.transpose(percent_unc_weighted_average)), axis=-1)
        df = pd.DataFrame(outfile,columns=["foil_index", "Average Current (nA)", "Uncertainty in Average Current (nA)", "% Uncertainty"])
        path = os.getcwd() + '/generatedfiles/beamcurrent/true_weighted_average_beam_current/'
        # df = df.iloc[::-1].reset_index(drop=True)
        print(df)
        if csv_filename:
            print("should save to: " + path + csv_filename)
            df.to_csv(path + csv_filename + '.csv')

        # return weighted_average[::-1], unc_weighted_average[::-1 ] #returning reversed lists
        return weighted_average, unc_weighted_average #returning reversed lists
        
    def flux_weighted_average_energy(self, element='Ta'):
        energy, unc_left, unc_right = weighted_average_beam_energy(element=element, stack=self.stack)
        return energy, unc_left, unc_right
    
    def plot_weighted_average(self, color='darkred', note='wabc'):
        energy, unc_left, unc_right = weighted_average_beam_energy(element='Ta')
        if self.stack == 'stack_30_MeV':
            energy = energy[7:]; unc_left=unc_left[7:]; unc_right = unc_right[7:]
        else:
            energy = energy[0:7]; unc_left=unc_left[0:7]; unc_right = unc_right[0:7]
        wabc, unc_wabc = self.weighted_average_beam_current()
        wabc = np.nan_to_num(wabc); unc_wabc = np.nan_to_num(unc_wabc)
        # print(wabc)
        plt.errorbar(energy, wabc, marker='*', color=color,linewidth=0.0001,
            xerr=[unc_left, unc_right], yerr=unc_wabc, elinewidth=1.0, capthick=3.0, capsize=3.0,
            label=note, linestyle='none')
        
    def plot_beam_current_isotope(self, element, isotope_params, color, label, symbol='.'):
        params = getattr(self, isotope_params)
        # print(params())
        # TODO what is std eob acitivyt?
        eob_activitiy, std_eob_activity, decay_constant, areal_density, unc_areal_density, flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section = params()
        irr_time = np.ones(self.number_of_compartments)*self.irradiation_time
        unc_irr_time = np.ones(self.number_of_compartments)*self.unc_irradiation_time
        beam_current = self.beam_current(eob_activitiy, areal_density, decay_constant, irr_time, flux_weighted_average_cross_section)
        unc_beam_current = self.std_beam_current(beam_current, eob_activitiy, std_eob_activity, areal_density, unc_areal_density, irr_time, unc_irr_time, decay_constant, 0.01, flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section)
        energy, unc_left, unc_right = self.flux_weighted_average_energy(element)
        mask = beam_current>0
        # print(element)
        # print(beam_current)
        plt.errorbar(energy[mask], beam_current[mask], color=color, marker=symbol, ls = '',linewidth=0.001, xerr=[unc_left[mask], unc_right[mask]], yerr=np.abs(unc_beam_current[mask]), elinewidth=0.8, capthick=0.8, capsize=3.0,label=label)

    def beam_current(self, A0, areal_density, decay_constant, irr_time, reaction_integral):
        # print(A0.shape)
        # print("A0")
        # print(A0)
        # print("areal_density")
        # print(areal_density)
        # print("decay_constant")
        # print(decay_constant)
        # print("irr_time")
        # print(irr_time)
        # print("reaction_integral")
        # print(reaction_integral)
        return (elementary_charge*1e9 * A0) / (areal_density * (1 - np.exp(-decay_constant*irr_time)) * reaction_integral)

    def std_beam_current(self, beam_current, eob_activity, unc_eob_activity, areal_density, unc_areal_density, irradiation_time, unc_irradiation_time,
                         decay_constant, unc_decay_constant, flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section):
        return beam_current * np.sqrt(
              (unc_eob_activity / eob_activity)**2
            + (unc_areal_density/areal_density)**2 
            + (unc_irradiation_time/irradiation_time)**2 
            # + (unc_decay_constant/decay_constant)**2 
            + (unc_flux_weighted_average_cross_section/flux_weighted_average_cross_section)**2)
        
    def dIdA0(self, A0, areal_density, decay_constant, t_irradiation, reaction_integral):
        delta_x = 1E-8 * A0
        return ((self.beam_current(A0 + (delta_x/2), areal_density, decay_constant, t_irradiation, reaction_integral) - self.beam_current(A0 - (delta_x/2), areal_density, decay_constant, t_irradiation, reaction_integral)) / delta_x)
    def dIdRhoDr(self, A0, areal_density, decay_constant, t_irradiation, reaction_integral):
        delta_x = 1E-8 * areal_density
        return ((self.beam_current(A0, areal_density + (delta_x/2), decay_constant, t_irradiation, reaction_integral) - self.beam_current(A0, areal_density - (delta_x/2), decay_constant, t_irradiation, reaction_integral)) / delta_x)
    def dIdLambda(self, A0, areal_density, decay_constant, t_irradiation, reaction_integral):
        delta_x = 1E-8 * decay_constant
        return ((self.beam_current(A0, areal_density, decay_constant + (delta_x/2), t_irradiation, reaction_integral) - self.beam_current(A0, areal_density, decay_constant - (delta_x/2), t_irradiation, reaction_integral)) / delta_x)
    def dIdTIrradiation(self, A0, areal_density, decay_constant, t_irradiation, reaction_integral):
        delta_x = 1E-8 * t_irradiation
        return ((self.beam_current(A0, areal_density, decay_constant, t_irradiation + (delta_x/2), reaction_integral) - self.beam_current(A0, areal_density, decay_constant, t_irradiation - (delta_x/2), reaction_integral)) / delta_x)
    def dIdIntegral(self, A0, areal_density, decay_constant, t_irradiation, reaction_integral):
        delta_x = 1E-8 * reaction_integral
        return ((self.beam_current(A0, areal_density, decay_constant, t_irradiation, reaction_integral + (delta_x/2)) - self.beam_current(A0, areal_density, decay_constant, t_irradiation, reaction_integral - (delta_x/2))) / delta_x)

    def unc_beam_current_approx(self, A0, rho_dr, decay_constant, t_irradiation, reaction_integral, unc_A0, unc_rho_dr, unc_decay_constant, unc_t_irradiation, unc_reaction_integral):
        approx_error = np.sqrt(np.power(self.dIdA0(A0, rho_dr, decay_constant, t_irradiation, reaction_integral) * unc_A0,2) +
        np.power(self.dIdRhoDr(A0, rho_dr, decay_constant, t_irradiation, reaction_integral) * unc_rho_dr,2) +
        np.power(self.dIdLambda(A0, rho_dr, decay_constant, t_irradiation, reaction_integral) * unc_decay_constant,2) +
        np.power(self.dIdTIrradiation(A0, rho_dr, decay_constant, t_irradiation, reaction_integral) * unc_t_irradiation,2) +
        np.power(self.dIdIntegral(A0, rho_dr, decay_constant, t_irradiation, reaction_integral) * unc_reaction_integral,2))
        return approx_error

    def parameters(self):
        parameters = {}
        for isotope in self.list_of_isotopes:
                params = getattr(self, isotope)
                parameters[isotope] = params()
        return parameters
    
    def reshaping_parameters(self):
        irradiation_time, unc_irradiation_time = get_irradiation_time(self.stack)
        matrix_A0, matrix_sigma_A0, matrix_lambda_, matrix_mass_density, matrix_sigma_mass_density, matrix_reaction_integral, matrix_uncertainty_integral, matrix_irr_time, matrix_sigma_irr_time = self.empty_matrix()
        list_of_parameters =  np.array(list(self.parameters().values()), dtype=object)
   
        n = len(list_of_parameters)
        A0 = np.zeros(n); sigma_A0 = np.zeros(n)
        decay_constant = np.zeros(n)
        mass_density = np.zeros(n); unc_mass_density =  np.zeros(n)
        reaction_integral = np.zeros(n); unc_reaction_integral = np.zeros(n)
        matrix_irr_time = np.transpose(np.ones(self.number_of_compartments)*irradiation_time)
        matrix_sigma_irr_time = np.transpose(np.ones(self.number_of_compartments)*unc_irradiation_time)
        for i in range(len(list_of_parameters)):
            A0 = list_of_parameters[i,0]
            sigma_A0 = list_of_parameters[i,1]
            lambda_ = list_of_parameters[i,2]
            mass_density = list_of_parameters[i,3]
            sigma_mass_density = list_of_parameters[i,4]
            reaction_integral = list_of_parameters[i,5]
            uncertainty_integral = list_of_parameters[i,6]
            
            matrix_lambda_[:, i] = lambda_
            matrix_A0[:,i] = A0
            matrix_sigma_A0[:,i] = sigma_A0
            matrix_mass_density[:,i] = mass_density
            matrix_sigma_mass_density[:,i] = sigma_mass_density
            matrix_reaction_integral[:,i]  = reaction_integral
            matrix_uncertainty_integral[:,i]  = uncertainty_integral
        
        rows = matrix_A0.shape[0]
        columns = matrix_A0.shape[1]
        for i in range(rows):
            for j in range(columns):
                if np.isinf(matrix_sigma_A0[i,j]):
                    matrix_sigma_A0[i,j]=0

        # print(matrix_A0)
        # print("***")
        # print(matrix_sigma_A0)
        # print("***")
        # print(matrix_lambda_)
        # print("***")
        # print(matrix_mass_density)
        # print("***")
        # print(matrix_sigma_mass_density)
        # print("***")
        # print(matrix_reaction_integral)
        # print("***")
        # print(matrix_uncertainty_integral)
        # print("***")
        # print(matrix_irr_time)
        # print("***")
        # print(matrix_sigma_irr_time)
        # print("***")
        return matrix_A0, matrix_sigma_A0, matrix_lambda_, matrix_mass_density, matrix_sigma_mass_density, matrix_reaction_integral, matrix_uncertainty_integral, matrix_irr_time, matrix_sigma_irr_time

    def empty_matrix(self):
        number_of_reactions = np.sum(self.monitor_reactions_per_foil)
        matrix_A0 = np.zeros(( self.number_of_compartments, number_of_reactions ))
        matrix_sigma_A0 = np.zeros(( self.number_of_compartments, number_of_reactions))
        matric_decay_constant = np.zeros(( self.number_of_compartments, number_of_reactions))
        matrix_mass_density = np.zeros(( self.number_of_compartments, number_of_reactions))
        matrix_sigma_mass_density = np.zeros(( self.number_of_compartments, number_of_reactions))
        matrix_reaction_integral = np.zeros(( self.number_of_compartments, number_of_reactions))
        matrix_uncertainty_integral = np.zeros(( self.number_of_compartments, number_of_reactions))
        matrix_irr_time = np.zeros(( self.number_of_compartments, 1))
        matrix_sigma_irr_time = np.zeros(( self.number_of_compartments, 1))
        return matrix_A0, matrix_sigma_A0, matric_decay_constant, matrix_mass_density, matrix_sigma_mass_density, matrix_reaction_integral, matrix_uncertainty_integral, matrix_irr_time, matrix_sigma_irr_time

    def assemble_parameters(self, stack, element, isotope, independent, wa):
        foils, areal_density, unc_areal_density = areal_density_from_files(element, stack)
        eob_activitiy, std_eob_activity = eob_activity_manually(element, isotope, independent=independent, stack=stack)
        energy, flux = wa.get_flux_energy_stack(element)
        flux_weighted_average_energy, unc_energy_left, unc_energy_right = wa.flux_weighted_average_energy(energy, flux)
        flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section = wa.monitor_flux_weighted_average_cross_section(element, isotope)
        decay_constant = ci.Isotope(isotope).decay_const()
        return eob_activitiy, std_eob_activity, decay_constant, areal_density, unc_areal_density, flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section

    def params_Cu_65Zn(self):
        return self.assemble_parameters(self.stack, 'Cu', '65ZN', None, self.wa)
    def params_Cu_63Zn(self):
        return self.assemble_parameters(self.stack, 'Cu', '63ZN', None, self.wa)
    def params_Cu_62Zn(self):
        return self.assemble_parameters(self.stack, 'Cu', '62ZN', None, self.wa)
    def params_Cu_61Cu(self):
        return self.assemble_parameters(self.stack, 'Cu', '61CU', None, self.wa)
    def params_Cu_58Co(self):
        return self.assemble_parameters(self.stack, 'Cu', '58CO', None, self.wa)
    def params_Cu_56Co(self):
        return self.assemble_parameters(self.stack, 'Cu', '56CO', None, self.wa)
    def params_Ni_57Ni(self):
        return self.assemble_parameters(self.stack, 'Ni', '57NI', None, self.wa)
    def params_Ni_56Ni(self):
        return self.assemble_parameters(self.stack, 'Ni', '56NI', None, self.wa)
    def params_Ni_58Co(self):
        return self.assemble_parameters(self.stack, 'Ni', '58CO', None, self.wa)
    def params_Ni_57Co(self):
        return self.assemble_parameters(self.stack, 'Ni', '57CO', None, self.wa)
    def params_Ni_56Co(self):
        return self.assemble_parameters(self.stack, 'Ni', '56CO', None, self.wa)
    def params_Ni_55Co(self):
        return self.assemble_parameters(self.stack, 'Ni', '55CO', None, self.wa)



# flux_30_after = 'TaSn_stack_30MeV_dp_1.010%_fluxes.csv' ; full_stack_30 = 'TaSn_stack_30MeV_dp_1.010%.csv'
# flux_55_after = 'TaSn_stack_55MeV_dp_1.020%_fluxes.csv' ; full_stack_55 = 'TaSn_stack_55MeV_dp_1.010%.csv'
# flux_55_after2 = 'TaSn_stack_55MeV_dp_1.010%_fluxes.csv' ; full_stack_55 = 'TaSn_stack_55MeV_dp_1.010%.csv'
# flux_30_before = 'TaSn_stack_30MeV_dp_1.000%_fluxes.csv'; full_stack_30 = 'TaSn_stack_30MeV_dp_1.000%.csv'
# flux_55_before = 'TaSn_stack_55MeV_dp_1.000%_fluxes.csv'; full_stack_55 = 'TaSn_stack_55MeV_dp_1.000%.csv'

# wa_55_before = WeightedAverageFlux(stack_fluxes = flux_55_before, full_stack = full_stack_55, stack = 'stack_55_MeV')
# wa_55_after = WeightedAverageFlux(stack_fluxes = flux_55_after, full_stack = full_stack_55, stack = 'stack_55_MeV')
# wa_55_after2 = WeightedAverageFlux(stack_fluxes = flux_55_after2, full_stack = full_stack_55, stack = 'stack_55_MeV')
# wa_30_before = WeightedAverageFlux(stack_fluxes = flux_30_before, full_stack = full_stack_30, stack = 'stack_30_MeV')
# wa_30_after = WeightedAverageFlux(stack_fluxes = flux_30_after, full_stack = full_stack_30, stack = 'stack_30_MeV')

# # wabc = WeightedAverageBeamCurrent(7, [6,6], wa_30_after, 'stack_30_MeV', ['params_Cu_65Zn', 'params_Cu_63Zn', 'params_Ni_57Ni'] )
# params_all = [
#     'params_Cu_65Zn', 'params_Cu_63Zn', 'params_Cu_62Zn', 'params_Cu_61Cu', 'params_Cu_58Co', 'params_Cu_56Co', 
#     'params_Ni_57Ni', 'params_Ni_56Ni', 'params_Ni_58Co', 'params_Ni_57Co', 'params_Ni_56Co', 'params_Ni_55Co'
#     ]

# params_55_stack = [
#     'params_Cu_65Zn', 'params_Cu_63Zn', 'params_Cu_62Zn', 'params_Cu_61Cu', 'params_Cu_58Co',
#     'params_Ni_57Ni', 'params_Ni_56Ni', 'params_Ni_58Co', 'params_Ni_57Co', 'params_Ni_56Co', 'params_Ni_55Co'
# ]

# params_30_stack = [
#     'params_Cu_65Zn', 'params_Cu_63Zn', 'params_Cu_62Zn', 'params_Cu_61Cu',
#     'params_Ni_57Ni', 'params_Ni_56Ni', 'params_Ni_58Co', 'params_Ni_57Co', 'params_Ni_56Co', 'params_Ni_55Co'
# ]

# wabc_30_after = WeightedAverageBeamCurrent(7, 2, np.array([4,6]), wa_30_after, 'stack_30_MeV', params_30_stack )

# wabc_55_after = WeightedAverageBeamCurrent(7, 2, np.array([6,6]), wa_55_after2, 'stack_55_MeV', params_all )


# wabc_55_after.average_beam_current('wabc_55MeV_dp1_01.csv')
# wabc_30_after.average_beam_current('wabc_30MeV_dp1_01.csv')
