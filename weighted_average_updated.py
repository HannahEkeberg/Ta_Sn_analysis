import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge
from get_variables import *
import curie as ci


flux_30_after = 'TaSn_stack_30MeV_dp_1.010%_fluxes.csv' ; full_stack_30 = 'TaSn_stack_30MeV_dp_1.010%.csv'
flux_55_after = 'TaSn_stack_55MeV_dp_1.020%_fluxes.csv' ; full_stack_55 = 'TaSn_stack_55MeV_dp_1.010%.csv'
flux_55_after2 = 'TaSn_stack_55MeV_dp_1.010%_fluxes.csv' ; full_stack_55 = 'TaSn_stack_55MeV_dp_1.010%.csv'
flux_30_before = 'TaSn_stack_30MeV_dp_1.000%_fluxes.csv'; full_stack_30 = 'TaSn_stack_30MeV_dp_1.000%.csv'
flux_55_before = 'TaSn_stack_55MeV_dp_1.000%_fluxes.csv'; full_stack_55 = 'TaSn_stack_55MeV_dp_1.000%.csv'

wa_55_before = WeightedAverageFlux(stack_fluxes = flux_55_before, full_stack = full_stack_55, stack = 'stack_55_MeV')
wa_55_after = WeightedAverageFlux(stack_fluxes = flux_55_after, full_stack = full_stack_55, stack = 'stack_55_MeV')
wa_55_after2 = WeightedAverageFlux(stack_fluxes = flux_55_after2, full_stack = full_stack_55, stack = 'stack_55_MeV')
wa_30_before = WeightedAverageFlux(stack_fluxes = flux_30_before, full_stack = full_stack_30, stack = 'stack_30_MeV')
wa_30_after = WeightedAverageFlux(stack_fluxes = flux_30_after, full_stack = full_stack_30, stack = 'stack_30_MeV')



def get_parameters(stack, element, isotope, independent, wa):
    foils, areal_density, unc_areal_density = areal_density_from_files(element, stack)
    eob_activitiy, std_eob_activity = eob_activity_manually(element, isotope, independent=independent, stack=stack)
    energy, flux = wa.get_flux_energy_stack(element)
    flux_weighted_average_energy, unc_energy_left, unc_energy_right = wa.flux_weighted_average_energy(energy, flux)
    flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section = wa.monitor_flux_weighted_average_cross_section(element, isotope)
    decay_constant = ci.Isotope(isotope).decay_const()
    return eob_activitiy, std_eob_activity, decay_constant, areal_density, unc_areal_density, flux_weighted_average_cross_section, unc_flux_weighted_average_cross_section


def reshaping_parameters(stack, wa):
    params_Cu_65Zn = get_parameters(stack, 'Cu', '65ZN', True, wa)
    params_Cu_63Zn = get_parameters(stack, 'Cu', '63ZN', True, wa)
    params_Cu_62Zn = get_parameters(stack, 'Cu', '62ZN', True, wa)
    # params_Cu_58CO = get_parameters(stack, 'Cu', '58CO', False, wa)
    # params_Cu_56CO = get_parameters(stack, 'Cu', '56CO', False, wa)
    params_Ni_57NI = get_parameters(stack, 'Ni', '57NI', False, wa)

    irradiation_time = 3600; unc_irradiation_time = 1  # s

    number_of_foils = 7; number_of_reactions = 4
    matrix_A0 = np.zeros((number_of_foils,number_of_reactions))
    matrix_sigma_A0 = np.zeros((number_of_foils,number_of_reactions))
    matrix_lambda_ = np.zeros((number_of_foils,number_of_reactions))
    matrix_mass_density = np.zeros((number_of_foils,number_of_reactions))
    matrix_sigma_mass_density = np.zeros((number_of_foils,number_of_reactions))
    matrix_reaction_integral = np.zeros((number_of_foils,number_of_reactions))
    matrix_uncertainty_integral = np.zeros((number_of_foils,number_of_reactions))
    matrix_irr_time = np.zeros((number_of_foils,1))
    matrix_sigma_irr_time = np.zeros((number_of_foils,1))

    # list_of_params = [params_Cu_65Zn, params_Cu_63Zn, params_Cu_62Zn, params_Cu_58CO, params_Cu_56CO, params_Ni_57NI]
    list_of_params = [params_Cu_65Zn, params_Cu_63Zn, params_Cu_62Zn, params_Ni_57NI]
    # list_of_params = [params_Cu_65Zn, params_Cu_63Zn, params_Cu_62Zn, params_Cu_58CO, params_Ni_57NI]
    # list_of_params = np.array((list_of_params))

    list_of_params = np.array([ params_Cu_65Zn, params_Cu_63Zn, params_Cu_62Zn, params_Ni_57NI], dtype=object)
    # list_of_params = np.array([ params_Cu_65Zn, params_Cu_63Zn, params_Cu_62Zn, params_Cu_58CO, params_Cu_56CO, params_Ni_57NI], dtype=object)
    #print(list_of_params[0,2])
    #print(list_of_params.shape)
    #print(type(list_of_params))

    n = len(list_of_params)
    A0 = np.zeros(n)
    sigma_A0 = np.zeros(n)
    lambda_ = np.zeros(n)
    mass_density = np.zeros(n)
    sigma_mass_density =  np.zeros(n)
    reaction_integral = np.zeros(n)
    uncertainty_integral = np.zeros(n)
    matrix_irr_time = np.transpose(np.ones(number_of_foils)*3600)
    matrix_sigma_irr_time = np.transpose(np.ones(number_of_foils)*3)

    shape_cols = (10,)
    for i in range(len(list_of_params)):
        A0 = list_of_params[i,0]
        sigma_A0 = list_of_params[i,1]
        lambda_ = list_of_params[i,2]
        mass_density = list_of_params[i,3]
        sigma_mass_density = list_of_params[i,4]
        reaction_integral = list_of_params[i,5]
        uncertainty_integral = list_of_params[i,6]
        #irr_time = list_of_params[i,7]
        #sigma_irr_time = list_of_params[i,8]

        try:
            matrix_lambda_[:, i] = lambda_
            #matrix_irr_time[i] = irr_time
            #matrix_sigma_irr_time[i] = sigma_irr_time
            matrix_A0[:,i] = A0
            matrix_sigma_A0[:,i] = sigma_A0
            matrix_mass_density[:,i] = mass_density

            matrix_sigma_mass_density[:,i] = sigma_mass_density
            matrix_reaction_integral[:,i]  = reaction_integral
            matrix_uncertainty_integral[:,i]  = uncertainty_integral
        except:
            #print("Shape problem with Fe_56Co ")
            A0 = np.pad(A0, (0, 7), 'constant')
            sigma_A0 = np.pad(sigma_A0, (0, 7), 'constant')
            mass_density = np.pad(mass_density, (0, 7), 'constant')
            sigma_mass_density = np.pad(sigma_mass_density, (0, 7), 'constant')
            reaction_integral = np.pad(reaction_integral, (0, 7), 'constant')
            uncertainty_integral = np.pad(uncertainty_integral, (0, 7), 'constant')
            matrix_A0[:,i] = A0
            matrix_sigma_A0[:,i] = sigma_A0
            matrix_mass_density[:,i] = mass_density
            matrix_sigma_mass_density[:,i] = sigma_mass_density
            matrix_reaction_integral[:,i]  = reaction_integral
            matrix_uncertainty_integral[:,i]  = uncertainty_integral

        
    ### NEED TO REMOVE inf from sigma_A0 (which was caused by changing A0 for 56Ni to 0 in foil 4-10.
    ### Gave inf problems in zero division.  )
    rows = matrix_A0.shape[0]
    cols = matrix_A0.shape[1]
    for i in range(rows):
        for j in range(cols):
            if np.isinf(matrix_sigma_A0[i,j]):
                matrix_sigma_A0[i,j]=0
    #print(matrix_sigma_A0)
            #print(matrix_sigma_A0[i,j])

    # print(matrix_sigma_A0)
    # print(matrix_lambda_.shape)
    # print(matrix_irr_time)
    # print(matrix_sigma_irr_time)
    # print(matrix_A0)
    # print(matrix_sigma_A0)
    # print(matrix_mass_density)
    # # print(sigma_matrix_mass_density)
    # print(matrix_reaction_integral)
    # print(matrix_uncertainty_integral)


    return matrix_A0, matrix_sigma_A0, matrix_lambda_, matrix_mass_density, matrix_sigma_mass_density, matrix_reaction_integral, matrix_uncertainty_integral, matrix_irr_time, matrix_sigma_irr_time


# # Hannah:
# number_of_monitor_foils = 2
# monitor_reactions_per_foil = np.array([5, 1])

### Read in numbers of decays from csv file

def Average_BeamCurrent(A0, sigma_A0, mass_density, sigma_mass_density, lambda_, reaction_integral, uncertainty_integral, irr_time, sigma_irr_time, csv_filename='averaged_currents.csv', save_csv=False):
    number_of_monitor_foils = 2
    monitor_reactions_per_foil = np.array([3, 1])
    # monitor_reactions_per_foil = np.array([4, 1])
    # monitor_reactions_per_foil = np.array([5, 1])
    def decomment(csvfile):
        for row in csvfile:
            raw = row.split('#')[0].strip()
            if raw: yield raw

    def read_csv(name_of_csv_file):
        results = []
        with open(name_of_csv_file) as csvfile:
            reader = csv.reader(decomment(csvfile))
            for row in reader:
                results.append(row)

        return np.asarray(results, dtype=float)

    def beam_current(A0, rho_dr, lambdas, t_irradiation, reaction_integral):
        return (elementary_charge*1e9 * A0) / (rho_dr * (1 - np.exp(-lambdas*t_irradiation)) * reaction_integral)

	# Numerical partial derivatives
    def dIdA0(A0, rho_dr, lambdas, t_irradiation, reaction_integral):
        delta_x = 1E-8 * A0
        return ((beam_current(A0 + (delta_x/2), rho_dr, lambdas, t_irradiation, reaction_integral) - beam_current(A0 - (delta_x/2), rho_dr, lambdas, t_irradiation, reaction_integral)) / delta_x)
    def dIdRhoDr(A0, rho_dr, lambdas, t_irradiation, reaction_integral):
        delta_x = 1E-8 * rho_dr
        return ((beam_current(A0, rho_dr + (delta_x/2), lambdas, t_irradiation, reaction_integral) - beam_current(A0, rho_dr - (delta_x/2), lambdas, t_irradiation, reaction_integral)) / delta_x)
    def dIdLambda(A0, rho_dr, lambdas, t_irradiation, reaction_integral):
        delta_x = 1E-8 * lambdas
        return ((beam_current(A0, rho_dr, lambdas + (delta_x/2), t_irradiation, reaction_integral) - beam_current(A0, rho_dr, lambdas - (delta_x/2), t_irradiation, reaction_integral)) / delta_x)
    def dIdTIrradiation(A0, rho_dr, lambdas, t_irradiation, reaction_integral):
        delta_x = 1E-8 * t_irradiation
        return ((beam_current(A0, rho_dr, lambdas, t_irradiation + (delta_x/2), reaction_integral) - beam_current(A0, rho_dr, lambdas, t_irradiation - (delta_x/2), reaction_integral)) / delta_x)
    def dIdIntegral(A0, rho_dr, lambdas, t_irradiation, reaction_integral):
        delta_x = 1E-8 * reaction_integral
        return ((beam_current(A0, rho_dr, lambdas, t_irradiation, reaction_integral + (delta_x/2)) - beam_current(A0, rho_dr, lambdas, t_irradiation, reaction_integral - (delta_x/2))) / delta_x)

	# Approximate uncertainties in beam current
    def sigma_I_approximate(A0, rho_dr, lambdas, t_irradiation, reaction_integral, unc_A0, unc_rho_dr, unc_lambdas, unc_t_irradiation, unc_reaction_integral):
        approx_error = np.sqrt(np.power(dIdA0(A0, rho_dr, lambdas, t_irradiation, reaction_integral) * unc_A0,2) +
        np.power(dIdRhoDr(A0, rho_dr, lambdas, t_irradiation, reaction_integral) * unc_rho_dr,2) +
        np.power(dIdLambda(A0, rho_dr, lambdas, t_irradiation, reaction_integral) * unc_lambdas,2) +
        np.power(dIdTIrradiation(A0, rho_dr, lambdas, t_irradiation, reaction_integral) * unc_t_irradiation,2) +
        #np.power( unc_reaction_integral,2))
        np.power(dIdIntegral(A0, rho_dr, lambdas, t_irradiation, reaction_integral) * unc_reaction_integral,2))
		# approx_error = 0
		# approx_error = np.power(dIdA0(A0, rho_dr, lambdas, t_irradiation, reaction_integral),2)
        return approx_error

    number_of_monitor_reactions = 0
    #print('yo')
    submatrix_lower_indices = np.zeros(number_of_monitor_foils)
    submatrix_upper_indices = np.zeros(number_of_monitor_foils)

	#### PARAMETERS IN THE FUNCTION des19_BeamCurrent.py
	#A0, dA0, mass_density, sigma_mass_density, lambda_, reaction_integral, uncertainty_integral, irr_time, sigma_irr_time

	# Load in activation data
	#
	# All in nuclei / cm^2
	#                       Cu       Ti
    areal_density = mass_density
    uncertainty_areal_density = sigma_mass_density

	# All in Bq
	#                      Sc46    V48    Zn62    Zn63
    EoB_activities = A0
    uncertainty_EoB_activities = sigma_A0
	# Normalized integral(sigma * dPhidE * dE)
	#
    reaction_integral = reaction_integral
    #unc_rxn_integral = uncertainty_integral #rxn_int * 	percent_rn_uncertainties      #rxn = reactions
    unc_rxn_integral = uncertainty_integral * reaction_integral#rxn_int * 	percent_rn_uncertainties      #rxn = reactions
	# All in 1/s
    lambdas = lambda_
	# print(type(lambdas))
	# All in s
    t_irradiation = irr_time
    uncertainty_t_irradiation = sigma_irr_time


    for i in range(0, number_of_monitor_foils):
        submatrix_lower_indices[i] = number_of_monitor_reactions
        number_of_monitor_reactions += monitor_reactions_per_foil[i]
        submatrix_upper_indices[i] = number_of_monitor_reactions

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
    output_foil_index = []
    output_mu = []
    output_unc_mu = []
    output_percent_unc = []


	# Get correlation submtarix for each monitor foil - n x n, where n= # of reactions per foil
    for i in range(0, number_of_monitor_foils): # Monitor reactions per foil [3,3,1]
        submatrix = np.ones((monitor_reactions_per_foil[i], monitor_reactions_per_foil[i]))
        corr_areal_density[submatrix_lower_indices[i]:submatrix_upper_indices[i], submatrix_lower_indices[i]:submatrix_upper_indices[i]] = submatrix
        corr_reaction_integral[submatrix_lower_indices[i]:submatrix_upper_indices[i], submatrix_lower_indices[i]:submatrix_upper_indices[i]] = 0.3*submatrix


	# Ensure all diagonal elements are still ones
    np.fill_diagonal(corr_lambda,1)
    np.fill_diagonal(corr_areal_density,1)
    np.fill_diagonal(corr_reaction_integral,1)
    np.fill_diagonal(corr_t_irradiation,1)
    np.fill_diagonal(corr_EoB_activities,1)


	# Loop over all beam positions
    number_of_energies = len(areal_density)
	# print(number_of_energies)
	# Test mode!!!!
	# number_of_energies = 1

	# Hold curents as we go along...
    currents = np.zeros((number_of_energies,number_of_monitor_reactions))
    unc_currents = np.zeros((number_of_energies,number_of_monitor_reactions))
    function_dictionary = {'0':dIdA0, '1':dIdRhoDr, '2':dIdLambda, '3':dIdTIrradiation, '4':dIdIntegral}

    for i_energy in range(0, number_of_energies):
        #print('i_energy: ',i_energy)
		# Get nonzero entries in A0:
        nonzero_indices = np.nonzero(EoB_activities[i_energy,:])
        ad = areal_density[i_energy,:]
        unc_ad = uncertainty_areal_density[i_energy,:]
        A0 = EoB_activities[i_energy,:]
        unc_A0 = uncertainty_EoB_activities[i_energy,:]
        rxn_int = reaction_integral[i_energy,:]
        delta_t = np.ones(number_of_monitor_reactions) *t_irradiation[i_energy]
        unc_delta_t = np.ones(number_of_monitor_reactions) *uncertainty_t_irradiation[i_energy]
        unc_rxn_int = unc_rxn_integral[i_energy,:]
        loop_lambdas = lambdas
        uncertainty_lambdas = loop_lambdas * 0.001

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
            # print('disjoint indices: ', disjoint_indices)

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

        # print('beam_current inputs: ', A0, ad, loop_lambdas[i_energy,:], delta_t, rxn_int)
        temp_currents =  beam_current(A0, ad, loop_lambdas[i_energy,:], delta_t, rxn_int)
        currents[i_energy, :] =  temp_currents
        # print('temp_currents: ', temp_currents)
        # print('unc_beam_current inputs: ', A0, ad, loop_lambdas[i_energy,:], delta_t, rxn_int, unc_A0, unc_ad, uncertainty_lambdas[i_energy,:], unc_delta_t, unc_rxn_int)
        unc_temp_currents = sigma_I_approximate(A0, ad, loop_lambdas[i_energy,:], delta_t, rxn_int, unc_A0, unc_ad, uncertainty_lambdas[i_energy,:], unc_delta_t, unc_rxn_int)
        unc_currents[i_energy,:] = unc_temp_currents

        value_array = np.array([A0, ad, loop_lambdas[0], delta_t, rxn_int])
        uncertainty_array = np.array([unc_A0, unc_ad, uncertainty_lambdas[0], unc_delta_t, unc_rxn_int])
        correlation_array = np.array([loop_corr_EoB_activities, loop_corr_areal_density, loop_corr_lambda, loop_corr_t_irradiation, loop_corr_reaction_integral])


		# Set up covariance matrix for current energy position
        cov = np.zeros((len(nonzero_indices[0]),len(nonzero_indices[0])))

		# NaN handling - replace range(0,number_of_monitor_reactions) with indices of nonzero elements of A0?
		# Fill correlation matrices
        for i_index,i_element in enumerate(nonzero_indices[0]):
            for j_index,j_element in enumerate(nonzero_indices[0]):
                for dict_index,dict_key in enumerate(function_dictionary):
                    dIdxi = function_dictionary[dict_key](A0[i_element], ad[i_element], loop_lambdas[0,i_element], delta_t[i_element], rxn_int[i_element])
                    dIdxj = function_dictionary[dict_key](A0[j_element], ad[j_element], loop_lambdas[0,j_element], delta_t[j_element], rxn_int[j_element])
                    cov[i_index,j_index] += dIdxi * uncertainty_array[int(dict_key),i_element] *  correlation_array[int(dict_key),i_index,j_index] *  uncertainty_array[int(dict_key),j_element] * dIdxj

		# print("Final covariance matrix: \n", cov)
        inverted_covariance = np.linalg.inv(cov)
        numerator = 0.0
        denominator = 0.0

        for i_index,i_element in enumerate(nonzero_indices[0]):
            for j_index,j_element in enumerate(nonzero_indices[0]):
                numerator += temp_currents[j_element] * inverted_covariance[i_index,j_index]
                denominator +=  inverted_covariance[i_index,j_index]

        weighted_average_current = numerator/denominator
        uncertainty_weighted_average_current = np.sqrt(1.0/denominator)


        ####print("weighted_average_current: ",weighted_average_current, " +/- ",uncertainty_weighted_average_current, " nA     (", 100*uncertainty_weighted_average_current/weighted_average_current ," %)")

		# Append values for current energy
        output_foil_index.append(i_energy)
        output_mu.append(weighted_average_current)
        output_unc_mu.append(uncertainty_weighted_average_current)
        output_percent_unc.append(100*uncertainty_weighted_average_current/weighted_average_current)

	# Save final results to csv
    outfile = np.stack((np.transpose(output_foil_index),np.transpose(output_mu),np.transpose(output_unc_mu),np.transpose(output_percent_unc)), axis=-1)
    df = pd.DataFrame(outfile,columns=["foil_index", "Average Current (nA)", "Uncertainty in Average Current (nA)", "% Uncertainty"])
    print(df)
    # csv_outname = 'WABC_' + csv_filename[10:-11] + '.csv'
    path = os.getcwd() + '/generatedfiles/beamcurrent/true_weighted_average_beam_current/'
    csv_outname = path + csv_filename 
    print(csv_outname)
    
    
    # if save_csv==True:
        # np.savetxt("./{}".format(csv_outname), outfile, delimiter=",", header="Foil Index, Average Current (nA), Uncertainty in Average Current (nA), % Uncertainty")
        # np.savetxt("./{}".format(csv_outname), outfile, delimiter=",", header="Foil Index, Average Current (nA), Uncertainty in Average Current (nA), % Uncertainty")
    output_foil_index2 = np.array(output_foil_index) - 0.2

    return output_mu[::-1], output_unc_mu[::-1 ] #returning reversed lists


# matrix_A0, matrix_sigma_A0, matrix_lambda_, matrix_mass_density, matrix_sigma_mass_density, matrix_reaction_integral, matrix_uncertainty_integral, matrix_irr_time, matrix_sigma_irr_time = reshaping_parameters('stack_55_MeV', wa_55_after)
matrix_A0, matrix_sigma_A0, matrix_lambda_, matrix_mass_density, matrix_sigma_mass_density, matrix_reaction_integral, matrix_uncertainty_integral, matrix_irr_time, matrix_sigma_irr_time = reshaping_parameters('stack_30_MeV', wa_30_after)
# matrix_A0, matrix_sigma_A0, matrix_lambda_, matrix_mass_density, matrix_sigma_mass_density, matrix_reaction_integral, matrix_uncertainty_integral, matrix_irr_time, matrix_sigma_irr_time = reshaping_parameters('stack_55_MeV', wa_55_after)
I, dI = Average_BeamCurrent(matrix_A0, matrix_sigma_A0, matrix_mass_density, matrix_sigma_mass_density, matrix_lambda_, matrix_reaction_integral, matrix_uncertainty_integral, matrix_irr_time, matrix_sigma_irr_time, csv_filename='averaged_currents_55_MeV_stack.csv', save_csv=True)
I, dI = Average_BeamCurrent(matrix_A0, matrix_sigma_A0, matrix_mass_density, matrix_sigma_mass_density, matrix_lambda_, matrix_reaction_integral, matrix_uncertainty_integral, matrix_irr_time, matrix_sigma_irr_time, csv_filename='averaged_currents_30_MeV_stack.csv', save_csv=True)
print(I)
[1249, 365, 369, 177, 178, 180, 597, 469, 470, 472, 1242, 604]
[1047, 1048, 543, 544, 1064, 1067, 1076, 1079, 568, 57, 58, 1083, 1081, 1086, 576, 578, 71, 73, 588, 79, 81, 595, 83, 85, 88, 601, 90, 92, 94, 1119, 1120, 96, 97, 100, 101, 102, 103, 104, 108, 109, 1135, 113, 1138, 117, 119, 121, 1146, 124, 1149, 1151, 1153, 132, 1157, 1188, 1189, 690, 691, 1213, 704, 706, 1221, 1223, 712, 715, 717, 719, 1233, 722, 725, 727, 1240, 729, 731, 732, 1246, 735, 736, 737, 738, 739, 743, 744, 748, 752, 754, 756, 759, 250, 251, 769, 261, 264, 271, 274, 276, 278, 281, 311, 312, 335, 850, 851, 343, 346, 864, 866, 356, 872, 362, 875, 363, 877, 366, 879, 882, 886, 889, 891, 894, 895, 898, 899, 900, 901, 905, 907, 911]