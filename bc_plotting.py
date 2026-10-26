from beamcurrent import *
from flux_stack import *
from weighted_average_updated import *

################ BEAM CURRENT ##########################


mon_30MeV = [
    ('Cu', '63ZN'),
    ('Cu', '65ZN'),
    ('Cu', '62ZN'),
    ('Cu', '61CU'), #testing
    ('Cu', '58CO'), # only 2-3 last points 
    ('Cu', '56CO'), # not valid for 30 mev
    ('Ni', '56NI'),
    ('Ni', '57NI'), # too high testing
    ('Ni', '55CO'),
    ('Ni', '56CO'),
    ('Ni', '57CO'),
    ('Ni', '58CO'),
    ]

mon_55MeV = [
    ('Cu', '63ZN'),
    ('Cu', '65ZN'),
    ('Cu', '62ZN'),
    ('Cu', '61CU'),
    ('Cu', '58CO'),
    ('Cu', '56CO'), # as expected, too close to threshold -> too high values
    ('Ni', '56NI'),
    ('Ni', '57NI'),
    ('Ni', '55CO'),
    ('Ni', '56CO'),
    ('Ni', '57CO'), # wayyy too low.... 
    ('Ni', '58CO'), # too high..... 
    ]
flux_30_after = 'TaSn_stack_30MeV_dp_1.010%_fluxes.csv' ; full_stack_30 = 'TaSn_stack_30MeV_dp_1.010%.csv'
flux_30_after2 = 'TaSn_stack_30MeV_dp_1.020%_fluxes.csv' ; full_stack_30 = 'TaSn_stack_30MeV_dp_1.010%.csv'
flux_55_after = 'TaSn_stack_55MeV_dp_1.020%_fluxes.csv' ; full_stack_55 = 'TaSn_stack_55MeV_dp_1.010%.csv'
flux_55_after2 = 'TaSn_stack_55MeV_dp_1.010%_fluxes.csv' ; full_stack_55 = 'TaSn_stack_55MeV_dp_1.010%.csv'
flux_30_before = 'TaSn_stack_30MeV_dp_1.000%_fluxes.csv'; full_stack_30 = 'TaSn_stack_30MeV_dp_1.000%.csv'
flux_55_before = 'TaSn_stack_55MeV_dp_1.000%_fluxes.csv'; full_stack_55 = 'TaSn_stack_55MeV_dp_1.000%.csv'

wa_55_before = WeightedAverageFlux(stack_fluxes = flux_55_before, full_stack = full_stack_55, stack = 'stack_55_MeV')
wa_55_after = WeightedAverageFlux(stack_fluxes = flux_55_after, full_stack = full_stack_55, stack = 'stack_55_MeV')
wa_55_after2 = WeightedAverageFlux(stack_fluxes = flux_55_after2, full_stack = full_stack_55, stack = 'stack_55_MeV')
wa_30_before = WeightedAverageFlux(stack_fluxes = flux_30_before, full_stack = full_stack_30, stack = 'stack_30_MeV')
wa_30_after = WeightedAverageFlux(stack_fluxes = flux_30_after, full_stack = full_stack_30, stack = 'stack_30_MeV')
wa_30_after2 = WeightedAverageFlux(stack_fluxes = flux_30_after2, full_stack = full_stack_30, stack = 'stack_30_MeV')

bc_30_before = BeamCurrent(wa_30_before, stack = 'stack_30_MeV', monitor_reactions=mon_30MeV)
bc_30_after = BeamCurrent(wa_30_after, stack = 'stack_30_MeV', monitor_reactions=mon_30MeV)
bc_30_after2 = BeamCurrent(wa_30_after2, stack = 'stack_30_MeV', monitor_reactions=mon_30MeV)
bc_55_before = BeamCurrent(wa_55_before, stack = 'stack_55_MeV', monitor_reactions=mon_55MeV)
bc_55_after = BeamCurrent(wa_55_after, stack = 'stack_55_MeV', monitor_reactions=mon_55MeV)
bc_55_after2 = BeamCurrent(wa_55_after2, stack = 'stack_55_MeV', monitor_reactions=mon_55MeV)


params_all = [
    'params_Cu_65Zn', 'params_Cu_63Zn', 'params_Cu_62Zn', 'params_Cu_61Cu', 'params_Cu_58Co', 'params_Cu_56Co', 
    'params_Ni_57Ni', 'params_Ni_56Ni', 'params_Ni_58Co', 'params_Ni_57Co', 'params_Ni_56Co', 'params_Ni_55Co'
    ]

params_55_stack = [
    'params_Cu_65Zn', 'params_Cu_63Zn', 'params_Cu_62Zn', 'params_Cu_61Cu', 'params_Cu_58Co',
    'params_Ni_57Ni', 'params_Ni_56Ni', 'params_Ni_58Co', 'params_Ni_57Co', 'params_Ni_56Co', 'params_Ni_55Co'
]

params_30_stack = [
    'params_Cu_65Zn', 'params_Cu_63Zn', 'params_Cu_62Zn', 'params_Cu_61Cu',
    'params_Ni_56Ni', 'params_Ni_58Co', 'params_Ni_57Co', 'params_Ni_56Co', 'params_Ni_55Co'
]

params_old_monitors = [
    'params_Cu_65Zn', 'params_Cu_63Zn', 'params_Cu_62Zn', 'params_Cu_58Co', 'params_Cu_56Co', 'params_Ni_57Ni', 
]
params_old_monitors_30 = [
    'params_Cu_65Zn', 'params_Cu_63Zn', 'params_Cu_62Zn', 'params_Cu_58Co','params_Ni_57Ni', 
]

# wabc_30_after = WeightedAverageBeamCurrent(7, 2, np.array([3,5]), wa_30_after, 'stack_30_MeV', params_30_stack )#.weighted_average_beam_current('new_params_30')
# wabc_30_after1 = WeightedAverageBeamCurrent(7, 2, np.array([6,6]), wa_30_after, 'stack_30_MeV', params_all)#.weighted_average_beam_current('new_all_params_30')
# wabc_30_after2 = WeightedAverageBeamCurrent(7, 2, np.array([5, 1]), wa_30_after, 'stack_30_MeV', params_old_monitors)#.weighted_average_beam_current('new_old_params_30')

# wabc_30_old_dp1_01_iaea = WeightedAverageBeamCurrent(7, 2, np.array([5, 1]), WeightedAverageFlux(stack_fluxes = flux_30_after, full_stack = full_stack_30, stack = 'stack_30_MeV', iaea_version='2017'), 'stack_30_MeV', params_old_monitors)#.weighted_average_beam_current('old_old_params_30')
# wabc_55_old_dp1_01_iaea = WeightedAverageBeamCurrent(7, 2, np.array([5, 1]), WeightedAverageFlux(stack_fluxes = flux_55_after2, full_stack = full_stack_30, stack = 'stack_55_MeV', iaea_version='2017'), 'stack_55_MeV', params_old_monitors)#.weighted_average_beam_current('old_old_params_55')

wabc_55_after = WeightedAverageBeamCurrent(7, 2, np.array([5,6]), wa_55_before, 'stack_55_MeV', params_55_stack)#.weighted_average_beam_current('new_params_55')
# wabc_55_after1 = WeightedAverageBeamCurrent(7, 2, np.array([6,6]), wa_55_after, 'stack_55_MeV', params_all)#.weighted_average_beam_current('new_all_params_55')
# wabc_55_after2 = WeightedAverageBeamCurrent(7, 2, np.array([5,1]), wa_55_after, 'stack_55_MeV', params_old_monitors)#.weighted_average_beam_current('new_old_params_55')

# wabc_30_after.plot_beam_current_isotope('Cu', 'params_Cu_65Zn', color='peru', label='Cu 65Zn')
# wabc_30_after.plot_beam_current_isotope('Cu', 'params_Cu_63Zn', color='blue', label='Cu 63Zn')
# wabc_30_after.plot_beam_current_isotope('Cu', 'params_Cu_62Zn', color='yellow', label='Cu 62Zn')
# wabc_30_after.plot_beam_current_isotope('Cu', 'params_Cu_61Cu', color='darkred', label='Cu 61Cu')
# wabc_30_after.plot_beam_current_isotope('Cu', 'params_Cu_56Co', color='#a64d79', label='Cu 56Co')
# wabc_30_after.plot_beam_current_isotope('Cu', 'params_Cu_58Co', color='#6aa84f', label='Cu 58Co')
# wabc_30_after.plot_beam_current_isotope('Ni', 'params_Ni_57Ni', color='#3d85c6', label='Ni 57Ni')
# wabc_30_after.plot_beam_current_isotope('Ni', 'params_Ni_56Ni', color='orange', label='Ni 56Ni')
# wabc_30_after.plot_beam_current_isotope('Ni', 'params_Ni_58Co', color='green', label='Ni 58Co')
# wabc_30_after.plot_beam_current_isotope('Ni', 'params_Ni_57Co', color='orchid', label='Ni 57Co')
# wabc_30_after.plot_beam_current_isotope('Ni', 'params_Ni_56Co', color='hotpink', label='Ni 56Co')
# wabc_30_after.plot_beam_current_isotope('Ni', 'params_Ni_55Co', color='royalblue', label='Ni 55Co')
# wabc_30_after.plot_weighted_average('darkred', 'params 30 stack')
# wabc_30_after1.plot_weighted_average('orange', 'params all')
# wabc_30_old_dp1_01_iaea.plot_weighted_average('#a64d79', 'params old moni')
# plt.legend()
# plt.show()

wabc_55_after.plot_beam_current_isotope('Cu', 'params_Cu_65Zn', color='peru', label='Cu 65Zn')
wabc_55_after.plot_beam_current_isotope('Cu', 'params_Cu_63Zn', color='blue', label='Cu 63Zn')
wabc_55_after.plot_beam_current_isotope('Cu', 'params_Cu_62Zn', color='yellow', label='Cu 62Zn')
wabc_55_after.plot_beam_current_isotope('Cu', 'params_Cu_61Cu', color='darkred', label='Cu 61Cu')
wabc_55_after.plot_beam_current_isotope('Cu', 'params_Cu_56Co', color='#a64d79', label='Cu 56Co')
wabc_55_after.plot_beam_current_isotope('Cu', 'params_Cu_58Co', color='#6aa84f', label='Cu 58Co')
wabc_55_after.plot_beam_current_isotope('Ni', 'params_Ni_57Ni', color='#3d85c6', label='Ni 57Ni')
wabc_55_after.plot_beam_current_isotope('Ni', 'params_Ni_56Ni', color='orange', label='Ni 56Ni')
wabc_55_after.plot_beam_current_isotope('Ni', 'params_Ni_58Co', color='green', label='Ni 58Co')
wabc_55_after.plot_beam_current_isotope('Ni', 'params_Ni_57Co', color='orchid', label='Ni 57Co')
wabc_55_after.plot_beam_current_isotope('Ni', 'params_Ni_56Co', color='hotpink', label='Ni 56Co')
wabc_55_after.plot_beam_current_isotope('Ni', 'params_Ni_55Co', color='royalblue', label='Ni 55Co')
# wabc_55_after1.plot_weighted_average('orange', 'params all')
# wabc_55_old_dp1_01_iaea.plot_weighted_average('#a64d79', 'params old moni')
wabc_55_after.plot_weighted_average('darkred', 'params 55 stack')
plt.legend()
plt.show()



# wabc_55_after.average_beam_current()
# wabc_30_after.average_beam_current()

# save_bc = os.getcwd() + '/generatedfiles/beamcurrent/figures/' 

# bc_30_before.plot_all('30 MeV before variance minimization dp 1.00')
# plt.ylim(0,250)
# plt.legend(fontsize=5) 
# plt.savefig(save_bc + '30MeV_after_dp1.00.pdf')
# plt.show()

# bc_30_after.plot_all('30 MeV after variance minimization dp 1.01')
# wabc_30_after.plot_weighted_average()
# plt.ylim(0,250)
# plt.legend(fontsize=5) 
# # plt.savefig(save_bc + '30MeV_after_dp1.01.pdf')
# plt.show()

# bc_30_after2.plot_all('30 MeV after variance minimization dp 1.02')
# plt.ylim(0,250)
# plt.legend(fontsize=5) 
# plt.savefig(save_bc + '30MeV_after_dp1.02.pdf')
# plt.show()


# bc_55_before.plot_all('55 MeV before variance minimization dp 1.00')
# plt.ylim(0,220)
# plt.legend(fontsize=5) 
# plt.savefig(save_bc + '55MeV_after_dp1.00.pdf')
# plt.show()


# bc_55_after.plot_all('55 MeV after variance minimization dp 1.02')
# plt.ylim(0,210)
# plt.legend(fontsize=5) 
# # plt.savefig(save_bc + '55MeV_after_dp1.02.pdf')
# plt.show()

# bc_55_after2.plot_all('55 MeV after variance minimization dp 1.01')
# wabc_55_after.plot_weighted_average()
# plt.ylim(0,210)
# plt.legend(fontsize=5) 
# # plt.savefig(save_bc + '55MeV_after_dp1.01.pdf')
# plt.show()


#####################################################################################################
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)


def plot_weighted_average(dp, stack, note, color, params, number_of_compartments, number_of_monitor_foils, monitor_reactions_per_foil, csv_filename):
    # np.array([5,6])
    if stack == 'stack_55_MeV':
        flux = 'TaSn_stack_55MeV_dp_' + dp + '%_fluxes.csv'
    else:
        flux = 'TaSn_stack_30MeV_dp_' + dp + '%_fluxes.csv'
    full_stack = ''
    # flux = 'TaSn_stack_30MeV_dp_1.010%_fluxes.csv' ; full_stack = 'TaSn_stack_30MeV_dp_1.010%.csv'
    wa = WeightedAverageFlux(stack_fluxes = flux, full_stack = full_stack, stack = stack)
    # energy, flux = wa.get_flux_energy_stack(element='Ta')
    # print(flux)
    wabc = WeightedAverageBeamCurrent(number_of_compartments, number_of_monitor_foils, monitor_reactions_per_foil, wa, stack, params)
    # wabc_old = WeightedAverageBeamCurrent(number_of_compartments, number_of_monitor_foils, monitor_reactions_per_foil, wa_old, stack, params_old_monitors)
    wabc.weighted_average_beam_current(csv_filename)
    # wabc_old.weighted_average_beam_current('old')
    plot_params(params, wabc, color)
    wabc.plot_weighted_average(color, note)

    # wabc_old.plot_weighted_average(color='orange', note='iaea 17')


def plot_iaea_17(dp, stack, note, color, params, number_of_compartments, number_of_monitor_foils, monitor_reactions_per_foil, csv_filename):
    if stack == 'stack_55_MeV':
        flux = 'TaSn_stack_55MeV_dp_' + dp + '%_fluxes.csv'
    else:
        flux = 'TaSn_stack_30MeV_dp_' + dp + '%_fluxes.csv'
    full_stack = ''
    wa_old = WeightedAverageFlux(stack_fluxes = flux, full_stack = full_stack, stack = stack, iaea_version='2017')
    wabc_old = WeightedAverageBeamCurrent(number_of_compartments, number_of_monitor_foils, monitor_reactions_per_foil, wa_old, stack, params_old_monitors)
    wabc_old.weighted_average_beam_current(csv_filename)
    wabc_old.plot_weighted_average(color, note)

def plot_params(params, wabc, color):
    for p in params:
        if p == 'params_Cu_65Zn':
            # wabc.plot_beam_current_isotope('Cu', 'params_Cu_65Zn', color=color, label='Cu 65Zn')
            wabc.plot_beam_current_isotope('Cu', 'params_Cu_65Zn', color='#0072B2', label='Cu 65Zn', symbol='v')
        if p == 'params_Cu_63Zn':
            # wabc.plot_beam_current_isotope('Cu', 'params_Cu_63Zn', color=color, label='Cu 63Zn')
            wabc.plot_beam_current_isotope('Cu', 'params_Cu_63Zn', color= '#D55E00', label='Cu 63Zn', symbol='*')
        if p == 'params_Cu_62Zn':
            # wabc.plot_beam_current_isotope('Cu', 'params_Cu_62Zn', color=color, label='Cu 62Zn')
            wabc.plot_beam_current_isotope('Cu', 'params_Cu_62Zn', color='#009E73', label='Cu 62Zn', symbol='s')
        if p == 'params_Cu_61Cu':
            wabc.plot_beam_current_isotope('Cu', 'params_Cu_61Cu', color= '#CC79A7', label='Cu 61Cu', symbol='x')
            # wabc.plot_beam_current_isotope('Cu', 'params_Cu_61Cu', color=color, label='Cu 61Cu')
        # if p == 'params_Cu_56Co':
        #     # wabc.plot_beam_current_isotope('Cu', 'params_Cu_56Co', color=color, label='Cu 56Co')
        #     wabc.plot_beam_current_isotope('Cu', 'params_Cu_56Co', color='#a64d79', label='Cu 56Co')
        if p == 'params_Cu_58Co':
            # wabc.plot_beam_current_isotope('Cu', 'params_Cu_58Co', color=color, label='Cu 58Co')
            wabc.plot_beam_current_isotope('Cu', 'params_Cu_58Co', color='#E69F00', label='Cu 58Co', symbol='X')
        if p == 'params_Ni_57Ni':
            # wabc.plot_beam_current_isotope('Ni', 'params_Ni_57Ni', color=color, label='Ni 57Ni')
            wabc.plot_beam_current_isotope('Ni', 'params_Ni_57Ni', color='#56B4E9', label='Ni 57Ni', symbol='P')
        if p == 'params_Ni_56Ni':
            # wabc.plot_beam_current_isotope('Ni', 'params_Ni_56Ni', color=color, label='Ni 56Ni')
            wabc.plot_beam_current_isotope('Ni', 'params_Ni_56Ni', color= '#8B4513', label='Ni 56Ni', symbol='d')
        if p == 'params_Ni_58Co':
            # wabc.plot_beam_current_isotope('Ni', 'params_Ni_58Co', color=color, label='Ni 58Co')
            wabc.plot_beam_current_isotope('Ni', 'params_Ni_58Co', color= '#7B61A8', label='Ni 58Co', symbol='D')
        if p == 'params_Ni_57Co':
            # wabc.plot_beam_current_isotope('Ni', 'params_Ni_57Co', color=color, label='Ni 57Co')
            wabc.plot_beam_current_isotope('Ni', 'params_Ni_57Co', color='#2200b1', label='Ni 57Co')
        if p == 'params_Ni_56Co':
            # wabc.plot_beam_current_isotope('Ni', 'params_Ni_56Co', color=color, label='Ni 56Co')
            wabc.plot_beam_current_isotope('Ni', 'params_Ni_56Co', color='#009E9A', label='Ni 56Co', symbol='o')
        if p == 'params_Ni_55Co':
            # wabc.plot_beam_current_isotope('Ni', 'params_Ni_55Co', color=color, label='Ni 55Co')
            wabc.plot_beam_current_isotope('Ni', 'params_Ni_55Co', color='#C44E52', label='Ni 55Co', symbol='^')


def save_figs_30MeV(dp_str):
    dp_0 = '1.000'
    for dp in dp_str:
        # plot_weighted_average(dp, 'stack_30_MeV', 'dp: ' + dp, 'darkred', params_all, 7, 2, [6,6], 'test.csv')
        plot_weighted_average(dp, 'stack_30_MeV', 'weighted average', '#7E0031', params_30_stack, 7, 2, [4,6], 'test_30.csv')
        # plot_iaea_17(dp, 'stack_30_MeV', 'iaea 17: dp: ' + dp, 'orange', params_old_monitors, 7, 2, [5,1], 'test_30_old_iaea.csv')
        # plot_weighted_average(dp_0, 'stack_30_MeV', 'iaea 25: dp: ' + dp_0, '#CCDAC2', params_30_stack, 7, 2, [6,6], 'test_30_og.csv')
        # plot_weighted_average(dp, 'stack_30_MeV', 'dp: ' + dp, 'darkred', params_30_stack, 7, 2, [3,5], 'test.csv')
        plt.legend(fontsize=6, loc='best', ncols=3)
        plt.title('30 MeV stack - ' + dp)
        plt.ylim(0,200)
        plt.savefig('generatedfiles/varianceminimization/iaea25/beam_current_figs_new/allcurrents_stack_30_MeV_' + dp + '.pdf')
        # plt.show()
        plt.close()

def save_figs_50MeV(dp_str):
    dp_0 = '1.000'
    for dp in dp_str:
        plot_weighted_average(dp, 'stack_55_MeV', 'weighted average', '#7E0031', params_all, 7, 2, [6,6], 'test_55.csv')
        # plot_iaea_17(dp, 'stack_55_MeV', 'dp: ' + 'iaea 17: dp: ' + dp, 'orange', params_old_monitors, 7, 2, [5,1], 'test_55_old_iaea.csv')
        # plot_iaea_17(dp, 'stack_30_MeV', 'iaea 17: dp: ' + dp, 'orange', params_old_monitors, 7, 2, [5,1], 'test_30_old_iaea.csv')
        # plot_weighted_average(dp_0, 'stack_55_MeV', 'iaea 25: dp: ' + dp_0, '#CCDAC2', params_all, 7, 2, [6,6], 'test_55_og.csv')
        # plot_weighted_average(dp, 'stack_30_MeV', 'dp: ' + dp, 'darkred', params_30_stack, 7, 2, [3,5], 'test.csv')
        plt.title('50 MeV stack - ' + dp)
        plt.ylim(0,200)
        plt.legend(fontsize=6, loc='best', ncols=3)
        plt.savefig('generatedfiles/varianceminimization/iaea25/beam_current_figs_new/allcurrents_stack_55_MeV_' + dp + '.pdf')
        # plt.show()
        plt.close()


def plot_only_wabc():
    dp = '1.000'
    plot_weighted_average(dp, 'stack_55_MeV', 'dp: ' + dp, '#7E0031', params_all, 7, 2, [6,6], 'test_55_second_plot.csv')
    dp = '1.016'
    plot_weighted_average(dp, 'stack_55_MeV', 'dp: ' + dp, '#CCDAC2', params_all, 7, 2, [6,6], 'test_55_second_plot.csv')

    plt.show()
    dp = '1.000'
    plot_weighted_average(dp, 'stack_30_MeV', 'dp: ' + dp, '#7E0031',  params_30_stack, 7, 2, [4,6], 'test_30_second_plot.csv')
    dp = '1.006'
    plot_weighted_average(dp, 'stack_30_MeV', 'dp: ' + dp, '#CCDAC2',  params_30_stack, 7, 2, [4,6], 'test_30_second_plot.csv')
    plt.show()


# % 1.006 og 1.019  for 30 MeV stack
# % 1.015 og 1.016  for 55 MeV stack

# plot_weighted_average('1.000', 'stack_30_MeV', 'dp 1.00', 'darkred', params_30_stack, 7, 2, [3,5], 'test.csv')
# plot_weighted_average('1.000', 'stack_55_MeV', 'dp 1.00', 'darkred', params_all, 7, 2, [6,6], 'test.csv')
# plot_weighted_average('1.010', 'stack_30_MeV', 'dp 1.01', 'orange', params_30_stack, 7, 2, [3,5], 'test.csv')
# plot_weighted_average('1.020', 'stack_30_MeV', 'dp 1.02', '#a64d79', params_30_stack, 7, 2, [3,5], 'test.csv')
# plot_weighted_average('1.030', 'stack_30_MeV', 'dp 1.03', '#a64d79', params_30_stack, 7, 2, [3,5], 'test.csv')
# plot_weighted_average('1.040', 'stack_30_MeV', 'dp 1.04', '#a64d79', params_30_stack, 7, 2, [3,5], 'test.csv')
# plt.legend()
# plt.title('30 MeV stack')
# plt.show()


plot_weighted_average('1.010', 'stack_55_MeV', 'dp 1.00', 'orange', params_all, 7, 2, [6,6], 'test.csv')
plt.legend()
plt.show()
# plot_weighted_average('1.020', 'stack_55_MeV', 'dp 1.02', '#a64d79', params_all, 7, 2, [6,6], 'test.csv')
# 1. Hent ut alle eksisterende "handles" (grafikk) og "labels" (tekst)


# plot_weighted_average('1.000', 'stack_30_MeV', 'dp 1.000', 'darkred', params_all, 7, 2, [6,6], 'test.csv')
# # plot_weighted_average('1.000', 'stack_30_MeV', 'dp 1.000', 'darkred', params_30_stack, 7, 2, [3,5], 'test.csv')
# # plot_weighted_average('1.000', 'stack_55_MeV', 'dp 1.000', 'darkred', params_all, 7, 2, [6,6], 'test.csv')
# handles, labels = plt.gca().get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# plt.legend(by_label.values(), by_label.keys(), ncol=3)
# plt.title('55 MeV stack')
# plt.show()

# plot_weighted_average('1.006', 'stack_30_MeV', 'dp 1.006', 'dodgerblue', params_30_stack, 7, 2, [3,5], 'test.csv')
# plot_weighted_average('1.016', 'stack_55_MeV', 'dp 1.016', 'darkred', params_all, 7, 2, [6,6], 'test.csv')
# handles, labels = plt.gca().get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# plt.legend(by_label.values(), by_label.keys(), ncol=3)
# plt.title('55 MeV stack')
# plt.show()







dp_array = np.array([0.990, 0.991, 1.003, 1.001, 1.021, 0.998, 1.013, 
                    1.024, 1.070, 0.996, 1.050, 1.014, 0.900, 1.012, 
                    1.110, 1.200, 1.007, 1.027, 1.040, 0.997, 1.008, 
                    1.016, 1.022, 1.011, 1.025, 0.999, 0.994, 1.005, 
                    1.009, 1.026, 1.080, 1.090, 0.970, 1.023, 1.029, 
                    1.002, 1.019, 1.020, 1.006, 1.010, 0.993, 1.100, 
                    0.950, 1.000, 1.015, 1.030, 0.992, 1.018, 0.800, 
                    1.017, 1.028, 1.060, 1.004, 0.980, 0.960, 0.995])
dp_str = []
for dp in dp_array:
    dp_str.append(f"{dp:.3f}")

# save_figs_30MeV(dp_str)
# save_figs_50MeV(dp_str)



# plot_only_wabc()


## 30 MeV 

# plot_weighted_average(dp='1.000', stack='stack_30_MeV', note='1.000', color='#7E0031', params=params_30_stack, number_of_compartments=7, number_of_monitor_foils=2, monitor_reactions_per_foil = [4,6], csv_filename=None)
# # plt.title('dp')
# plt.savefig('bc_Ni_57Ni.pdf')
# plt.show()


# plot_weighted_average(dp='1.006', stack='stack_30_MeV', note='1.006', color='#CCDAC2', params=params_30_stack, number_of_compartments=7, number_of_monitor_foils=2, monitor_reactions_per_foil = [4,6], csv_filename=None)
# plot_weighted_average(dp='1.000', stack='stack_30_MeV', note='weighted average beam current', color='#7E0031', params=params_30_stack, number_of_compartments=7, number_of_monitor_foils=2, monitor_reactions_per_foil = [4,6], csv_filename='33_no_varmin')


# plot_weighted_average(dp='1.006', stack='stack_30_MeV', note='weighted average beam current', color='#7E0031', params=params_30_stack, number_of_compartments=7, number_of_monitor_foils=2, monitor_reactions_per_foil = [4,6], csv_filename=None)
# plt.legend(fontsize=7, loc='best', ncols=3)
# plt.savefig('30MeV_1.006.pdf')
# plt.show()


# plot_weighted_average(dp='1.016', stack='stack_55_MeV', note='1.016', color='#CCDAC2', params=params_55_stack, number_of_compartments=7, number_of_monitor_foils=2, monitor_reactions_per_foil = [5,6], csv_filename=None)
# plt.legend(fontsize=7, loc='best', ncols=3)
# plt.show()



# plot_weighted_average(dp='1.000', stack='stack_30_MeV', note='1.000', color='#7E0031', params=params_all, number_of_compartments=7, number_of_monitor_foils=2, monitor_reactions_per_foil = [6,6], csv_filename=None)
# plot_weighted_average(dp='1.000', stack='stack_55_MeV', note='1.000', color='#7E0031', params=params_all, number_of_compartments=7, number_of_monitor_foils=2, monitor_reactions_per_foil = [6,6], csv_filename='55_no_varmin')
# plot_weighted_average(dp='1.016', stack='stack_55_MeV', note='1.016', color='#CCDAC2', params=params_all, number_of_compartments=7, number_of_monitor_foils=2, monitor_reactions_per_foil = [6,6], csv_filename=None)
# plot_weighted_average(dp='1.016', stack='stack_55_MeV', note='weighted average beam current', color='#7E0031', params=params_all, number_of_compartments=7, number_of_monitor_foils=2, monitor_reactions_per_foil = [6,6], csv_filename=None)
# plt.savefig('55MeV_1.016.pdf')
# plt.savefig('bc55_Cu_56Co.pdf')
# plt.title('dp 1.000')
# plt.show()
# plt.legend(fontsize=6, loc='best', ncols=3)
# plt.title('dp 1.016')
# plt.show()

# plot_iaea_17(dp='1.020', stack='stack_55_MeV', note='', color='#7E0031', params=params_old_monitors, number_of_compartments=7, number_of_monitor_foils=2, monitor_reactions_per_foil = [5,1], csv_filename='iaea17_55MeV_dp1.020')
# plot_iaea_17(dp='1.010', stack='stack_30_MeV', note='', color='#7E0031', params=params_old_monitors_30, number_of_compartments=7, number_of_monitor_foils=2, monitor_reactions_per_foil = [4,1], csv_filename='iaea17_30MeV_dp1.010')




no_30 = 'generatedfiles/beamcurrent/true_weighted_average_beam_current/33_no_varmin.csv'
yes_30 = 'generatedfiles/beamcurrent/true_weighted_average_beam_current/iaea25_30MeV_dp1.006.csv'

no_55 = 'generatedfiles/beamcurrent/true_weighted_average_beam_current/55_no_varmin.csv'
yes_55 = 'generatedfiles/beamcurrent/true_weighted_average_beam_current/iaea25_55MeV_dp1.016.csv'


def calculate_perc_increase(df_with, df_no):
    df1 = pd.read_csv(df_with)
    df2 = pd.read_csv(df_no)
    bc1 = df1['Average Current (nA)'].values
    bc2 = df2['Average Current (nA)'].values

    perc_change = (bc1 - bc2) / bc2 * 100
    for p in perc_change:
        if p >= 0:
            print(f"Økning på {p:.2f}%")
        else:
            print(f"Nedgang på {abs(p):.2f}%")
    return perc_change
        # perc_increase = (bc1[i] - bc2[i]) / bc2[i] * 100
        # print(str(perc_increase) + '%')

# print("******Stack 30*******")
# calculate_perc_increase(yes_30, no_30)
# print("******Stack 55*******")
# calculate_perc_increase(yes_55, no_55)