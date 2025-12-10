import numpy as np
import curie as ci
import pandas as pd
import os

# x_kapton = 0.013
# x_silicone = 0.013
# ad_degrader_a = 599.   #2.24 mm
# ad_degrader_b = 415.0   #1.55 mm
# ad_degrader_c = 261.5   #0.97 mm
# # ad_degrader_d = 599.0
# ad_degrader_e = 68.3   #0.256 mm
# ad_degrader_h = 33.8
# ad_be_backing = 4.425     #23.9130435 microns


df_cu = pd.read_csv('./generatedfiles/arealdensity/areal_density_Cu.csv')
df_ni = pd.read_csv('./generatedfiles/arealdensity/areal_density_Ni.csv')
df_ta = pd.read_csv('./generatedfiles/arealdensity/areal_density_Ta.csv')
df_sn = pd.read_csv('./generatedfiles/arealdensity/areal_density_Sn.csv')

ad_cu = df_cu['arealdensity g/cm2']*1e3; ad_ni = df_ni['arealdensity g/cm2']*1e3; ad_ta = df_ta['arealdensity g/cm2']*1e3; ad_sn = df_sn['arealdensity g/cm2']*1e3
save_stack_dir = './generatedfiles/stack/'



def variance_min(listOfPercentages):
    x_kapton = 0.013
    x_silicone = 0.013
    ad_degrader_a = 599.   #2.24 mm
    ad_degrader_b = 415.0   #1.55 mm
    ad_degrader_c = 261.5   #0.97 mm
    # ad_degrader_d = 599.0
    ad_degrader_e = 68.3   #0.256 mm
    ad_degrader_h = 33.8
    ad_be_backing = 4.425     #23.9130435 microns

    i = 0.075 # 5 percent decrease
    
    ad_degrader_a = i * ad_degrader_a
    ad_degrader_b = i * ad_degrader_b
    ad_degrader_e = i * ad_degrader_e
    return ad_degrader_a, ad_degrader_b, ad_degrader_c, ad_degrader_e

ad_degrader_a, ad_degrader_b, ad_degrader_c, ad_degrader_e = variance_min([])

def stack_55(E_p):
    stack = [
        #compartment01
       {'compound':'Ni', 'name':'Ni01', 'ad':ad_ni[0]},  # ad in mg/cm^2
       {'compound':'Sn', 'name':'Sn01', 'ad':ad_sn[0]},
       {'compound':'Ta', 'name':'Ta01', 'ad':ad_ta[0]},
       {'compound':'Cu', 'name':'Cu01', 'ad':ad_cu[0]},
       {'compound':'Al', 'name':'Al_degrader_A1', 'ad':ad_degrader_a},
       # compartment02
       {'compound':'Ni', 'name':'Ni02', 'ad':ad_ni[1]}, 
       {'compound':'Sn', 'name':'Sn02', 'ad':ad_sn[1]},
       {'compound':'Ta', 'name':'Ta02', 'ad':ad_ta[1]},
       {'compound':'Cu', 'name':'Cu02', 'ad':ad_cu[1]},
       {'compound':'Al', 'name':'Al_degrader_B2', 'ad':ad_degrader_b},
       # compartment03
       {'compound':'Ni', 'name':'Ni03', 'ad':ad_ni[2]},
       {'compound':'Sn', 'name':'Sn03', 'ad':ad_sn[2]},
       {'compound':'Ta', 'name':'Ta03', 'ad':ad_ta[2]},
       {'compound':'Cu', 'name':'Cu03', 'ad':ad_cu[2]},
       {'compound':'Al', 'name':'Al_degrader_C5', 'ad':ad_degrader_c},
       # compartment04
       {'compound':'Ni', 'name':'Ni04', 'ad':ad_ni[3]},
       {'compound':'Sn', 'name':'Sn04', 'ad':ad_sn[3]},
       {'compound':'Ta', 'name':'Ta04', 'ad':ad_ta[3]},
       {'compound':'Cu', 'name':'Cu04', 'ad':ad_cu[3]},
       {'compound':'Al', 'name':'Al_degrader_E8_E1', 'ad':2*ad_degrader_e},
       # compartment05
       {'compound':'Ni', 'name':'Ni05', 'ad':ad_ni[4]},
       {'compound':'Sn', 'name':'Sn05', 'ad':ad_sn[4]},
       {'compound':'Ta', 'name':'Ta05', 'ad':ad_ta[4]},
       {'compound':'Cu', 'name':'Cu05', 'ad':ad_cu[4]},
       {'compound':'Al', 'name':'Al_degrader_E4_E3_E2', 'ad':3*ad_degrader_e},
       # compartment06
       {'compound':'Ni', 'name':'Ni06', 'ad':ad_ni[5]},
       {'compound':'Sn', 'name':'Sn06', 'ad':ad_sn[5]},
       {'compound':'Ta', 'name':'Ta06', 'ad':ad_ta[5]},
       {'compound':'Cu', 'name':'Cu06', 'ad':ad_cu[5]},
       {'compound':'Al', 'name':'Al_degrader_E7_E6_E5', 'ad':3*ad_degrader_e},
       # compartment07
       {'compound':'Ni', 'name':'Ni07', 'ad':ad_ni[6]},
       {'compound':'Sn', 'name':'Sn07', 'ad':ad_sn[6]},
       {'compound':'Ta', 'name':'Ta07', 'ad':ad_ta[6]},
       {'compound':'Cu', 'name':'Cu07', 'ad':ad_cu[6]},
    ]
    full_stack = ci.Stack(stack, E0=E_p, particle='p', dE0=0.55, N=1E4, max_steps=100)
    return full_stack

def stack_30(E_p):
    stack = [
        #compartment08
       {'compound':'Ni', 'name':'Ni08', 'ad':ad_ni[7]},
       {'compound':'Sn', 'name':'Sn08', 'ad':ad_sn[7]},
       {'compound':'Ta', 'name':'Ta08', 'ad':ad_ta[7]},
       {'compound':'Cu', 'name':'Cu08', 'ad':ad_cu[7]},
       {'compound':'Al', 'name':'Al_degrader_E2_E1', 'ad':2*ad_degrader_e},
       # compartment09
       {'compound':'Ni', 'name':'Ni09', 'ad':ad_ni[8]},
       {'compound':'Sn', 'name':'Sn09', 'ad':ad_sn[8]},
       {'compound':'Ta', 'name':'Ta09', 'ad':ad_ta[8]},
       {'compound':'Cu', 'name':'Cu09', 'ad':ad_cu[8]},
       {'compound':'Al', 'name':'Al_degrader_E4_E3', 'ad':2*ad_degrader_e},
       # compartment10
       {'compound':'Ni', 'name':'Ni10', 'ad':ad_ni[9]},
       {'compound':'Sn', 'name':'Sn10', 'ad':ad_sn[9]},
       {'compound':'Ta', 'name':'Ta10', 'ad':ad_ta[9]},
       {'compound':'Cu', 'name':'Cu10', 'ad':ad_cu[9]},
       {'compound':'Al', 'name':'Al_degrader_E6_E5', 'ad':2*ad_degrader_e},
       # compartment11
       {'compound':'Ni', 'name':'Ni11', 'ad':ad_ni[10]}, 
       {'compound':'Sn', 'name':'Sn11', 'ad':ad_sn[10]},
       {'compound':'Ta', 'name':'Ta11', 'ad':ad_ta[10]},
       {'compound':'Cu', 'name':'Cu11', 'ad':ad_cu[10]},
       {'compound':'Al', 'name':'Al_degrader_E8_E7', 'ad':2*ad_degrader_e},
       # compartment12
       {'compound':'Ni', 'name':'Ni12', 'ad':ad_ni[11]}, 
       {'compound':'Sn', 'name':'Sn12', 'ad':ad_sn[11]},
       {'compound':'Ta', 'name':'Ta12', 'ad':ad_ta[11]},
       {'compound':'Cu', 'name':'Cu12', 'ad':ad_cu[11]},
       {'compound':'Al', 'name':'Al_degrader_E9', 'ad':ad_degrader_e},
       # compartment13
       {'compound':'Ni', 'name':'Ni13', 'ad':ad_ni[12]}, 
       {'compound':'Sn', 'name':'Sn13', 'ad':ad_sn[12]},
       {'compound':'Ta', 'name':'Ta13', 'ad':ad_ta[12]},
       {'compound':'Cu', 'name':'Cu13', 'ad':ad_cu[12]},
       {'compound':'Al', 'name':'Al_degrader_E10', 'ad':ad_degrader_e},
       # compartment14
       {'compound':'Ni', 'name':'Ni14', 'ad':ad_ni[13]},
       {'compound':'Sn', 'name':'Sn14', 'ad':ad_sn[13]},
       {'compound':'Ta', 'name':'Ta14', 'ad':ad_ta[13]},
       {'compound':'Cu', 'name':'Cu14', 'ad':ad_cu[13]},
    #    {'compound':'Al', 'name':'Al_degrader_E14', 'ad':ad_degrader_e},
    ]
    full_stack = ci.Stack(stack, E0=E_p, particle='p', dE0=0.55, N=1E4, max_steps=100)
    return full_stack

# def getAverageStackEnergy(stack, foil, numberOfFoils):
#     # st = stack(beamEnergy)
#     foils = []; E = []
#     for i in range(numberOfFoils):
#         foilNumber =  foil + '0' + str(i+1)
#         foils.append(foilNumber)
#         av_E = np.mean(stack.get_flux(foilNumber)[0])
#         E.append(int(av_E))
#     data = [foils, E]
#     return data


# def areal_density_Ga(properties):
#     #properties = {'Ga08': ['diameter', 'weight'], 'Ga09': ['diameter', 'weight']}
#     areal_densities = {}
#     for key,value in properties.items():
#         diameter = value[0]/10; weight = value[1]*1000 # mm-->cm, g-->mg
#         radius = diameter/2.0
#         area = np.pi * radius**2
#         areal_density = (weight/area) # mg/cm^2
#         areal_densities[key] = areal_density

#     return areal_densities


# areal_density_Ga = areal_density_Ga(
#     {'Ga1': [22.30, 0.0615],'Ga2' : [22.05, 0.2379], 'Ga3': [22.11, 0.4085], 'Ga4' : [22.05, 0.2154],'Ga5' : [16.14, 0.1067],'Ga6': [16.07, 0.0438], 'Ga7': [21.70, 0.7679], 'Ga8' : [22.59, 0.0814] ,
#         'Ga08_frame': [21.92, 0.3843],'Ga09_frame': [21.96, 0.0859],'Ga10': [22.10, 0.2713], 'Ga11': [22.46, 0.1014],'Ga12': [22.23, 0.6328],'Ga13': [1, 0], 'Ga14': [22.19, 0.2169]}
#     )



# ad_bare_kapton = 3.546473206 #mg/cm^2
# ad_poly = 6.668393219
# ad_kapton_dot = 8.169806242

# print(areal_density_Ga)

# subtract_from_Ga = [2*ad_bare_kapton, 2*ad_bare_kapton, 2*ad_bare_kapton, 2*ad_bare_kapton, 2*ad_bare_kapton, 2*ad_bare_kapton, 0, 1*ad_kapton_dot+1*ad_poly, 0, 0, 0, 0,0,0]

# print(subtract_from_Ga)
# # subtracted_ad = {}

# # for key,value in areal_density_Ga.items():
# #     subtracted_ad[key] = value - 2*ad_bare_kapton
    







# E_55 = 55; E_30 = 30
# E_Ni_55 = getAverageStackEnergy(stack_55, 'Ni', 7)
# E_Ni_30 = getAverageStackEnergy(stack_30, 'Ni', 7)
# E_Cu_55 = getAverageStackEnergy(stack_55, 'Cu', 7)
# E_Cu_30 = getAverageStackEnergy(stack_30, 'Cu', 7)
# print(E_Ni_55[1], E_Ni_30[1])
# print(E_Cu_55[1], E_Cu_30[1])

# print(stack_55.summarize())
# stack_30.plot('Ga')
# print(stack_30.summarize())
# stack_30.saveas('lbnl_TaSn_stack_30MeV.db')

# print(stack_55.summarize())
# print(type(stack_55))
# stack_55.saveas('analysis_TaSn_stack_55MeV.db')

# full_stack = pd.concat([stack_55, stack_30], ignore_index=True)
# full_stack.saveas('lbnl_TaSn_full_stack_55_30.db')



stack_55 = stack_55(55)
stack_30 = stack_30(30)
# print(stack_30.summarize())
# stack_55.saveas(save_stack_dir + 'TaSn_stack_55MeV_-5%.csv')
stack_55.saveas(save_stack_dir + 'TaSn_stack_55MeV_-25%.csv')
# stack_30.saveas(save_stack_dir + 'TaSn_stack_30MeV_-5%.csv')
stack_30.saveas(save_stack_dir + 'TaSn_stack_30MeV_-25%.csv')

