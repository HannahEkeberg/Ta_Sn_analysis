import numpy as np
import curie as ci
import pandas as pd

x_kapton = 0.013
x_silicone = 0.013
ad_degrader_a = 599.   #2.24 mm
ad_degrader_b = 415.0   #1.55 mm
ad_degrader_c = 261.5   #0.97 mm
# ad_degrader_d = 599.0
ad_degrader_e = 68.3   #0.256 mm
ad_degrader_h = 33.8
ad_be_backing = 4.425     #23.9130435 microns




def stack_55(E_p):
    stack = [
        #compartment01
       {'compound':'Ni', 'name':'Ni01', 't':0.025},
       {'compound':'Sn', 'name':'Sn01', 't':0.050},
       {'compound':'Ta', 'name':'Ta01', 't':0.022},
       {'compound':'Cu', 'name':'Cu01', 't':0.025},
       {'compound':'Al', 'name':'Al_degrader_A1', 'ad':ad_degrader_a},
       # compartment02
       {'compound':'Ni', 'name':'Ni02', 't':0.025}, 
       {'compound':'Sn', 'name':'Sn02', 't':0.048},
       {'compound':'Ta', 'name':'Ta02', 't':0.022},
       {'compound':'Cu', 'name':'Cu02', 't':0.025},
       {'compound':'Al', 'name':'Al_degrader_B2', 'ad':ad_degrader_b},
       # compartment03
       {'compound':'Ni', 'name':'Ni03', 't':0.025}, 
       {'compound':'Sn', 'name':'Sn03', 't':0.055},
       {'compound':'Ta', 'name':'Ta03', 't':0.020},
       {'compound':'Cu', 'name':'Cu03', 't':0.025},
       {'compound':'Al', 'name':'Al_degrader_C5', 'ad':ad_degrader_c},
       # compartment04
       {'compound':'Ni', 'name':'Ni04', 't':0.025}, 
       {'compound':'Sn', 'name':'Sn04', 't':0.012},
       {'compound':'Ta', 'name':'Ta04', 't':0.021},
       {'compound':'Cu', 'name':'Cu04', 't':0.025},
       {'compound':'Al', 'name':'Al_degrader_E8_E1', 'ad':2*ad_degrader_e},
       # compartment05
       {'compound':'Ni', 'name':'Ni05', 't':0.025}, 
       {'compound':'Sn', 'name':'Sn05', 't':0.012},
       {'compound':'Ta', 'name':'Ta05', 't':0.050},
       {'compound':'Cu', 'name':'Cu05', 't':0.025},
       {'compound':'Al', 'name':'Al_degrader_E4_E3_E2', 'ad':3*ad_degrader_e},
       # compartment06
       {'compound':'Ni', 'name':'Ni06', 't':0.025}, 
       {'compound':'Sn', 'name':'Sn06', 't':0.013},
       {'compound':'Ta', 'name':'Ta06', 't':0.021},
       {'compound':'Cu', 'name':'Cu06', 't':0.025},
       {'compound':'Al', 'name':'Al_degrader_E7_E6_E5', 'ad':3*ad_degrader_e},
       # compartment07
       {'compound':'Ni', 'name':'Ni07', 't':0.025}, 
       {'compound':'Sn', 'name':'Sn07', 't':0.013},
       {'compound':'Ta', 'name':'Ta07', 't':0.021},
       {'compound':'Cu', 'name':'Cu07', 't':0.025},
    ]
    full_stack = ci.Stack(stack, E0=E_p, particle='p', dE0=0.55, N=1E4, max_steps=100)
    return full_stack

def stack_30(E_p):
    stack = [
        #compartment08
       {'compound':'Ni', 'name':'Ni08', 't':0.025}, 
       {'compound':'Sn', 'name':'Sn08', 't':0.012},
       {'compound':'Ta', 'name':'Ta08', 't':0.021},
       {'compound':'Cu', 'name':'Cu08', 't':0.025},
       {'compound':'Al', 'name':'Al_degrader_E2_E1', 'ad':2*ad_degrader_e},
       # compartment09
       {'compound':'Ni', 'name':'Ni09', 't':0.025}, 
       {'compound':'Sn', 'name':'Sn09', 't':0.012},
       {'compound':'Ta', 'name':'Ta09', 't':0.021},
       {'compound':'Cu', 'name':'Cu09', 't':0.025},
       {'compound':'Al', 'name':'Al_degrader_E4_E3', 'ad':2*ad_degrader_e},
       # compartment10
       {'compound':'Ni', 'name':'Ni10', 't':0.025}, 
       {'compound':'Sn', 'name':'Sn10', 't':0.012},
       {'compound':'Ta', 'name':'Ta10', 't':0.022},
       {'compound':'Cu', 'name':'Cu10', 't':0.025},
       {'compound':'Al', 'name':'Al_degrader_E6_E5', 'ad':2*ad_degrader_e},
       # compartment11
       {'compound':'Ni', 'name':'Ni11', 't':0.025}, 
       {'compound':'Sn', 'name':'Sn11', 't':0.014},
       {'compound':'Ta', 'name':'Ta11', 't':0.021},
       {'compound':'Cu', 'name':'Cu11', 't':0.025},
       {'compound':'Al', 'name':'Al_degrader_E8_E7', 'ad':2*ad_degrader_e},
       # compartment12
       {'compound':'Ni', 'name':'Ni12', 't':0.025}, 
       {'compound':'Sn', 'name':'Sn12', 't':0.014},
       {'compound':'Ta', 'name':'Ta12', 't':0.021},
       {'compound':'Cu', 'name':'Cu12', 't':0.025},
       {'compound':'Al', 'name':'Al_degrader_E9', 'ad':ad_degrader_e},
       # compartment13
       {'compound':'Ni', 'name':'Ni13', 't':0.025}, 
       {'compound':'Sn', 'name':'Sn13', 't':0.011},
       {'compound':'Ta', 'name':'Ta13', 't':0.021},
       {'compound':'Cu', 'name':'Cu13', 't':0.025},
       {'compound':'Al', 'name':'Al_degrader_E10', 'ad':ad_degrader_e},
       # compartment14
       {'compound':'Ni', 'name':'Ni14', 't':0.025},
       {'compound':'Sn', 'name':'Sn14', 't':0.011},
       {'compound':'Ta', 'name':'Ta14', 't':0.021},
       {'compound':'Cu', 'name':'Cu14', 't':0.025},
    #    {'compound':'Al', 'name':'Al_degrader_E14', 'ad':ad_degrader_e},
    ]
    full_stack = ci.Stack(stack, E0=E_p, particle='p', dE0=0.55, N=1E4, max_steps=100)
    return full_stack

def getAverageStackEnergy(stack, foil, numberOfFoils):
    # st = stack(beamEnergy)
    foils = []; E = []
    for i in range(numberOfFoils):
        foilNumber =  foil + '0' + str(i+1)
        foils.append(foilNumber)
        av_E = np.mean(stack.get_flux(foilNumber)[0])
        E.append(int(av_E))
    data = [foils, E]
    return data


def areal_density_Ga(properties):
    #properties = {'Ga08': ['diameter', 'weight'], 'Ga09': ['diameter', 'weight']}
    areal_densities = {}
    for key,value in properties.items():
        diameter = value[0]/10; weight = value[1]*1000 # mm-->cm, g-->mg
        radius = diameter/2.0
        area = np.pi * radius**2
        areal_density = (weight/area) # mg/cm^2
        areal_densities[key] = areal_density

    return areal_densities


areal_density_Ga = areal_density_Ga(
    {'Ga1': [22.30, 0.0615],'Ga2' : [22.05, 0.2379], 'Ga3': [22.11, 0.4085], 'Ga4' : [22.05, 0.2154],'Ga5' : [16.14, 0.1067],'Ga6': [16.07, 0.0438], 'Ga7': [21.70, 0.7679], 'Ga8' : [22.59, 0.0814] ,
        'Ga08_frame': [21.92, 0.3843],'Ga09_frame': [21.96, 0.0859],'Ga10': [22.10, 0.2713], 'Ga11': [22.46, 0.1014],'Ga12': [22.23, 0.6328],'Ga13': [1, 0], 'Ga14': [22.19, 0.2169]}
    )



ad_bare_kapton = 3.546473206 #mg/cm^2
ad_poly = 6.668393219
ad_kapton_dot = 8.169806242

print(areal_density_Ga)

subtract_from_Ga = [2*ad_bare_kapton, 2*ad_bare_kapton, 2*ad_bare_kapton, 2*ad_bare_kapton, 2*ad_bare_kapton, 2*ad_bare_kapton, 0, 1*ad_kapton_dot+1*ad_poly, 0, 0, 0, 0,0,0]

print(subtract_from_Ga)
# subtracted_ad = {}

# for key,value in areal_density_Ga.items():
#     subtracted_ad[key] = value - 2*ad_bare_kapton
    







# E_55 = 55; E_30 = 30
stack_55 = stack_55(55)
# stack_30 = stack_30(30)
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
stack_55.saveas('lbnl_TaSn_stack_55MeV.db')


# full_stack = pd.concat([stack_55, stack_30], ignore_index=True)
# full_stack.saveas('lbnl_TaSn_full_stack_55_30.db')