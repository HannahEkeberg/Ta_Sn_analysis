import numpy as np
import curie as ci
from scipy.constants import N_A
import pandas as pd
import os

def areal_density_squarefoils(side1, side2, thickness, weight):
    """
    Should return g/cm2
    """
    side1_mean = np.mean(side1) / 10 # mm --> cm
    side2_mean = np.mean(side2) / 10 # mm --> cm
    thickness_mean = np.mean(thickness) # mm
    weight_mean = np.mean(weight) # g

    side1_std = np.std(side1) / 10 # mm --> cm
    side2_std = np.std(side2) / 10 # mm --> cm
    thickness_std = np.std(thickness) # mm
    weight_std = np.std(weight) # g

    print(side1_mean, side2_mean, weight_mean, thickness_mean*1e3)

    areal_density_mean = weight_mean / (side1_mean * side2_mean) # g/cm^2
    areal_density_std = np.sqrt(
        (weight_std/weight_mean)**2 
        + (side1_std/side1_mean)**2 
        + (side2_std/side2_mean)**2 
    ) * areal_density_mean  # g/cm^2

    return areal_density_mean, areal_density_std

def areal_density_cirlefoils(diameter, thickness, weight):
    """
    Should return g/cm2
    """
    diameter_mean = np.mean(diameter) / 10 # mm --> cm
    thickness_mean = np.mean(thickness) # mm
    weight_mean = np.mean(weight) # g
    diameter_std = np.std(diameter) / 10 # mm --> cm
    thickness_std = np.std(thickness) # mm
    weight_std = np.std(weight) # g

    areal_density_mean = weight_mean / (np.pi * (diameter_mean/2)**2) # g/cm^2
    areal_density_std = np.sqrt(
        (weight_std/weight_mean)**2 
        + (diameter_std/diameter_mean)**2 
    ) * areal_density_mean  # g/cm^2

    return areal_density_mean, areal_density_std

def areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element):
    molar_mass = ci.Element(element).mass # g/mol
    # N_A ~ 1/mol
    numb_of_target_nuclei = (areal_density*N_A)/molar_mass #nuclei/cm^2
    unc_numb_of_target_nuclei = (unc_areal_density*N_A)/molar_mass
    return numb_of_target_nuclei, unc_numb_of_target_nuclei


def assembleDataframe(listOfFoilProperties, element):
    data = []
    for i in range(len(listOfFoilProperties)):
        data.append(listOfFoilProperties[i])
    df = pd.DataFrame(data, columns=['foil', 'nuclei/cm2', 'unc nuclei/cm2', 'arealdensity g/cm2', 'unc arealdensity g/cm2'])
    df.to_csv('./generatedfiles/arealdensity/areal_density_' + element + '.csv')


def Cu01():
    foilshape = 'square'
    element ='Cu'; foil = 'Cu01'
    side1 = [25.68,25.67,25.65,25.63] # mm
    side2 = [24.91,24.98,24.80,24.78] # mm
    thickness = [0.029,0.030,0.029,0.030] # mm
    weight = [0.1398,0.1398,0.1397] # mg
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density


def Cu02():
    foilshape = 'square'
    element ='Cu'; foil = 'Cu02'
    side1 = [25.32,25.36,25.33,25.34]
    side2 = [25.33,25.38,25.35,25.28]
    thickness = [0.030,0.030,0.030,0.030]
    weight = [0.1428,0.1427,0.1427]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Cu03():
    foilshape = 'square'
    element ='Cu'; foil = 'Cu03'
    side1 = [25.41,25.39,25.38,25.38]
    side2 = [25.20,25.28,25.35,25.36]
    thickness = [0.029,0.030,0.028,0.030]
    weight = [0.1430,0.1429,0.1429]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Cu04():
    foilshape = 'square'
    element ='Cu'; foil = 'Cu04'
    side1 = [25.35,25.37,25.37,25.37]
    side2 = [25.42,25.48,25.48,25.48]
    thickness = [0.029,0.029,0.029,0.029]
    weight = [0.1439,0.1440,0.1440]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Cu05():
    foilshape = 'square'
    element ='Cu'; foil = 'Cu05'
    side1 = [25.34,25.34,25.34,25.32] 
    side2 = [25.32,25.32,25.30,25.30]
    thickness = [0.030,0.030,0.029,0.030]
    weight = [0.1432,0.1430,0.1430] 
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Cu06():
    foilshape = 'square'
    element ='Cu'; foil = 'Cu06'
    side1 = [24.68,24.73,24.73,24.73]
    side2 = [25.32,25.36,25.41,25.45]
    thickness = [0.029,0.030,0.030,0.030]
    weight = [0.1396,0.1395,0.1394] 
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Cu07():
    foilshape = 'square'
    element ='Cu'; foil = 'Cu07'
    side1 = [25.43,25.41,25.38,25.36]
    side2 = [25.46,25.47,25.42,25.41]
    thickness = [0.030,0.031,0.029,0.030]
    weight = [0.1438,0.1438,0.1438] 
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Cu08():
    foilshape = 'square'
    element ='Cu'; foil = 'Cu08'
    side1 = [25.17,25.16,25.21,25.20]
    side2 = [25.18,25.18,25.18,25.18]
    thickness = [0.028,0.031,0.029,0.030]
    weight = [0.1407,0.1408,0.1408] 
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Cu09():
    foilshape = 'square'
    element ='Cu'; foil = 'Cu09'
    side1 = [25.11,25.16,25.17,25.18]
    side2 = [24.99,24.96,25.01,25.02]
    thickness = [0.029,0.028,0.030,0.029]
    weight = [0.1390,0.1389,0.1388] 
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Cu10():
    foilshape = 'square'
    element ='Cu'; foil = 'Cu10'
    side1 = [25.38,25.39,25.40,25.42]
    side2 = [25.31,25.33,25.26,25.28]
    thickness = [0.030,0.030,0.029,0.028]
    weight = [0.1430,0.1430,0.1431] 
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Cu11():
    foilshape = 'square'
    element ='Cu'; foil = 'Cu11'
    side1 = [24.71,24.80,24.91,24.93]
    side2 = [24.88,24.90,24.91,25.01]
    thickness = [0.031,0.030,0.030,0.028]
    weight = [0.1375,0.1375,0.1375] 
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Cu12():
    foilshape = 'square'
    element ='Cu'; foil = 'Cu12'
    side1 = [25.07,25.12,25.01,25.07]
    side2 = [24.89,24.92,24.94,25.00]
    thickness = [0.029,0.030,0.030,0.030]
    weight = [0.1388,0.1385,0.1386] 
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Cu13():
    foilshape = 'square'
    element ='Cu'; foil = 'Cu13'
    side1 = [24.94,24.95,24.97,24.97]
    side2 = [24.67,24.75,24.82,24.76]
    thickness = [0.030,0.029,0.029,0.030]
    weight = [0.1369,0.1369,0.1370] 
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Cu14():
    foilshape = 'square'
    element ='Cu'; foil = 'Cu14'
    side1 = [25.05,25.04,25.03,25.04]
    side2 = [24.71,24.75,24.93,24.96]
    thickness = [0.028,0.030,0.029,0.029]
    weight = [0.1376,0.1376,0.1376] 
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ni01():
    foilshape = 'square'
    element ='Ni'; foil = 'Ni01'
    side1 = [25.03,25.02,25.03,25.08]
    side2 = [25.02,25.01,25.06,25.05]
    thickness = [0.028,0.03,0.028,0.031]
    weight = [0.1447,0.1446,0.1447,0.1447] 
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ni02():
    foilshape = 'square'
    element ='Ni'; foil = 'Ni02'
    side1 = [25.32,25.32,25.32,25.3]
    side2 = [25.11,25.09,25.11,25.11]
    thickness = [0.031,0.03,0.03,0.03]
    weight = [0.1467,0.1467,0.1468,0.1468] 
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ni03():
    foilshape = 'square'
    element ='Ni'; foil = 'Ni03'
    side1 = [25.02,25.02,25.03,25.08]
    side2 = [24.82,24.85,24.96,24.95]
    thickness = [0.029,0.031,0.029,0.029]
    weight = [0.1432,0.1431,0.1432,0.1431]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ni04():
    foilshape = 'square'
    element ='Ni'; foil = 'Ni04'
    side1 = [24.85,24.87,24.91,24.96]
    side2 = [25.02,25,24.99,24.99]
    thickness = [0.029,0.03,0.03,0.03]
    weight = [0.1431,0.1431,0.143,0.143]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ni05():
    foilshape = 'square'
    element ='Ni'; foil = 'Ni05'
    side1 = [25.18,25.14,25.18,25.19]
    side2 = [25.05,25.09,25.05,25.05]
    thickness = [0.03,0.03,0.03,0.03]
    weight = [0.1447,0.1447,0.1447,0.1448]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ni06():
    foilshape = 'square'
    element ='Ni'; foil = 'Ni06'
    side1 = [25.1,25.15,25.16,25.19]
    side2 = [24.87,24.88,24.91,24.89]
    thickness = [0.028,0.03,0.03,0.03]
    weight = [0.1416,0.1417,0.1417,0.1416]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ni07():
    foilshape = 'square'
    element ='Ni'; foil = 'Ni07'
    side1 = [24.98,24.99,25.01,24.99]
    side2 = [25.4,25.35,25.32,25.36]
    thickness = [0.029,0.031,0.03,0.029]
    weight = [0.1428,0.1428,0.1428,0.1427]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ni08():
    foilshape = 'square'
    element ='Ni'; foil = 'Ni08'
    side1 = [24.83,24.82,24.8,24.83]
    side2 = [25.31,25.4,25.46,25.56]
    thickness = [0.03,0.03,0.029,0.029]
    weight = [0.1461,0.1461,0.1461,0.1461]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ni09():
    foilshape = 'square'
    element ='Ni'; foil = 'Ni09'
    side1 = [25.08,25.07,25.07,25.07]
    side2 = [25.06,25.06,25.04,25.07]
    thickness = [0.029,0.03,0.029,0.029]
    weight = [0.1423,0.1424,0.1424,0.1424]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ni10():
    foilshape = 'square'
    element ='Ni'; foil = 'Ni10'
    side1 = [25.01,25,25.02,25.02]
    side2 = [25.16,25.12,25.17,25.16]
    thickness = [0.029,0.029,0.03,0.03]
    weight = [0.1453,0.1454,0.1454,0.1454]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ni11():
    foilshape = 'square'
    element ='Ni'; foil = 'Ni11'
    side1 = [24.97,24.97,24.95,24.95]
    side2 = [25.18,25.23,25.23,25.25]
    thickness = [0.03,0.029,0.03,0.028]
    weight = [0.1445,0.1444,0.1444,0.1444]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ni12():
    foilshape = 'square'
    element ='Ni'; foil = 'Ni12'
    side1 = [25.18,25.2,25.21,25.21]
    side2 = [25.03,25.03,25.02,25.02]
    thickness = [0.03,0.03,0.029,0.03]
    weight = [0.146,0.146,0.146,0.1461]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ni13():
    foilshape = 'square'
    element ='Ni'; foil = 'Ni13'
    side1 = [24.97,24.98,24.99,25.01]
    side2 = [25.17,25.19,25.23,25.25]
    thickness = [0.03,0.03,0.03,0.03]
    weight = [0.1459,0.1459,0.1458,0.1458]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ni14():
    foilshape = 'square'
    element ='Ni'; foil = 'Ni14'
    side1 = [25.15,25.16,25.11,25.12]
    side2 = [25.22,25.24,25.24,25.25]
    thickness = [0.03,0.029,0.028,0.03]
    weight = [0.1463,0.1464,0.1463,0.1463]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ta01():
    foilshape = 'square'
    element ='Ta'; foil = 'Ta01'
    side1 = [25.18,25.21,25.22,25.21]
    side2 = [25.27,25.25,25.22,25.09]
    thickness = [0.021,0.021,0.021,0.021]
    weight = [0.1737,0.1736,0.1737]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ta02():
    foilshape = 'square'
    element ='Ta'; foil = 'Ta02'
    side1 = [25.2,25.24,25.28,25.31]
    side2 = [25.04,25.13,25.19,25.32]
    thickness = [0.021,0.021,0.022,0.021]
    weight = [0.1741,0.1742,0.1742]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ta03():
    foilshape = 'square'
    element ='Ta'; foil = 'Ta03'
    side1 = [25.59,25.57,25.55,25.54]
    side2 = [24.96,25.02,25.07,25.12]
    thickness = [0.022,0.022,0.022,0.022]
    weight = [0.1774,0.1773,0.1773]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ta04():
    foilshape = 'square'
    element ='Ta'; foil = 'Ta04'
    side1 = [25.29,25.33,25.37,25.42]
    side2 = [25.25,25.23,25.25,25.2]
    thickness = [0.021,0.021,0.021,0.021]
    weight = [0.176,0.176,0.1761]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ta05():
    foilshape = 'square'
    element ='Ta'; foil = 'Ta05'
    side1 = [25.69,25.76,25.86,25.94]
    side2 = [25.55,25.6,25.67,25.72]
    thickness = [0.021,0.02,0.021,0.02]
    weight = [0.1799,0.18,0.18]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ta06():
    foilshape = 'square'
    element ='Ta'; foil = 'Ta06'
    side1 = [25.45,25.41,25.39,25.36]
    side2 = [25.68,25.66,25.5,25.38]
    thickness = [0.021,0.021,0.021,0.02]
    weight = [0.1777,0.1777,0.1777]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ta07():
    foilshape = 'square'
    element ='Ta'; foil = 'Ta07'
    side1 = [25.38,25.38,25.35,25.29]
    side2 = [25.18,25.15,25.13,25.06]
    thickness = [0.021,0.021,0.021,0.021]
    weight = [0.1755,0.1754,0.1754]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ta08():
    foilshape = 'square'
    element ='Ta'; foil = 'Ta08'
    side1 = [25.69,25.76,25.84,25.96]
    side2 = [25.46,25.42,25.4,25.37]
    thickness = [0.021,0.022,0.022,0.02]
    weight = [0.1786,0.1786,0.1786]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ta09():
    foilshape = 'square'
    element ='Ta'; foil = 'Ta09'
    side1 = [25.42,25.48,25.49,25.53]
    side2 = [25.42,25.35,25.3,25.33]
    thickness = [0.022,0.022,0.022,0.022]
    weight = [0.1768,0.1768,0.1769]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ta10():
    foilshape = 'square'
    element ='Ta'; foil = 'Ta10'
    side1 = [24.71,24.8,24.88,25.03]
    side2 = [25.28,25.27,25.27,25.23]
    thickness = [0.02,0.02,0.02,0.02]
    weight = [0.1706,0.1706,0.1705]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ta11():
    foilshape = 'square'
    element ='Ta'; foil = 'Ta11'
    side1 = [25.1,25.05,25,24.93]
    side2 = [25.34,25.41,25.49,25.58]
    thickness = [0.021,0.022,0.021,0.021]
    weight = [0.1733,0.1732,0.1732]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ta12():
    foilshape = 'square'
    element ='Ta'; foil = 'Ta12'
    side1 = [25.5,25.62,25.61,25.77]
    side2 = [25.58,25.64,25.67,25.72]
    thickness = [0.022,0.022,0.021,0.022]
    weight = [0.1821,0.1821,0.1821]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ta13():
    foilshape = 'square'
    element ='Ta'; foil = 'Ta13'
    side1 = [25.32,25.31,25.32,25.31]
    side2 = [25.26,25.24,25.24,25.18]
    thickness = [0.021,0.021,0.02,0.021]
    weight = [0.1746,0.1745,0.1747]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Ta14():
    foilshape = 'square'
    element ='Ta'; foil = 'Ta14'
    side1 = [25.41,25.39,25.35,25.26]
    side2 = [25.33,25.35,25.35,25.48]
    thickness = [0.022,0.022,0.021,0.021]
    weight = [0.176,0.1761,0.176]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density


def Sn01():
    foilshape = 'square'
    element ='Sn'; foil = 'Sn01'
    thickness = [0.050, 0.047, 0.051, 0.051]
    weight = [0.2250, 0.2249, 0.2249, 0.2249]
    side1 = [23.76, 23.98, 24.19, 24.41]
    side2 = [24.92, 24.94, 24.88, 24.94]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Sn02():
    foilshape = 'square'
    element ='Sn'; foil = 'Sn02'
    thickness = [0.048, 0.048, 0.049, 0.048]
    weight = [0.2422, 0.2423, 0.2423, 0.2422]
    side1 = [24.87, 24.91, 24.96, 25.08]
    side2 = [26.04, 26.09, 26.04, 26.18]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Sn03():
    foilshape = 'square'
    element ='Sn'; foil = 'Sn03'
    thickness = [0.054, 0.054, 0.055, 0.055]
    weight = [0.2345, 0.2344, 0.2344, 0.2344]
    side1 = [24.87, 24.99, 25.00, 24.99]
    side2 = [25.30, 25.39, 25.26, 25.22]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Sn04():
    foilshape = 'square'
    element ='Sn'; foil = 'Sn04'
    thickness = [0.055, 0.055, 0.054, 0.055]
    weight = [0.2293, 0.2293, 0.2292, 0.2293]
    side1 = [24.61, 25.05, 24.88, 24.98]
    side2 = [24.98, 24.93, 24.86, 24.88]
    areal_density, unc_areal_density = areal_density_squarefoils(side1, side2, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Sn05():
    foilshape = 'circle'
    element ='Sn'; foil = 'Sn05'
    thickness = [0.012, 0.012, 0.012, 0.013]
    weight = [0.0391, 0.0391, 0.0392, 0.0392]
    diameter = [24.96, 24.89, 24.97, 25.00]
    areal_density, unc_areal_density = areal_density_cirlefoils(diameter, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Sn06():
    foilshape = 'circle'
    element ='Sn'; foil = 'Sn06'
    thickness = [0.013, 0.013, 0.014, 0.013]
    weight = [0.0381, 0.0382, 0.0382, 0.0383]
    diameter = [24.87, 24.97, 24.92, 25.02]
    areal_density, unc_areal_density = areal_density_cirlefoils(diameter, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Sn07():
    foilshape = 'circle'
    element ='Sn'; foil = 'Sn07'
    thickness = [0.013, 0.013, 0.013, 0.013]
    weight = [0.0388, 0.0387, 0.0388, 0.388]
    diameter = [24.95, 25.01, 24.91, 24.99]
    areal_density, unc_areal_density = areal_density_cirlefoils(diameter, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Sn08():
    foilshape = 'circle'
    element ='Sn'; foil = 'Sn08'
    thickness = [0.012, 0.012, 0.012, 0.012]
    weight = [0.0393, 0.0393, 0.0394, 0.0394]
    diameter = [24.98, 24.99, 24.94, 24.90]
    areal_density, unc_areal_density = areal_density_cirlefoils(diameter, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Sn09():
    foilshape = 'circle'
    element ='Sn'; foil = 'Sn09'
    thickness = [0.012, 0.012, 0.012, 0.012]
    weight = [0.0384, 0.0384, 0.0384, 0.0383]
    diameter = [24.98, 24.99, 24.98, 24.98]
    areal_density, unc_areal_density = areal_density_cirlefoils(diameter, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Sn10():
    foilshape = 'circle'
    element ='Sn'; foil = 'Sn10'
    thickness = [0.012, 0.012, 0.013, 0.013]
    weight = [0.0404, 0.0404, 0.0405, 0.0405]
    diameter = [25.02, 25.01, 25.02, 25.00]
    areal_density, unc_areal_density = areal_density_cirlefoils(diameter, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Sn11():
    foilshape = 'circle'
    element ='Sn'; foil = 'Sn11'
    thickness = [0.015, 0.014, 0.014, 0.014]
    weight = [0.0391, 0.0390, 0.0390, 0.0390]
    diameter = [24.89, 25.00, 24.98, 24.93]
    areal_density, unc_areal_density = areal_density_cirlefoils(diameter, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Sn12():
    foilshape = 'circle'
    element ='Sn'; foil = 'Sn12'
    diameter = [24.91,24.97,24.92,24.85]
    thickness = [0.013,0.013,0.014,0.012]
    weight = [0.039,0.0392,0.0391,0.0391]
    areal_density, unc_areal_density = areal_density_cirlefoils(diameter, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Sn13():
    foilshape = 'circle'
    element ='Sn'; foil = 'Sn13'
    diameter = [24.99,24.98,24.99,24.87]
    thickness = [0.011,0.011,0.011,0.011]
    weight = [0.0361,0.036,0.0359,0.0359]
    areal_density, unc_areal_density = areal_density_cirlefoils(diameter, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density

def Sn14():
    foilshape = 'circle'
    element ='Sn'; foil = 'Sn14'
    diameter = [24.99,24.97,24.97,25.00]
    thickness = [0.011,0.012,0.011,0.011]
    weight = [0.0374,0.0374,0.0373,0.0373]
    areal_density, unc_areal_density = areal_density_cirlefoils(diameter, thickness, weight)
    numb_of_target_nuclei, unc_numb_of_target_nuclei = areal_density_to_number_of_target_nuclei(areal_density, unc_areal_density, element)
    return foil, numb_of_target_nuclei, unc_numb_of_target_nuclei, areal_density, unc_areal_density



# assembleDataframe([Cu01(), Cu02(), Cu03(), Cu04(), Cu05(), Cu06(),Cu07(), Cu08(), Cu09(), Cu10(), Cu11(), Cu12(), Cu13(), Cu14()], 'Cu')
# assembleDataframe([Ni01(), Ni02(), Ni03(), Ni04(), Ni05(), Ni06(),Ni07(), Ni08(), Ni09(), Ni10(), Ni11(), Ni12(), Ni13(), Ni14()], 'Ni')
# assembleDataframe([Ta01(), Ta02(), Ta03(), Ta04(), Ta05(), Ta06(),Ta07(), Ta08(), Ta09(), Ta10(), Ta11(), Ta12(), Ta13(), Ta14()], 'Ta')
# assembleDataframe([Sn01(), Sn02(), Sn03(), Sn04(), Sn05(), Sn06(),Sn07(), Sn08(), Sn09(), Sn10(), Sn11(), Sn12(), Sn13(), Sn14()], 'Sn')

# Cu01()