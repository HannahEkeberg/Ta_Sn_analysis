from check_gammas import *

# MONITORS
gammas_65Zn = [1115.539]
gammas_63Zn = [449.93, 669.62, 742.25, 962.06, 1123.72, 1327.03, 1392.55, 1412.08, 1547.04, 2026.8, 2336.5]
gammas_62Zn = [246.95, 260.43, 349.6, 394.03, 548.35, 596.56, 637.41]
gammas_61Cu = [246.95, 282.956, 373.05, 529.169, 588.605, 656.008, 816.692, 841.211, 902.294, 1132.351, 1662.0]
gammas_58Co = [810.7593]
gammas_56Co = [1360.212, 1771.357, 2015.215, 2034.791]
gammas_57Ni = [127.164, 1377.63, 1757.55]
gammas_56Ni = [158.38, 269.5, 480.44, 749.95, 811.85, 1561.8]
gammas_57Co = [122.06065, 136.47356]
gammas_55Co = [931.1, 1316.6, 1370.0, 1408.5]



#Cu reactions: 
gammas_60Cu = [1791.6, 826.4]
gammas_57Ni = [127.164, 1377.63, 1757.55]
gammas_61Co = [917.5, 67.412]
gammas_61Co = [67.412]
gammas_56Mn = [846.7638, 1810.726]
gammas_56Mn = [1810.726]

# Ni reactions:
gammas_61Cu = [282.956, 373.05, 529.169, 588.605, 816.692, 841.211]
gammas_52Mn = [744.233, 935.544, 1246.278,1434.092]

def keep_selected(gammas,isotope, element):
    for foil in foils(element):
        try:
            keep_selected_gammas(foil=foil, isotope=isotope, gammas=gammas, unmark_first=True)
        except:
            print("file not found in foil " + foil)

# keep_selected(gammas_60Cu, '60CUg', 'Cu')
# keep_selected(gammas_57Ni, '57NIg', 'Cu')
# keep_selected(gammas_61Co, '61COg', 'Cu')
# keep_selected(gammas_56Mn, '56MNg', 'Cu')
# keep_selected(gammas_61Cu, '61CUg', 'Ni')
# keep_selected(gammas_52Mn, '52MNg', 'Ni')


# keep_selected(gammas_63Zn, '63ZNg', 'Cu')
# keep_selected(gammas_65Zn, '65ZNg', 'Cu')
# keep_selected(gammas_62Zn, '62ZNg', 'Cu')
# keep_selected(gammas_61Cu, '61CUg', 'Cu')
# keep_selected(gammas_58Co, '58COg', 'Cu')
# keep_selected(gammas_58Co, '58COg', 'Ni')
# keep_selected(gammas_56Co, '56COg', 'Cu')
# keep_selected(gammas_56Co, '56COg', 'Ni')
# keep_selected(gammas_57Ni, '57NIg', 'Ni')
# keep_selected(gammas_56Ni, '56NIg', 'Ni')
# keep_selected(gammas_57Co, '57COg', 'Ni')
# keep_selected(gammas_55Co, '55COg', 'Ni')