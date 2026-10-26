import matplotlib.pylab as plt


# For this, you'll need to look up the excitation energy for each of your isomers in Empire, 
# found in your defaults.lst file.   Try doing Control+F to look for "is an isomer", and grab the Ex
# once you find it for a given isomer
# 
# # Syntax:  'Ex': 'isotope'
# products= {'0.0895':'50-Sn-119', '0.1636':'51-Sb-122'}

# # Target isotope - needed if you have multiple sub-directories for each target isotope in an element
# target = '123Sb'

def Ta():
    #TA - non isomers
    # . Lacking: 179mW, 180mHf
    products_isomer = {'0.2219': '74-W - 179', '0.1000': '73-Ta-178',  '1.1415': '72-Hf -180', '0.3750': '72-Hf-179', '0.1238': '71-Lu-177'}
    # products = ['74-W -181', '74-W -179', '74-W -178', '74-W -177', '74-W -176',
    #            '73-Ta-180', '73-Ta-179', '73-Ta-178', '73-Ta-177', '73-Ta-176', '73-Ta-175',
    #            '72-Hf-175', '72-Hf-173',
    #         ]
    # Target isotope - needed if you have multiple sub-directories for each target isotope in an element
    target = 'Ta/181Ta'
    return products_isomer, target

def Cu(isotope):
    # #Cu - non isomers. Lacking: 
    products_isomer = {'0.0249': '27-Co- 58', '0.0586': '27-Co- 60'}
    products= ['30-Zn- 65', '30-Zn- 63', '30-Zn- 62', '29-Cu- 61', '27-Co- 58', '27-Co- 56', # monitors
            '29-Cu- 64', '29-Cu- 60',
            '28-Ni- 57',
            '27-Co- 61', '27-Co- 60', '27-Co- 55',
            '26-Fe- 59',
            '25-Mn- 56', '25-Mn- 54'
            ]
    if isotope == '63Cu':
        target = 'Cu/63Cu'
    elif isotope == '65Cu':
        target = 'Cu/65Cu'
    else:
        raise Exception("No valid target")
    return products_isomer, target

def Ni(isotope):
    # #Ni - non isomers. Lacking: 
    products_isomer = {'0.0249': '27-Co- 58', '0.0586': '27-Co- 60','0.3777': '25-Mn- 52'}
    products= ['28-Ni- 58', '27-Co- 57', '27-Co- 56', '27-Co- 55', '28-Ni- 56', '28-Ni- 57', # monitors
            '29-Cu- 64', '29-Cu- 61',
            '25-Mn- 54', '25-Mn- 52',
            '24-Cr- 51',
            ]
    if isotope == '58Ni':
        target = 'Ni/58Ni'
    elif isotope == '60Ni':
        target = 'Ni/60Ni'
    elif isotope == '61Ni':
        target = 'Ni/61Ni'
    elif isotope == '62Ni':
        target = 'Ni/62Ni'
    elif isotope == '64Ni':
        target = 'Ni/64Ni'
    else:
        raise Exception("No valid target")
    return products_isomer, target

# products, target = Ta()
# products, target = Cu('63Cu')
# products, target = Cu('65Cu')
# products, target = Ni('58Ni')
# products, target = Ni('60Ni')
# products, target = Ni('61Ni')
# products, target = Ni('62Ni')
# products, target = Ni('64Ni')


with open('./'+target+'/defaults.lst', 'r') as f:

    # energies = [2.0]
    energies = []
    cross_sections = dict.fromkeys(products.keys())
    reading = False

    for line in f:
        if ' incident energy  ' in line and '(CMS)' not in line and 'C.M.' not in line :
            # print(line)
            # current_energy = float(line[18:28].replace('D',"E"))
            current_energy = float(line[40:].replace('D',"E").replace("MeV (LAB)\n",""))
            # print(current_energy)
            if len(energies) == 0 or current_energy != energies[-1]:
                energies.append(current_energy)
        elif 'is an isomer' in line and 'Level of energy' in line:
            energy = line.split()[3]
            #print(energy)
            if energy in cross_sections.keys():
                # print(line)
                if cross_sections[energy] is None:
                    cross_sections[energy] = []
                while len(cross_sections[energy])<len(energies)-1:
                    cross_sections[energy].append(0.0)
                cross_sections[energy].append(float(line.split()[11]))





for isotope in products:
    if cross_sections[isotope] is None:
        print('nothing for ', isotope)
        continue
    #print('plotting ', isotope)
    # print(energies)
    # print(cross_sections[isotope])
    try:
        plt.plot(energies,cross_sections[isotope])
        plt.xlabel('Proton Energy (MeV)')
        plt.ylabel('Cross Section (mb)')
        plt.title(products[isotope]+" M")
        plt.savefig(target+'/plots/'+products[isotope]+"M_empire.png")
        plt.close()
    except ValueError:
        continue
    with open('./'+target+'/plots/'+products[isotope].replace(" ","")+"M_empire.txt",'w') as f:
        for i in range(len(energies)):
            f.write(str(energies[i])+"\t"+str(cross_sections[isotope][i])+"\n")
