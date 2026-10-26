import matplotlib.pylab as plt
import os 


# A list of all isotopes you want extracted for plotting
# products= ['52-Te-116', '52-Te-117', '52-Te-118', '52-Te-119', '52-Te-120', '52-Te-121', '52-Te-122', 
#            '51-Sb-115', '51-Sb-116', '51-Sb-117', '51-Sb-118', '51-Sb-119', '51-Sb-120', '51-Sb-121', 
#            '50-Sn-112', '50-Sn-113', '50-Sn-114', '50-Sn-115', '50-Sn-116', '50
# -Sn-117', '50-Sn-118', '50-Sn-119', '50-Sn-120' #, 
#            # etc etc
#            ]

def Ta():
    #TA - non isomers
    # . Lacking: 179mW, 180mHf
    products= [
            #'74-W -181', '74-W -179', '74-W -178', '74-W -177', '74-W -176', ' 74-W -175',
            #    '73-Ta-180', '73-Ta-179', '73-Ta-178', '73-Ta-177', '73-Ta-176', '73-Ta-175',
            #    '72-Hf-175', '72-Hf-173',
            ]
    # Target isotope - needed if you have multiple sub-directories for each target isotope in an element
    target = 'Ta/181Ta'
    return products, target

def Cu(isotope):
    # #Cu - non isomers. Lacking: 

    products= [
            # '29-Cu- 64', '30-Zn- 60', '29-Cu- 57', '27-Co -61', '27-Co -60', '27-Co- 55', '28-Ni- 55', '29-Cu- 55', '26-Fe- 59',
            '29-Cu- 57',
            #'30-Zn- 65', '30-Zn- 63', '30-Zn- 62', '29-Cu- 61', '27-Co- 58', '27-Co- 56', # monitors
            #'29-Cu- 64', '29-Cu- 60',
            #'28-Ni- 57',
            #'27-Co- 61', '27-Co- 60', '27-Co- 55',
            #'26-Fe- 59',
            #'25-Mn- 56', '25-Mn- 54'
            ]
    if isotope == '63Cu':
        target = 'Cu/63Cu'
    elif isotope == '65Cu':
        target = 'Cu/65Cu'
    else:
        raise Exception("No valid target")
    return products, target

def Ni(isotope):
    # #Ni - non isomers. Lacking: 
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
    return products, target

# products, target = Ta()
# products, target = Cu('63Cu')
products, target = Cu('65Cu')
# products, target = Ni('58Ni')
# products, target = Ni('60Ni')
# products, target = Ni('61Ni')
# products, target = Ni('62Ni')
# products, target = Ni('64Ni')



path = './'+target+'/'
os.makedirs(path + 'plots', exist_ok=True)

with open('./'+target+'/defaults.lst', 'r') as f:



    # energies = [2.0]
    energies = []
    cross_sections = dict.fromkeys(products)
    reading = False

    i=1

    for line in f:
        # print(i)
        # # if ' Incident energy  ' in line and '(CMS)' not in line :
        if ' incident energy  ' in line and '(CMS)' not in line and 'C.M.' not in line:
            print(line)
            # current_energy = float(line[18:28].replace('D',"E"))
            current_energy = float(line[40:].replace('D',"E").replace("MeV (LAB)\n",""))
            # print(current_energy)
            if len(energies) == 0 or current_energy != energies[-1]:
                energies.append(current_energy)
        elif 'production cross section' in line and 'reac:' in line:
            isotope = line[2:11]
            print(isotope)
            if isotope == '23-V - 47':
                print(line)
            if isotope in cross_sections.keys():
                if cross_sections[isotope] is None:
                    cross_sections[isotope] = []
                while len(cross_sections[isotope])<len(energies)-1:
                    cross_sections[isotope].append(0.0)
                cross_sections[isotope].append(float(line[38:49]))
        i=i+1


# print(cross_sections)


for isotope in products:
    # print(isotope)
    if cross_sections[isotope] is None:
        print('nothing for ', isotope)
        continue
    #print('plotting ', isotope)
    try:
        # print(energies)
        plt.plot(energies,cross_sections[isotope])
        plt.xlabel('Proton Energy (MeV)')
        plt.ylabel('Cross Section (mb)')
        plt.title(isotope)
        plt.savefig(target+'/plots/'+isotope+"_empire.png")
        plt.close()
        # plt.show()
    except ValueError:
        continue
    # with open('.plots/'+isotope.replace(" ","")+"_empire.txt",'w') as f:
    with open('./'+target+'/plots/'+isotope.replace("-","")+"_empire.txt",'w') as f:
        for i in range(len(energies)):
            f.write(str(energies[i])+"\t"+str(cross_sections[isotope][i])+"\n")
