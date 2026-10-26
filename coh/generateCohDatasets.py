from OrderCoh import *

# def generateCohDatasetsForIr():
#     targetIsotopes = ['191Ir', '193Ir']
#     noIsomers = ['078-188Pt', '078-189Pt', '078-191Pt', '077-188Ir']
#     groundstatesAndIsomers = ['078-193Pt', '077-189Ir', '077-190Ir', '077-192Ir', '077-194Ir']
#     for target in targetIsotopes:
#         coh = OrderCoh(target, 'Ir')
#         coh.createFolder()
#         coh.unifyFilesForParsing()
#         coh.readCoh(noIsomers)
#         coh.readCohIsomers(groundstatesAndIsomers)


# def generateCohDatasetsForFe():
#     FeTargetIsotopes = ['54Fe', '56Fe', '57Fe', '58Fe']
#     Fe_noIsomers = ['027-057Co', '027-056Co', '027-055Co','027-058Co', '026-059Fe', '025-056Mn', '025-054Mn', '024-051Cr', '023-048V','026-053Fe']
#     Fe_groundStatesAndIsomes = ['027-058Co', '025-052Mn', '026-053Fe']
#     for target in FeTargetIsotopes:
#         cohFe = OrderCoh(target, 'Fe')
#         cohFe.createFolder()
#         cohFe.unifyFilesForParsing()
#         cohFe.readCoh(Fe_noIsomers)
#         cohFe.readCohIsomers(Fe_groundStatesAndIsomes)

def generateCohDatasetsForNi():
    targetIsotopes = ['58Ni', '60Ni', '61Ni', '62Ni', '64Ni']
    noIsomers = ['029-064Cu','029-061Cu', '029-060Cu', '028-057Ni', '028-056Ni', '027-060Co', '027-057Co', '027-058Co', '027-056Co', '027-055Co', '026-059Fe', '025-054Mn', '025-052Mn', '025-051Mn', '024-051Cr']
    groundStatesAndIsomes = ['027-060Co', '027-058Co', '025-052Mn']
    for target in targetIsotopes:
        coh = OrderCoh(target, 'Ni')
        coh.createFolder()
        coh.unifyFilesForParsing()
        coh.readCoh(noIsomers)
        coh.readCohIsomers(groundStatesAndIsomes)

def generateCohDatasetsForCu():
    targetIsotopes = ['63Cu', '65Cu']
    # noIsomers = ['029-064Cu','029-061Cu','026-059Fe', '027-060Co', '027-061Co', '027-057Co', '027-056Co','028-065Ni', '030-065Zn', '030-063Zn', '030-062Zn', '025-054Mn', '025-056Mn', '029-060Cu', '029-057Cu', '028-055Ni', '028-055Cu', ]
    noIsomers = ['027-055Co']
    # groundStatesAndIsomes = ['028-065Ni']
    # groundStatesAndIsomes = ['027-060Co', '027-058Co', '025-052Mn']
    groundStatesAndIsomes = ['028-057Ni']
    # groundStatesAndIsomes = ['028-057Ni']
    for target in targetIsotopes:
        coh = OrderCoh(target, 'Cu')
        coh.createFolder()
        coh.unifyFilesForParsing()
        coh.readCoh(noIsomers)
        # coh.readCohIsomers(groundStatesAndIsomes)

def generateCohDatasetsForTa():
    target_isotopes = ['181Ta', '180Ta']
    groundstate_isotopes = ['074-181W', '074-178W', '074-177W', '074-176W', '074-175W', '073-179Ta', '073-177Ta', '073-176Ta', '073-175Ta', '073-174Ta', '073-173Ta', '072-175Hf', '072-173Hf']
    groundstates_and_isomer_isotopes = ['074-179W', '073-180Ta', '073-178Ta', '072-180Hf', '071-177Lu', '071-178Lu', '071-179Lu']
    for target in target_isotopes:
        coh = OrderCoh(target, 'Ta')
        coh.createFolder()
        coh.unifyFilesForParsing()
        coh.readCoh(groundstate_isotopes)
        coh.readCohIsomers(groundstates_and_isomer_isotopes)




# generateCohDatasetsForTa()
generateCohDatasetsForCu()
# generateCohDatasetsForNi()
