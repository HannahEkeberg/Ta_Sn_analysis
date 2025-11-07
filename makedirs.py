import os 

currentPath = os.getcwd()

myGeneratedFiles = currentPath + '/generatedFiles'
calibrationFiles = myGeneratedFiles + '/calibrationFiles'
experimentSpectra = 'experimentSpectra'
calibrationSpectra = 'calibrationSpectra'
peakData = myGeneratedFiles + '/peakData'

def checkOrMakeDirectory(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
        print(f"Directory '{dir}' was successfully created")
    else:
        print(f"Directory '{dir}' already exists.")


checkOrMakeDirectory(myGeneratedFiles)
checkOrMakeDirectory(calibrationFiles)
checkOrMakeDirectory(peakData)
