import os 

currentPath = os.getcwd()

myGeneratedFiles = currentPath + '/generatedfiles'
calibrationFiles = myGeneratedFiles + '/calibrationfiles'
experimentSpectra = 'experimentspectra'
calibrationSpectra = 'calibrationspectra'
peakData = myGeneratedFiles + '/peakdata'
peakData_data = myGeneratedFiles + '/peakdata/data' 
peakData_figures = myGeneratedFiles + '/peakdata/figures' 


def checkOrMakeDirectory(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
        print(f"Directory '{dir}' was successfully created")
    else:
        print(f"Directory '{dir}' already exists.")



# checkOrMakeDirectory(myGeneratedFiles)
# checkOrMakeDirectory(calibrationFiles)
# checkOrMakeDirectory(peakData)

# checkOrMakeDirectory(peakData_data)
# checkOrMakeDirectory(peakData_figures)

