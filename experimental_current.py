import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

irradiation_55MeV = './experimental_current/stack_55.csv'
irradiation_30MeV = './experimental_current/stack_30.csv'

df1 = pd.read_csv(irradiation_55MeV)
df1 = pd.read_csv(irradiation_30MeV)
current = df1['Current (nA)'].values
timepoints = df1['Time'].values


threshold = 20

mask = current > threshold

# if not np.any(mask):
#     raise ValueError("Ingen punkter over terskel – sjekk threshold")

start_idx = np.argmax(mask)
end_idx = len(mask) - np.argmax(mask[::-1])

current_trimmed = current[start_idx:end_idx]
time_trimmed = timepoints[start_idx:end_idx]
plt.plot(time_trimmed, current_trimmed)

# plt.plot(timepoints, current)
plt.locator_params(axis='x', nbins=5) 
# plt.xlabel('')
plt.show()