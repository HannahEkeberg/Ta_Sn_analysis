import pandas as pd
import matplotlib.pyplot as plt

irradiation_55MeV = './experimental_current/stack_55.csv'
irradiation_30MeV = './experimental_current/stack_30.csv'

df1 = pd.read_csv(irradiation_30MeV)
current = df1['Current (nA)'].values
timepoints = df1['Time'].values
plt.plot(timepoints, current)
plt.xlabel('')
plt.show()