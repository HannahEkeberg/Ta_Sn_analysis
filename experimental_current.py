import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

irradiation_55MeV = './experimental_current/stack_55.csv'
irradiation_30MeV = './experimental_current/stack_30.csv'

df_30 = pd.read_csv(irradiation_30MeV)
df_55 = pd.read_csv(irradiation_55MeV)

def plot(df):
    current = df['Current (nA)'].values
    timepoints = df['Time'].values
    threshold = 20
    mask = current > threshold

    # if not np.any(mask):
    #     raise ValueError("Ingen punkter over terskel – sjekk threshold")

    start_idx = np.argmax(mask)
    end_idx = len(mask) - np.argmax(mask[::-1])
    buffer=50
    end_idx+=buffer

    current_trimmed = current[start_idx:end_idx]
    time_trimmed = timepoints[start_idx:end_idx]
    plt.plot(time_trimmed, current_trimmed)

    # plt.plot(timepoints, current)
    plt.locator_params(axis='x', nbins=5) 
    # plt.xlabel('')
    plt.show()

plot(df_55)