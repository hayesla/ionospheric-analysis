import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
from sunpy.time import parse_time
from sunpy import timeseries as ts
import datetime


filey = 'birr_mag_20131028_000001.txt'
magno = pd.read_csv(filey, delim_whitespace=True, skiprows=1,
					names=['Date','Time', 'Index', 'Bx', 'By', 'Bz', 'E1', 'E2', 'E3', 'E4', 'T(FG)', 'T(E)', 'volts'])

magno_time = [datetime.datetime.strptime(magno.iloc[i]['Date'] + ' ' +magno.iloc[i]['Time'], '%d/%m/%Y %H:%M:%S') for i in range(len(magno))]

bx = pd.Series(np.array(magno['Bx']), index=magno_time)
by = pd.Series(np.array(magno['By']), index=magno_time)
bz = pd.Series(np.array(magno['Bz']), index=magno_time)

h = np.sqrt(np.array(bx)**2 + np.array(by)**2)
H = pd.Series(h, index = magno_time)