import matplotlib.pyplot as plt 
import numpy as np 
from sunpy.time import parse_time
import datetime
from sunpy import timeseries as ts 
import pandas as pd 


def euve_to_series(file):
	data = pd.read_csv(file, comment=';', names=['date', 'time', 'count', 'flag'], delim_whitespace=True)
	tt = data['date'] + ' ' +  data['time']
	t_index = [datetime.datetime.strptime(x, '%d-%b-%Y %H:%M:%S.%f') for x in tt]
	ser = pd.Series(data['count'].values, index=t_index)
	ser = ser.replace(-99999.0, np.nan)
	return ser

def mag_to_series(file):

	magno = pd.read_csv(file, delim_whitespace=True, skiprows=1,
						names=["Date","Time", "Index", "Bx", "By", "Bz"], usecols=[0, 1, 2, 3, 4, 5])


	magno_time = [datetime.datetime.strptime(magno.iloc[i]["Date"] + " " +magno.iloc[i]["Time"], "%d/%m/%Y %H:%M:%S") for i in range(len(magno))]

	bx = pd.Series(np.array(magno["Bx"]), index=magno_time)
	by = pd.Series(np.array(magno["By"]), index=magno_time)
	bz = pd.Series(np.array(magno["Bz"]), index=magno_time)
	return bx, by, bz

def sid_to_series(file):

	sid = pd.read_csv(file, comment="#", names=["times", "data"])
	tt = parse_time(sid["times"]).datetime

	ser = pd.Series(sid["data"].values, index=tt)
	return ser


