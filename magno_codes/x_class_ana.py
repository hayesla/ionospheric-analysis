from sunpy.time import parse_time
from sunpy import timeseries as ts
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib import dates
from pathlib import Path
import glob
from download_BIRR_data import search_data_and_download, find_files
import datetime

goes_data_dir = "/Users/laurahayes/QPP/stats_study/TEBBS/goes_rawdata/"
euve_data_dir = "./euve_data/"
magno_data_dir = "./magno_files/"

x_flares = pd.read_csv("x_class_w_files.csv")

def euve_to_series(file):
	data = pd.read_csv(file, comment=';', names=['date', 'time', 'count', 'flag'], delim_whitespace=True)
	tt = data['date'] + ' ' +  data['time']
	t_index = [parse_time(x) for x in tt]
	return pd.Series(data['count'].value, index=t_index)


i = 0

file_euve = glob.glob(euve_data_dir + parse_time(x_flares.iloc[i]['event_date']).strftime("*%Y%m%d.txt"))[0]
file_magno = glob.glob(magno_data_dir + parse_time(x_flares.iloc[i]['event_date']).strftime("*%Y%m%d*.txt"))[0]
goes_file = glob.glob(goes_data_dir + parse_time(x_flares.iloc[i]['event_date']).strftime("go15%Y%m%d.fits"))[0]


data_euve = pd.read_csv(file_euve, comment=';', names=['date', 'time', 'count', 'flag'], delim_whitespace=True)