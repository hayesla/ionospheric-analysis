from sunpy.time import parse_time
import pandas as pd 
import matplotlib.pyplot as plt 
from pathlib import Path
from download_BIRR_data import search_data_and_download, find_files
import sys
sys.path.append('..')
from goes_event_list import get_goes_event_list 


def get_flarelist():
	"""
	getting flare list that overlaps with the Milligan et al. 2020 paper. 
	Only looking for flares > M5 goes list.
	""" 
	t_start = '2011-02-01 00:00'
	t_end = '2016-06-06 00:00'
	get_goes_event_list(t_start, t_end, filename=Path.cwd().joinpath('goes_flare_list_m5.csv'), goes_class_filter='M5.0')


goes_flares = pd.read_csv('goes_flare_list_m5.csv')
goes_flares['peak_times_hours'] = [x.hour for x in pd.to_datetime(goes_flares['peak_time'])]
daytime_flares = goes_flares[(goes_flares['peak_times_hours']>8) & (goes_flares['peak_times_hours']<20)]

events_to_download = daytime_flares['event_date'].unique()

logg = []
for i in range(len(events_to_download)):
	print(i)
	try:
		log = search_data_and_download(events_to_download[i], path=Path.cwd().joinpath('magno_files'))
	except:
		log = '404 err probs'


	logg.append(log)