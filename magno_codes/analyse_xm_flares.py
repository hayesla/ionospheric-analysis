from sunpy.time import parse_time
from sunpy import timeseries as ts
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib import dates
from pathlib import Path
import glob
from download_BIRR_data import search_data_and_download, find_files
import datetime
import sys
sys.path.append('..')
from goes_event_list import get_goes_event_list 

goes_data_dir = '/Users/laurahayes/QPP/stats_study/TEBBS/goes_rawdata/'

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

def download_the_data():
	logg = []
	for i in range(len(events_to_download)):
		print(i)
		try:
			log = search_data_and_download(events_to_download[i], path=Path.cwd().joinpath('magno_files'))
		except:
			log = '404 err probs'


		logg.append(log)
		return log

files_downloaded = glob.glob('./magno_files/*txt')

# check which flares have downloaded files
index = []
for i in range(len(daytime_flares)):
	tt = parse_time(daytime_flares['event_date'].iloc[i]).strftime('%Y%m%d')
	ind = False
	for f in files_downloaded:
		if tt in f:
			ind=True
			continue
	index.append(ind)


flare_w_files = daytime_flares[index]

save_dir = '/Users/laurahayes/ionospheric_work/ionospheric-analysis/magno_codes/plots_magno/'
def plot_test(i):


	tt = parse_time(events_to_download[i]).strftime('%Y%m%d')
	files_magno = glob.glob('./magno_files/*{:s}*'.format(tt))
	if len(files_magno)==0:
		print('No magnetometer data')
		#break


	goes_file = goes_data_dir + 'go15' + tt + '.fits'
	if not Path(goes_file).exists():
		print('No goes data')
		#break

	goes_data = ts.TimeSeries(goes_file)
	gl = goes_data.data['xrsb']
	gs = goes_data.data['xrsa']
	flares_ind = np.where(daytime_flares['event_date'].isin([events_to_download[i]])==True)[0]
	flares = daytime_flares.iloc[flares_ind]

	filey = files_magno[0]

	magno = pd.read_csv(filey, delim_whitespace=True, skiprows=1,
						names=['Date','Time', 'Index', 'Bx', 'By', 'Bz', 'E1', 'E2', 'E3', 'E4', 'T(FG)', 'T(E)', 'volts'])

	magno_time = [datetime.datetime.strptime(magno.iloc[i]['Date'] + ' ' +magno.iloc[i]['Time'], '%d/%m/%Y %H:%M:%S') for i in range(len(magno))]

	bx = pd.Series(np.array(magno['Bx']), index=magno_time)
	by = pd.Series(np.array(magno['By']), index=magno_time)
	bz = pd.Series(np.array(magno['Bz']), index=magno_time)

	h = np.sqrt(np.array(bx)**2 + np.array(by)**2)
	H = pd.Series(h, index = magno_time)

	fig, ax = plt.subplots(2, sharex=True, figsize=(10, 8))

	ax1 = ax[0]
	ax2 = ax[1]

	ax1.plot(gl, color='r', label='1-8 $\mathrm{\AA}$')
	ax1.plot(gs, color='b', label='0.5-4 $\mathrm{\AA}$')	
	ax1.set_ylim(1e-9, 1e-3)
	ax1.set_yscale('log')
	ax1.tick_params(which='both', direction='in', right=True, top=True)
	ax1.set_ylabel('Flux (Wm$^{-2}$)')
	ax1.legend(loc='upper right')


	ax2.plot(magno_time, magno['Bx'], label='Bx', color='k')
	ax2.plot(np.nan, color='grey', label='By')
	ax2.plot(np.nan, color='green', label='Bx')
	ax2.legend(loc='lower right')
	ax3 = ax2.twinx()
	ax3.plot(magno_time, magno['By'], label='By', color='grey')
	ax4 = ax2.twinx()
	ax4.plot(magno_time, magno['Bz'], label='Bz', color='green')

	ax2.set_xlabel('Time {:s} UT'.format(events_to_download[i]))
	ax2.set_xlim(events_to_download[i] + ' 00:00', events_to_download[i] + ' 23:59')
	ax2.xaxis.set_major_locator(dates.HourLocator(interval=3))
	ax2.xaxis.set_minor_locator(dates.HourLocator(interval=1))
	ax2.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
	ax2.tick_params(which='both', direction='in', right=True, top=True)
	for f in flares['peak_time']:
		ax1.axvline(parse_time(f).datetime, color='k', ls='dashed')
		ax2.axvline(parse_time(f).datetime, color='k', ls='dashed')
	ax1.grid()
	ax2.grid()
	plt.tight_layout()
	plt.subplots_adjust(hspace=0.05)
	plt.savefig(save_dir + parse_time(events_to_download[i]).strftime('%Y%m%d.png'), dpi=200)
	plt.close()


# for i in range(len(events_to_download)):
# 	try:
# 		plot_test(i)
# 	except:
# 		print(i)


def save_x_class():
	goes_flares = pd.read_csv('goes_flare_list_m5.csv')
	goes_flares['class_short'] = [x[0] for x in goes_flares['goes_class']]
	goes_flares['peak_times_hours'] = [x.hour for x in pd.to_datetime(goes_flares['peak_time'])]

	daytime_flares = goes_flares[(goes_flares['peak_times_hours']>8) & (goes_flares['peak_times_hours']<20)]
	x_flares = daytime_flares[daytime_flares['class_short']=='X']

	files_downloaded = glob.glob('./magno_files/*txt')

	# check which flares have downloaded files
	index = []
	for i in range(len(x_flares)):
		tt = parse_time(x_flares['event_date'].iloc[i]).strftime('%Y%m%d')
		ind = False
		for f in files_downloaded:
			if tt in f:
				ind=True
				continue
		index.append(ind)


	final_flares = x_flares[index]
	final_flares.reset_index(drop=True, inplace=True)
	final_flares.to_csv('x_class_w_files.csv', index_label=True)
