import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
import glob
import datetime
from matplotlib import dates
from astropy import units as u
import time 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sunpy.time import parse_time
# all the files.
files = glob.glob("/Users/laurahayes/ionospheric_work/vlf_data_all_birr/sid_alll/*.csv")
files.sort()

# get all the dates from the files.
dates_files = [datetime.datetime.strptime(f.split("/")[-1][0:8], "%Y%m%d") for f in files]
# get unique dates - some days have multiple files.
unique_dates = list(set(dates_files))
unique_dates.sort()

# find all the days between the start and end date of available days to 
# allow us to populate nans for days with no data.
all_days = []
day1 = unique_dates[0]
all_days.append(day1)
while day1 <= unique_dates[-1]:
	day1 = day1 + datetime.timedelta(days=1)
	all_days.append(day1)


def find_files_by_date(date):
	"""
	Find the files available for a particular date.
	Some days have multiple files that need to be 
	concatenated - this returns a list of files
	for that date. 

	Parameters
	----------
	date : ~ datetime.datetime

	Returns
	-------
	file_list : `list`

	"""
	date_str = date.strftime('%Y%m%d')
	file_list = []
	for f in files:
		if date_str in f:
			file_list.append(f)
	return file_list

def read_files(files):
	"""
	Reads the file(s) in the list `files` and returns
	a pandas dataframe. If the list is > 1 then this
	function concatenates the list data into one dataframe.

	Parameters
	----------
	files : `list`
		a list contains the filepath to the file(s)

	Returns
	-------
	pd.DataFrame

	"""
	if len(files) == 1:
		return pd.read_csv(files[0], comment='#', names=["time", "volts"])

	elif len(files)>1:
		df = []
		for f in files:
			data = pd.read_csv(f, comment='#', names=["time", "volts"])
			df.append(data)
		new_df = pd.concat(df)
		new_df = new_df.drop_duplicates(subset='time')
		new_df.reset_index(drop=True, inplace=True)
		return new_df

def make_empty_df(date):
	"""
	Returns an dataframe of all times of the day for a particular 
	date with nans for the 'volts'.

	Parameters
	----------
	date : `datetime.datetime`

	Returns
	-------
	pd.DataFrame

	"""
	dti = pd.date_range(date.strftime('%Y-%m-%d'), periods=17280, freq='5S')
	empty_arr = np.empty(len(dti))
	empty_arr[:] = np.nan
	df = pd.DataFrame({'time':np.array(dti.astype('str')), 'volts':empty_arr})
	return df


def make_all(date):

	empty_df = make_empty_df(date)

	files = find_files_by_date(date)
	if len(files) == 0:
		return empty_df
	else:
		df_data = read_files(files)
		df_merge = pd.merge(empty_df, df_data, on='time', how='left').drop_duplicates('time')
		return df_merge.rename(columns={'volts_y':'volts'}).fillna(method='ffill')


def make_full_arr():
	t1 = time.time()
	full_arr = []
	for i in range(0, len(all_days)):
		full_arr.append(make_all(all_days[i])['volts'].values)
	full_arr = np.array(full_arr)
	t2 = time.time()
	print(t2 - t1)

	full_array = np.stack(np.array(full_arr))
	np.save('full_array_data_birr', full_array)

def plot_all_data(cmap=plt.cm.viridis, filename='test_days.png'):

	plt.rcParams['font.family'] = 'Helvetica'
	full_array = np.load('full_array_data_birr.npy')
	#cmap = plt.cm.viridis
	cmap_colors = cmap(np.linspace(0,1,100))
	cmap.set_bad(cmap_colors[0])


	xaxis = parse_time(make_empty_df(all_days[0])['time']).datetime
	yaxis = all_days

	fig, ax = plt.subplots(figsize=(9, 9))
	im = ax.imshow(full_array, origin='lower', aspect='auto', cmap=cmap, 
			   extent=[dates.date2num(xaxis[0]), dates.date2num(xaxis[-1]),
			   		   dates.date2num(yaxis[0]), dates.date2num(yaxis[-1])])

	ax.axvline(xaxis[10081], color="w", ls="dashed", lw=0.5)
	ax.axvline(xaxis[11521], color="w", ls="dashed", lw=0.5)

	ax.xaxis_date()
	ax.yaxis_date()

	ax.xaxis.set_major_locator(dates.HourLocator(interval=3))
	ax.xaxis.set_minor_locator(dates.HourLocator(interval=1))

	ax.yaxis.set_major_locator(dates.MonthLocator(interval=6))
	ax.yaxis.set_minor_locator(dates.MonthLocator(interval=1))

	ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
	ax.yaxis.set_major_formatter(dates.DateFormatter('%Y-%b'))
	
	ax.set_xlabel('Time of day')
	ax.set_ylabel('Date')


	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="2%", pad=0.05)
	fig.add_axes(cax)
	cbar = fig.colorbar(im, cax=cax, orientation="vertical")
	cbar.set_label('Volts')

	plt.tight_layout()

	plt.savefig(filename, dpi=200)
	plt.close()

def plot_avg():
	cmap = plt.cm.viridis
	plt.rcParams['font.family'] = 'Helvetica'
	full_array = np.load('full_array_data_birr.npy')
	#cmap = plt.cm.viridis
	cmap_colors = cmap(np.linspace(0,1,100))
	cmap.set_bad(cmap_colors[0])


	xaxis = parse_time(make_empty_df(all_days[0])['time']).datetime
	yaxis = np.array(all_days)

	to_plot = np.mean(full_array[:, 10081:11521], axis=1)

	fig, ax = plt.subplots(figsize=(10, 5))
	#ax.scatter(yaxis[to_plot>-4.8], to_plot[to_plot>-4.8], marker='.', alpha=0.5)
	ax.scatter(yaxis, to_plot, marker='.', alpha=0.5, c=to_plot)
	ax.xaxis.set_major_locator(dates.MonthLocator(interval=6))
	ax.xaxis.set_minor_locator(dates.MonthLocator(interval=1))
	ax.set_xlabel("Time")
	ax.set_ylabel("VLF amplitude (volts)")
	ax.tick_params(which="both", direction="in")
	plt.tight_layout()
	plt.savefig("yearly_variations.png", dpi=200)


def omni_test():
	import heliopy.data.omni as omni
	cmap = plt.cm.viridis
	plt.rcParams['font.family'] = 'Helvetica'
	full_array = np.load('full_array_data_birr.npy')
	#cmap = plt.cm.viridis
	cmap_colors = cmap(np.linspace(0,1,100))
	cmap.set_bad(cmap_colors[0])


	xaxis = parse_time(make_empty_df(all_days[0])['time']).datetime
	yaxis = np.array(all_days)

	to_plot = np.mean(full_array[:, 10081:11521], axis=1)

	omni_data = pd.read_csv("omni_data.csv", parse_dates=True)
	omni_data_resample = omni_data.resample("1H").mean()
	# omni_data = omni.hro2_5min(yaxis[0], yaxis[-1])
	# omni_data = omni_data.to_dataframe()

	tstarty = pd.to_datetime("2015-01-01")
	tendy = pd.to_datetime("2016-03-01")

	ind_start = np.where(yaxis >= tstarty)[0][0]
	ind_end = np.where(yaxis <= tendy)[0][-1]


	omni_data = omni_data.truncate(tstarty, tendy)

	yaxis = yaxis[ind_start:ind_end]
	to_plot = to_plot[ind_start:ind_end]	

	fig, ax = plt.subplots(figsize=(10, 5))
	#ax.scatter(yaxis[to_plot>-4.8], to_plot[to_plot>-4.8], marker='.', alpha=0.5)
	ax.scatter(yaxis, to_plot, marker='.', alpha=0.8)#, c=to_plot)

	ax.set_ylim(-5, 3)
	ax.tick_params(axis="y", direction="in")
	ax.set_xlabel("Time (UT)")
	ax.set_ylabel("VLF amplitude (volts)")

	ax2 = ax.twinx()
	ax2.plot(omni_data["SYM_H"].resample("60T").mean(), lw=0.5, color="r")
	ax.set_zorder(ax2.get_zorder()+1)
	ax.patch.set_visible(False)
	ax2.set_ylim(50, -200)
	ax.spines['right'].set_color('red')
	ax2.tick_params(axis='y', colors='red', direction="in")
	ax2.yaxis.label.set_color('red')



def kp_index():

	kp_data = np.loadtxt("val_kindex.txt")
	t1 = datetime.datetime(2013,1,2,0,0,0)
	kp_times = [t1]
	for i in range(len(kp_data)):
		t1 = t1 + datetime.timedelta(hours=3)
		kp_times.append(t1)

	kp_df = pd.Series(kp_data, index=kp_times[:-1])
	kp_df = kp_df.truncate("2013-01-02", "2018-11-27")
	kp_df.replace([8, 9], np.nan, inplace=True)
	return kp_df


def test_kp_plot(kp=True):

	kp_df = kp_index()

	cmap = plt.cm.viridis
	plt.rcParams['font.family'] = 'Helvetica'
	full_array = np.load('full_array_data_birr.npy')
	#cmap = plt.cm.viridis
	cmap_colors = cmap(np.linspace(0,1,100))
	cmap.set_bad(cmap_colors[0])


	xaxis = parse_time(make_empty_df(all_days[0])['time']).datetime
	yaxis = np.array(all_days)

	to_plot = np.mean(full_array[:, 10081:11521], axis=1)

	fig, ax = plt.subplots(figsize=(10, 5))
	#ax.scatter(yaxis[to_plot>-4.8], to_plot[to_plot>-4.8], marker='.', alpha=0.5)
	ax.scatter(yaxis, to_plot, marker='.', alpha=0.5, c=to_plot)

	ax.set_ylim(-5, 3)
	ax.tick_params(axis="y", direction="in")
	ax.set_xlabel("Time (UT)")
	ax.set_ylabel("VLF amplitude (volts)")

	if kp:
		ax2 = ax.twinx()
		ax2.plot(kp_df[kp_df>3], color='r', lw=0.5, drawstyle="steps-mid", label="KP index > 3")
		ax.set_zorder(ax2.get_zorder()+1)
		ax.patch.set_visible(False)
		ax2.set_ylabel("KP index")
		ax.spines['right'].set_color('red')
		ax2.tick_params(axis='y', colors='red', direction="in")
		ax2.yaxis.label.set_color('red')
		ax2.set_ylim(4, 7.1)

	ax.set_xlim(yaxis[0], yaxis[-1])
	ax.xaxis.set_major_locator(dates.MonthLocator(interval=6))
	ax.xaxis.set_minor_locator(dates.MonthLocator(interval=1))


	plt.tight_layout()
	if kp:
		plt.savefig("yearly_variations_kp.png", dpi=200)
		plt.close()
	else:
		plt.savefig("yearly_variations_nokp.png", dpi=200)		
		plt.close()

def test_kp_plot_zoom(kp=True):

	kp_df = kp_index()

	cmap = plt.cm.viridis
	plt.rcParams['font.family'] = 'Helvetica'
	full_array = np.load('full_array_data_birr.npy')
	#cmap = plt.cm.viridis
	cmap_colors = cmap(np.linspace(0,1,100))
	cmap.set_bad(cmap_colors[0])


	xaxis = parse_time(make_empty_df(all_days[0])['time']).datetime
	yaxis = np.array(all_days)

	to_plot = np.mean(full_array[:, 10081:11521], axis=1)

	tstarty = pd.to_datetime("2015-01-01")
	tendy = pd.to_datetime("2016-03-01")

	ind_start = np.where(yaxis >= tstarty)[0][0]
	ind_end = np.where(yaxis <= tendy)[0][-1]


	kp_df = kp_df.truncate(tstarty, tendy)

	yaxis = yaxis[ind_start:ind_end]
	to_plot = to_plot[ind_start:ind_end]


	fig, ax = plt.subplots()
	#ax.scatter(yaxis[to_plot>-4.8], to_plot[to_plot>-4.8], marker='.', alpha=0.5)
	ax.scatter(yaxis, to_plot, marker='.', alpha=0.8)

	ax.set_ylim(-5, 3)
	ax.tick_params(axis="y", direction="in")
	ax.set_xlabel("Time (UT)")
	ax.set_ylabel("VLF amplitude (volts)")

	if kp:
		ax2 = ax.twinx()
		ax2.plot(kp_df[kp_df>3], color='r', lw=0.7, drawstyle="steps-mid", label="KP index > 3")
		ax.set_zorder(ax2.get_zorder()+1)
		ax.patch.set_visible(False)
		ax2.set_ylabel("KP index")
		ax.spines['right'].set_color('red')
		ax2.tick_params(axis='y', colors='red', direction="in")
		ax2.yaxis.label.set_color('red')
		ax2.set_ylim(4, 7.1)

	ax.set_xlim(yaxis[0], yaxis[-1])
	ax.xaxis.set_major_locator(dates.MonthLocator(interval=2))
	ax.xaxis.set_minor_locator(dates.MonthLocator(interval=1))

	
	plt.tight_layout()
	if kp:
		plt.savefig("2015_variations_kp.png", dpi=200)
		plt.close()
	else:
		plt.savefig("2015_variations_nokp.png", dpi=200)		
		plt.close()
