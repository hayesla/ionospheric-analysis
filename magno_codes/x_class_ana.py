from sunpy.time import parse_time
from sunpy import timeseries as ts
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib import dates
from pathlib import Path
import glob
from download_BIRR_data import search_data_and_download, find_files
import datetime
import numpy as np

goes_data_dir = "/Users/laurahayes/QPP/stats_study/TEBBS/goes_rawdata/"
euve_data_dir = "./euve_data/"
magno_data_dir = "./magno_files/"

x_flares = pd.read_csv("x_class_w_files.csv")

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


save_dir = "/Users/laurahayes/ionospheric_work/ionospheric-analysis/magno_codes/plots_x_flares/"
i = 0
for i in range(len(x_flares)):
	print("analyzing {:d}".format(i))
	try:
		file_euve = glob.glob(euve_data_dir + parse_time(x_flares.iloc[i]['event_date']).strftime("*%Y%m%d.txt"))[0]
		file_magno = glob.glob(magno_data_dir + parse_time(x_flares.iloc[i]['event_date']).strftime("*%Y%m%d*.txt"))[0]
		goes_file = glob.glob(goes_data_dir + parse_time(x_flares.iloc[i]['event_date']).strftime("go15%Y%m%d.fits"))[0]


		data_euve = euve_to_series(file_euve)
		data_goes = ts.TimeSeries(goes_file)
		data_mag_bx, data_mag_by, data_mag_bz = mag_to_series(file_magno)

		gl = data_goes.to_dataframe()['xrsb'].truncate(x_flares.iloc[i]['start_time'], x_flares.iloc[i]['end_time'])
		gs = data_goes.to_dataframe()['xrsa'].truncate(x_flares.iloc[i]['start_time'], x_flares.iloc[i]['end_time'])
		data_mag_bx = data_mag_bx.truncate(x_flares.iloc[i]['start_time'], x_flares.iloc[i]['end_time'])
		data_mag_by = data_mag_by.truncate(x_flares.iloc[i]['start_time'], x_flares.iloc[i]['end_time'])
		data_mag_bz = data_mag_bz.truncate(x_flares.iloc[i]['start_time'], x_flares.iloc[i]['end_time'])
		data_euve = data_euve.truncate(x_flares.iloc[i]['start_time'], x_flares.iloc[i]['end_time'])



		fig, ax = plt.subplots(2, sharex=True, figsize=(10, 8))

		ax1 = ax[0]
		ax2 = ax[1]

		ax1.plot(gl, color="r", label="1-8 $\mathrm{\AA}$")
		ax1.plot(gs, color="b", label="0.5-4 $\mathrm{\AA}$")	
		ax1.plot(np.nan, color='k', label=r'Ly$\alpha$')
		ax1.set_ylim(1e-9, 1e-3)
		ax1.set_yscale("log")
		ax1.tick_params(which="both", direction="in", right=True, top=True)
		ax1.set_ylabel("Flux (Wm$^{-2}$)")
		ax1.legend(loc="upper right")

		ax11 = ax1.twinx()
		ax11.plot(data_euve, color='k')
		ax11.set_ylabel("Flux (Wm$^{-2}$)")


		ax2.plot(data_mag_bx, label="Bx", color="k")
		ax2.plot(np.nan, color="grey", label="By")
		ax2.plot(np.nan, color="green", label="Bx")
		ax2.legend(loc="lower right")
		ax3 = ax2.twinx()
		ax3.plot(data_mag_by, label="By", color="grey")
		ax4 = ax2.twinx()
		ax4.plot(data_mag_bz, label="Bz", color="green")

		ax2.set_xlabel("Time {:s} UT".format(x_flares.iloc[i]['event_date']))
		ax2.set_xlim(parse_time(x_flares.iloc[i]['start_time']).datetime, parse_time(x_flares.iloc[i]['end_time']).datetime)
		# ax2.xaxis.set_major_locator(dates.HourLocator(interval=3))
		# ax2.xaxis.set_minor_locator(dates.HourLocator(interval=1))
		ax2.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
		ax2.tick_params(which="both", direction="in", right=True, top=True)

		ax1.axvline(parse_time(x_flares.iloc[i]['peak_time']).datetime, color="k", ls="dashed")
		ax2.axvline(parse_time(x_flares.iloc[i]['peak_time']).datetime, color="k", ls="dashed")

		ax1.grid()
		ax2.grid()
		plt.tight_layout()
		plt.subplots_adjust(hspace=0.05)
		plt.savefig(save_dir + "x_class_flares_{:d}_{:s}.png".format(i, x_flares.iloc[i]['start_time']), dpi=200)
		plt.close()
	except:
		print("{:d} didnt work! wah!".format(i))






