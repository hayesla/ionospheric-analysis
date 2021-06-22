import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 
import glob 
from sunpy.time import parse_time 
from sunpy import timeseries as ts 
import datetime
from matplotlib import dates

# sid_files and GOES files
sid_file_dir = "/Users/laurahayes/ionospheric_work/vlf_data_all_birr/sid_alll/*%Y%m%d*NAA*"
all_sid = glob.glob("/Users/laurahayes/ionospheric_work/vlf_data_all_birr/sid_alll/*NAA*"); all_sid.sort()
goes_file_dir = "/Users/laurahayes/QPP/stats_study/TEBBS/goes_rawdata/*%Y%m%d*.fits"

flare_list = pd.read_csv("flare_list_for_vlf.csv")
flare_list["event_starttime"] = pd.to_datetime(flare_list["event_starttime"])

def sid_to_series(file):

	sid = pd.read_csv(file, comment="#", names=["times", "data"])
	tt = parse_time(sid["times"]).datetime
	ser = pd.Series(sid["data"].values, index=tt)
	ser.sort_index(inplace=True)
	return ser


unique_flare_days = flare_list["unique_day"].unique()

i = 0
	# check if sid file for day exists, and if so get it
def plot_all(plot=False):	
	ind_file = []
	too_many = []
	for i in range(len(unique_flare_days)):
		print(i)

		sid_file = glob.glob(pd.to_datetime(unique_flare_days[i]).strftime(sid_file_dir))
		if len(sid_file)==0:
			print("no file for {:d}".format(i))
			continue
		elif len(sid_file)>1:
			print("too many files for {:d}".format(i))
			too_many.append(i)
			continue
		else:
			print("whoop a file for {:d}".format(i))
			sid_data = sid_to_series(sid_file[0]).truncate(unique_flare_days[i] + " 07:00", 
														   unique_flare_days[i] + " 21:00")
			ind_file.append(i)
		if plot:
			flares_for_day = flare_list[flare_list["unique_day"].isin([unique_flare_days[i]])]
			fig, ax = plt.subplots()
			ax.plot(sid_data)
			ax.plot(sid_data.resample("30s").mean(), color="grey")
			ax.set_ylim(-5, 5)
			ax.set_xlim(pd.to_datetime(unique_flare_days[i] + " 07:00"), pd.to_datetime(unique_flare_days[i] + " 21:00"))
			for j in range(len(flares_for_day)):
				plt.axvline(pd.to_datetime(flares_for_day.iloc[j]["event_peaktime"]), lw=0.5, color="grey")
				plt.text(pd.to_datetime(flares_for_day.iloc[j]["event_peaktime"]), j/2, flares_for_day.iloc[j]["goes_class"])
			ax.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
			ax.set_xlabel("Time {:s}".format(unique_flare_days[i]))



			plt.tight_layout()
			plt.savefig("./test_plots/flare_{:s}.png".format(unique_flare_days[i]))
			plt.close()
	return ind_file, too_many

def do_things():
	vals, too_many = plot_all()
	unique_flare_with_data = unique_flare_days[vals]
	new_df = flare_list[flare_list["unique_day"].isin([unique_flare_with_data[0]])]
	for i in range(1, len(unique_flare_with_data)):
		new_df = new_df.append(flare_list[flare_list["unique_day"].isin([unique_flare_with_data[i]])])

	new_df.reset_index(inplace=True)


def second_runthrough():
	gg = glob.glob("/Users/laurahayes/ionospheric_work/ionospheric-analysis/stats_study/plots_that_work/*.png")
	gg.sort()
	flare_unique_day = [x.split("/")[-1][6:16] for x in gg]

	new_df = flare_list[flare_list["unique_day"].isin([flare_unique_day[0]])]
	for i in range(1, len(flare_unique_day)):
		new_df = new_df.append(flare_list[flare_list["unique_day"].isin([flare_unique_day[i]])])

	new_df.reset_index(inplace=True)





	for i in range(len(new_df)):
		print(i)
		sid_file = glob.glob(pd.to_datetime(new_df["event_starttime"].iloc[i]).strftime(sid_file_dir))[0]

		new_ts = pd.to_datetime(new_df["event_starttime"].iloc[i])-datetime.timedelta(minutes=10)
		new_te = pd.to_datetime(new_df["event_endtime"].iloc[i])+datetime.timedelta(minutes=10)

		sid_data = sid_to_series(sid_file).truncate(new_ts, new_te)

		goes_file = glob.glob(pd.to_datetime(new_df["event_starttime"].iloc[i]).strftime(goes_file_dir))[0]
		goes = ts.TimeSeries(goes_file).truncate(new_ts, new_te)
		gl = goes.to_dataframe()["xrsb"]


		fig, ax = plt.subplots(2, sharex=True)
		ax[0].plot(gl)
		ax[1].plot(sid_data)

		for a in ax:
			a.axvline(pd.to_datetime(new_df["event_peaktime"].iloc[i]))
			a.axvline(pd.to_datetime(new_df["event_starttime"].iloc[i]))
			a.axvline(pd.to_datetime(new_df["event_endtime"].iloc[i]))

		ax[0].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
		tstart_str = pd.to_datetime(new_df["event_starttime"].iloc[i]).strftime("%Y-%m-%dT%H:%M")
		ax[1].set_xlabel(new_df["event_peaktime"].iloc[i])
		plt.tight_layout()

		plt.savefig("./all_flares/flare_{:s}.png".format(tstart_str))
		plt.close()


def final_runthrough():
	gg = glob.glob("/Users/laurahayes/ionospheric_work/ionospheric-analysis/stats_study/final_flares/*.png")
	gg.sort()
	flare_starttimes = [x.split("/")[-1][6:22] for x in gg]
	matchtimes = [pd.to_datetime(x).strftime("%Y-%m-%d %H:%M:%S") for x in flare_starttimes]

	final_df = flare_list[flare_list["event_starttime"].isin([matchtimes[0]])]
	for i in range(1, len(matchtimes)):
		final_df = final_df.append(flare_list[flare_list["event_starttime"].isin([matchtimes[i]])])

	final_df.reset_index(inplace=True)
	final_df.to_csv("final_paper_vlf_flares2.csv", index_label=False)