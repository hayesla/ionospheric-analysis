import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 
import glob 
from sunpy.time import parse_time 
from sunpy import timeseries as ts 
import datetime
from matplotlib import dates
from scipy.signal import savgol_filter
plt.ion()

vlf_flares = pd.read_csv("final_paper_vlf_flares.csv")
vlf_flares["event_starttime"] = pd.to_datetime(vlf_flares["event_starttime"])

def calc_amp(x):
	return 20*np.log10(x + 5) - 61 + 107

def sid_to_series(file, amp=False):

	sid = pd.read_csv(file, comment="#", names=["times", "data"])
	tt = parse_time(sid["times"]).datetime
	if amp:
		ser = pd.Series(calc_amp(sid["data"].values), index=tt)
	else:
		ser = pd.Series(sid["data"].values, index=tt)		
	ser.sort_index(inplace=True)
	return ser

# File patterns
sid_file_dir = "/Users/laurahayes/ionospheric_work/vlf_data_all_birr/sid_alll/*%Y%m%d*NAA*"
goes_file_dir = "/Users/laurahayes/QPP/stats_study/TEBBS/goes_rawdata/*%Y%m%d*.fits"


errors = []
results = []
save=True
for i in range(len(vlf_flares)):
	print(i)
	try:
		sid_file = glob.glob(vlf_flares.iloc[i]["event_starttime"].strftime(sid_file_dir))[0]
		new_ts = pd.to_datetime(vlf_flares["event_starttime"].iloc[i])-datetime.timedelta(minutes=5)
		new_te = pd.to_datetime(vlf_flares["event_endtime"].iloc[i])+datetime.timedelta(minutes=5)
		sid_data = sid_to_series(sid_file).truncate(new_ts, new_te)
	

		goes_file = glob.glob(pd.to_datetime(vlf_flares["event_starttime"].iloc[i]).strftime(goes_file_dir))[0]
		goes = ts.TimeSeries(goes_file).truncate(new_ts, new_te)
		gl = goes.to_dataframe()["xrsb"]
		gs = goes.to_dataframe()["xrsa"]
		gl_flare = gl.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])
		gs_flare = gs.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])
		
		window_sec =  (sid_data.index[1] - sid_data.index[0]).total_seconds()
		window = int((2*60)/window_sec)
		if window%2 ==0:
			window = window+1


		sid_resample = pd.Series(savgol_filter(sid_data, int(window), 3), index=sid_data.index)
		sid_resample_flare = sid_resample.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])
		p_vlf = np.max(sid_resample_flare)
		p_vlf2 = np.abs(np.max(sid_resample_flare) - sid_resample_flare[0])
		p_xray = vlf_flares["goes_class_val"].iloc[i]
		dt_value_gl = (sid_resample_flare.index[np.argmax(sid_resample_flare)] - gl_flare.index[np.argmax(gl_flare)]).total_seconds()
		dt_value_gs = (sid_resample_flare.index[np.argmax(sid_resample_flare)] - gs_flare.index[np.argmax(gs_flare)]).total_seconds()

		event = {}
		event["start_time_goes"] = vlf_flares.iloc[i]["event_starttime"]
		event["peak_flare_gl"] = np.max(gl_flare)
		event["peak_flare_gs"] = np.max(gs_flare)
		event["max_vlf"] = p_vlf
		event["abs_vlf"] = p_vlf2
		event["dt_value_gl"] = dt_value_gl
		event["dt_value_gs"] = dt_value_gs


		results.append(event)
	
		fig, ax = plt.subplots(2, sharex=True)
		ax[0].plot(gl, color="r", label="1-8$\mathrm{\AA}$")
		ax[0].set_ylabel("Flux (Wm$^{-2}$)")
		ax[0].legend(loc="upper left")
		ax[0].set_yscale("log")

		ax[1].plot(sid_data, color="grey", lw=0.5, label="raw")
		ax[1].plot(sid_resample, color="k", label="2 minute smooth")
		ax[1].legend(loc="upper left")

		for a in ax:
			a.axvline(gl_flare.index[np.argmax(gl_flare)], color="r", lw=0.8, ls="dashed")
			a.axvline(gs_flare.index[np.argmax(gs_flare)], color="b", lw=0.8, ls="dashed")

			a.axvline(pd.to_datetime(vlf_flares["event_starttime"].iloc[i]), color="k", lw=0.8, ls="dashed" )
			a.axvline(pd.to_datetime(vlf_flares["event_endtime"].iloc[i]), color="k", lw=0.8, ls="dashed")

			a.axvline(sid_resample_flare.index[np.argmax(sid_resample_flare)], color="g")

		ax[0].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
		tstart_str = pd.to_datetime(vlf_flares["event_starttime"].iloc[i]).strftime("%Y-%m-%dT%H:%M")
		ax[1].set_xlabel(vlf_flares["event_peaktime"].iloc[i])
		plt.tight_layout()
		plt.subplots_adjust(hspace=0.01)
		if save:
			plt.savefig("./final_stats_study_tests/flare2_{:d}_{:s}.png".format(i, tstart_str))
			plt.close()


	except:
		errors.append(i)
		print("issue with {:d}".format(i))

results = pd.DataFrame(results)

def plot_flare(new_df, i, save=False):
	"""
	Function to plot a flare from a pandas DataFrame.

	Parameters:
	----------
	new_df : ~`pd.DataFrame
		DataFrame with each row a flare
	i : ~`int`
		row index to plot

	"""
	sid_file = glob.glob(pd.to_datetime(new_df["event_starttime"].iloc[i]).strftime(sid_file_dir))[0]

	new_ts = pd.to_datetime(new_df["event_starttime"].iloc[i])-datetime.timedelta(minutes=5)
	new_te = pd.to_datetime(new_df["event_endtime"].iloc[i])+datetime.timedelta(minutes=5)

	sid_data = sid_to_series(sid_file).truncate(new_ts, new_te)

	sid_resample = pd.Series(savgol_filter(sid_data, 2*60+1, 3), index=sid_data.index)
	if len(sid_data)>300:
		sid_resample2 = pd.Series(savgol_filter(sid_data, 5*60+1, 3), index=sid_data.index)
	tmax_sid = sid_resample.index[np.argmax(sid_resample)]


	goes_file = glob.glob(pd.to_datetime(new_df["event_starttime"].iloc[i]).strftime(goes_file_dir))[0]
	goes = ts.TimeSeries(goes_file).truncate(new_ts, new_te)
	gl = goes.to_dataframe()["xrsb"]
	gs = goes.to_dataframe()["xrsa"]

	fig, ax = plt.subplots(2, sharex=True)
	ax[0].plot(gl, color="r", label="1-8$\mathrm{\AA}$")
	ax[0].plot(gs, color="b", label="0.5-4$\mathrm{\AA}$")
	ax[0].set_ylabel("Flux (Wm$^{-2}$)")
	ax[0].legend(loc="upper left")
	ax[0].set_yscale("log")

	ax[1].plot(sid_data, color="grey", lw=0.5, label="raw")
	ax[1].plot(sid_resample, color="k", label="2 minute smooth")
	if len(sid_data)>300:
		ax[1].plot(sid_resample2, color="g", label="5 minute smooth")	
	ax[1].legend(loc="upper left")

	for a in ax:
		a.axvline(pd.to_datetime(new_df["event_peaktime"].iloc[i]), color="k", lw=0.8, ls="dashed")
		a.axvline(pd.to_datetime(new_df["event_starttime"].iloc[i]), color="k", lw=0.8, ls="dashed" )
		a.axvline(pd.to_datetime(new_df["event_endtime"].iloc[i]), color="k", lw=0.8, ls="dashed")
		a.axvline(tmax_sid, color="r")
	ax[0].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
	tstart_str = pd.to_datetime(new_df["event_starttime"].iloc[i]).strftime("%Y-%m-%dT%H:%M")
	ax[1].set_xlabel(new_df["event_peaktime"].iloc[i])
	plt.tight_layout()
	plt.subplots_adjust(hspace=0.01)
	if save:
		plt.savefig("./final_stats_study_tests/flare2_{:d}_{:s}.png".format(i, tstart_str))
		plt.close()

# for i in range(len(vlf_flares)):
# 	plot_flare(vlf_flares, i, save=True)
# 	print(i)