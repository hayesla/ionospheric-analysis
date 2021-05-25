import pandas as pd 
import matplotlib.pyplot as plt 
plt.ion()
import numpy as np 
import glob 
from sunpy.time import parse_time 
from sunpy import timeseries as ts 

# sid_files
sid_file_dir = "/Users/laurahayes/ionospheric_work/vlf_data_all_birr/sid_alll/*%Y%m%d*NAA*"
# goes files
goes_file_dir = "/Users/laurahayes/QPP/stats_study/TEBBS/goes_rawdata/*%Y%m%d*.fits"


flares = pd.read_csv("flare_list.csv")
flares["start_time"] = pd.to_datetime(flares["start_time"])

def sid_to_series(file):

	sid = pd.read_csv(file, comment="#", names=["times", "data"])
	tt = parse_time(sid["times"]).datetime

	ser = pd.Series(sid["data"].values, index=tt)
	return ser

def get_peak_flux(x):
	if x[0] == "X":
		return float(x[1:])*1e-4
	elif x[0] == "M":
		return float(x[1:])*1e-5
	elif x[0] == "C":
		return float(x[1:])*1e-6
	else:
		print("not what expected")

flares["goes_flux"] = flares["goes_class"].map(get_peak_flux)


results = []
for i in range(len(flares)):
	print(i)
	try:
		sid_file = glob.glob(flares.iloc[i]["start_time"].strftime(sid_file_dir))[0]
		sid_data = sid_to_series(sid_file).truncate(flares.iloc[i]["start_time"], flares.iloc[i]["end_time"])

		sid_data = sid_data.resample("5S").mean()

		p_vlf = np.max(sid_data)
		p_vlf2 = np.abs(np.max(sid_data) - sid_data[0])
		p_xray = flares.iloc[i]["goes_flux"]
		dt_value = (sid_data.index[np.argmax(sid_data)] - parse_time(flares.iloc[i]["peak_time"]).datetime).total_seconds()

		event = {}
		event["start_time_goes"] = flares.iloc[i]["start_time"]
		event["peak_flare"] = p_xray
		event["max_vlf"] = p_vlf
		event["abs_vlf"] = p_vlf2
		event["dt_value"] = dt_value

		results.append(event)
	except:
		print("issue with {:d}".format(i))

results = pd.DataFrame(results)

fig, ax = plt.subplots()
ax.scatter(results["peak_flare"], results["abs_vlf"], marker='.', c=np.log10(results["peak_flare"]), cmap="magma_r")
ax.set_xscale("log")
ax.set_xlim(1e-6, 1e-3)
ax.set_ylim(0, 5)

ax.set_xlabel("GOES peak 1-8$\mathrm{\AA}$ (Wm$^{-2}$)")
ax.set_ylabel("VLF absolute increase (volts)")
ax.tick_params(which="both", direction="in")
plt.tight_layout()


def check_avail():
	for i in range(len(flares)):
		goes_file = glob.glob(flares.iloc[i]["start_time"].strftime(goes_file_dir))
		sid_file = glob.glob(flares.iloc[i]["start_time"].strftime(sid_file_dir))
		if len(goes_file)==0:
			print("no goes file", i)
		if len(sid_file)==0:
			print("no sid_file", i)
