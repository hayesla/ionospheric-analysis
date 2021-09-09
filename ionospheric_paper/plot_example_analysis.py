import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 
import glob 
from sunpy.time import parse_time 
from sunpy import timeseries as ts 
import datetime
from matplotlib import dates
from scipy.signal import savgol_filter
import matplotlib
import seaborn as sns 
sns.set_context("paper", font_scale=1.2)

"""
Script for figure 3 of VLF paper
"""

matplotlib.rcParams['xtick.direction'] = "in"
matplotlib.rcParams['ytick.direction'] = "in"
matplotlib.rcParams['xtick.minor.visible'] = True
matplotlib.rcParams['ytick.minor.visible'] = True
plt.rcParams['font.family'] = 'Helvetica'


vlf_flares = pd.read_csv("final_paper_vlf_flares2.csv")
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

def plot_flare(i):
    new_ts = pd.to_datetime(vlf_flares["event_starttime"].iloc[i])-datetime.timedelta(minutes=5)
    new_te = pd.to_datetime(vlf_flares["event_endtime"].iloc[i])+datetime.timedelta(minutes=5)

    # SID data
    sid_file = glob.glob(vlf_flares.iloc[i]["event_starttime"].strftime(sid_file_dir))[0]
    sid_data = sid_to_series(sid_file).truncate(new_ts, new_te)
    sid_data_db = sid_to_series(sid_file, amp=True).truncate(new_ts, new_te)

    # smoothing window defined in terms of cadence
    window_sec =  (sid_data.index[1] - sid_data.index[0]).total_seconds()
    window = int((2*60)/window_sec)
    if window%2 == 0:
        window = window+1

    sid_resample = pd.Series(savgol_filter(sid_data, int(window), 3), index=sid_data.index)
    sid_resample_flare = sid_resample.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])
    sid_resample_db = pd.Series(savgol_filter(sid_data_db, int(window), 3), index=sid_data_db.index)
    sid_resample_flare_db = sid_resample_db.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])                
    # GOES data
    goes_file = glob.glob(pd.to_datetime(vlf_flares["event_starttime"].iloc[i]).strftime(goes_file_dir))[0]
    goes = ts.TimeSeries(goes_file).truncate(new_ts, new_te)
    gl = goes.to_dataframe()["xrsb"]
    gs = goes.to_dataframe()["xrsa"]
    gl_flare = gl.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])
    gs_flare = gs.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])

    fig, ax = plt.subplots(2, figsize=(8,6), sharex=True)
    ax[0].plot(gl, color="tab:red", label="1-8 $\mathrm{\AA}$")
    ax[0].plot(gs, color="tab:blue", label="0.5-4 $\mathrm{\AA}$")
    ax[0].set_ylabel("Flux [Wm$^{-2}$]")
    ax[0].legend(loc="upper left")
    ax[0].set_yscale("log")

    ax[1].plot(sid_data_db - sid_data_db[0], label="VLF amp", color="grey", lw=0.7)
    ax[1].plot(sid_resample_db - sid_data_db[0], label="Smoothed VLF amp", color="k")
    ax[1].legend(loc="upper left")      

    for a in ax:
        a.axvline(gl_flare.index[np.argmax(gl_flare)], color="tab:red", lw=0.4)
        a.axvline(gs_flare.index[np.argmax(gs_flare)], color="tab:blue", lw=0.4)
        a.axvline(pd.to_datetime(vlf_flares["event_starttime"].iloc[i]), ls="dashed", color="grey")
        a.axvline(pd.to_datetime(vlf_flares["event_endtime"].iloc[i]), ls="dashed", color="grey")
        a.axvline(sid_resample_flare.index[np.argmax(sid_resample_flare)], color="k", lw=0.4)

    tstart_str = pd.to_datetime(vlf_flares["event_starttime"].iloc[i]).strftime("%Y-%m-%dT%H:%M")
    ax[1].set_xlabel("Time {:s} UT".format(pd.to_datetime(vlf_flares["event_starttime"].iloc[i]).strftime("%Y-%m-%d %H:%M")))
    ax[1].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
    ax[1].set_xlim(gl.index[0], gl.index[-1])
    ax[1].tick_params(which="both", direction="in")
    ax[0].tick_params(which="both", direction="in")
    ax[1].set_ylabel("VLF amplitude excess [db]")
    plt.tight_layout()
    ax[1].xaxis.set_major_locator(dates.MinuteLocator(byminute=[35, 40, 45, 50]))
    ax[1].xaxis.set_minor_locator(dates.MinuteLocator(interval=1))

    plt.subplots_adjust(hspace=0.01)
    plt.savefig("./final_paper_plots/example_flare_ana.png", dpi=300, facecolor="w", bbox_inches="tight")
    plt.close()




# plot the 2012-11-20 12:36 flare
plot_flare(15)

# plot the 2015-10-02 15:26 flare
# plot_flare(262)
