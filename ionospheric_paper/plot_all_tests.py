import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 
import glob 
from sunpy.time import parse_time 
from sunpy import timeseries as ts 
import datetime
from matplotlib import dates
from scipy.signal import savgol_filter
import scipy.stats
import plotly.express as px
import sunpy.map
from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter
from astropy.visualization import hist
import seaborn as sns
sns.set_context("paper")


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


def get_test(i):
    new_ts = pd.to_datetime(vlf_flares["event_starttime"].iloc[i])-datetime.timedelta(minutes=5)
    new_te = pd.to_datetime(vlf_flares["event_endtime"].iloc[i])+datetime.timedelta(minutes=5)

    # SID data
    sid_file = glob.glob(vlf_flares.iloc[i]["event_starttime"].strftime(sid_file_dir))[0]

    sid_data = sid_to_series(sid_file).truncate(new_ts, new_te)
    full_sid = sid_to_series(sid_file, amp=True)
    full_sid_volts = sid_to_series(sid_file, amp=False)
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
    goes = ts.TimeSeries(goes_file)
    gl = goes.truncate(new_ts, new_te).to_dataframe()["xrsb"]
    gs = goes.truncate(new_ts, new_te).to_dataframe()["xrsa"]
    gl_flare = gl.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])
    gs_flare = gs.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])


    peak_vlf = np.max(sid_resample_flare)
    peak_vlf2 = np.abs(np.max(sid_resample_flare) - sid_resample_flare[0])

    peak_vlf_db = np.max(sid_resample_flare_db)
    peak_vlf2_db = np.abs(np.max(sid_resample_flare_db) - sid_resample_flare_db[0])

    dt_value_gl = (sid_resample_flare.index[np.argmax(sid_resample_flare)] - gl_flare.index[np.argmax(gl_flare)]).total_seconds()
    dt_value_gs = (sid_resample_flare.index[np.argmax(sid_resample_flare)] - gs_flare.index[np.argmax(gs_flare)]).total_seconds()

    event = {}
    event["start_time_goes"] = vlf_flares.iloc[i]["event_starttime"]
    event["peak_flare_gl"] = np.max(gl_flare)
    event["peak_flare_gs"] = np.max(gs_flare)
    event["peak_flare_gl_bg"] = np.max(gl_flare - gl_flare[0])
    event["peak_flare_gs_bg"] = np.max(gs_flare - gs_flare[0])
    event["max_vlf"] = peak_vlf
    event["abs_vlf"] = peak_vlf2
    event["max_vlf_db"] = peak_vlf_db
    event["abs_vlf_db"] = peak_vlf2_db
    event["dt_value_gl"] = dt_value_gl
    event["dt_value_gs"] = dt_value_gs


    # plot in amps
    fig, ax = plt.subplots(4, figsize=(8, 8))
    ax[0].plot(goes.to_dataframe()["xrsb"], color="r")
    ax[0].plot(goes.to_dataframe()["xrsa"], color="b")
    ax[0].set_yscale("log")
    ax[1].plot(full_sid, color="k")

    for a in [ax[0], ax[1]]:
        a.axvline(new_ts, color="grey")
        a.axvline(new_te, color="grey")
        a.axvline(vlf_flares["event_starttime"].iloc[i], ls='dashed', color="grey")
        a.axvline(vlf_flares["event_endtime"].iloc[i], ls='dashed', color="grey")
        a.set_xlim(pd.to_datetime(full_sid.index[0].strftime("%Y-%m-%d 06:00")), pd.to_datetime(full_sid.index[0].strftime("%Y-%m-%d 21:00")))
        a.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))

    ax[2].plot(gl_flare - gl_flare[0], color="r")
    ax[2].plot(gs_flare - gs_flare[0], color="b")
    ax[2].set_yscale("log")
    
    ax[3].plot(sid_data_db.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])-sid_resample_flare_db[0], color="grey")
    ax[3].plot(sid_resample_flare_db - sid_resample_flare_db[0], color="k")
    

    for a in (ax[2], ax[3]):
        a.set_xlim(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])
        a.axvline(sid_resample_flare_db.index[np.argmax(sid_resample_flare_db)], color="k")
        a.axvline(gl_flare.index[np.argmax(gl_flare)], color="r")
        a.axvline(gs_flare.index[np.argmax(gs_flare)], color="b")
        a.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))

    ax[3].set_xlabel("Time {:s}".format(pd.to_datetime(vlf_flares["event_starttime"].iloc[i]).strftime("%Y-%m-%d %H:%M")))
    plt.tight_layout()
    tstart_str = pd.to_datetime(vlf_flares["event_starttime"].iloc[i]).strftime("%Y-%m-%dT%H:%M")
    plt.savefig("./final_test_plots/flare_db_{:s}.png".format(tstart_str))
    plt.close()


    # plot in volts
    fig, ax = plt.subplots(4, figsize=(8, 8))
    ax[0].plot(goes.to_dataframe()["xrsb"], color="r")
    ax[0].plot(goes.to_dataframe()["xrsa"], color="b")
    ax[0].set_yscale("log")
    
    ax[1].plot(full_sid_volts, color="k")

    for a in [ax[0], ax[1]]:
        a.axvline(new_ts, color="grey")
        a.axvline(new_te, color="grey")
        a.axvline(vlf_flares["event_starttime"].iloc[i], ls='dashed', color="grey")
        a.axvline(vlf_flares["event_endtime"].iloc[i], ls='dashed', color="grey")
        a.set_xlim(pd.to_datetime(full_sid_volts.index[0].strftime("%Y-%m-%d 06:00")), pd.to_datetime(full_sid_volts.index[0].strftime("%Y-%m-%d 21:00")))
        a.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))

    ax[2].plot(gl_flare - gl_flare[0], color="r")
    ax[2].plot(gs_flare - gs_flare[0], color="b")
    ax[2].set_yscale("log")
    
    ax[3].plot(sid_data.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])-sid_resample_flare[0], color="grey")
    ax[3].plot(sid_resample_flare - sid_resample_flare[0], color="k")
    

    for a in (ax[2], ax[3]):
        a.set_xlim(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])
        a.axvline(sid_resample_flare.index[np.argmax(sid_resample_flare)], color="k")
        a.axvline(gl_flare.index[np.argmax(gl_flare)], color="r")
        a.axvline(gs_flare.index[np.argmax(gs_flare)], color="b")
        a.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))

    ax[3].set_xlabel("Time {:s}".format(pd.to_datetime(vlf_flares["event_starttime"].iloc[i]).strftime("%Y-%m-%d %H:%M")))
    plt.tight_layout()
    tstart_str = pd.to_datetime(vlf_flares["event_starttime"].iloc[i]).strftime("%Y-%m-%dT%H:%M")
    plt.savefig("./final_test_plots/flare_volts_{:s}.png".format(tstart_str))
    plt.close()




for i in range(len(vlf_flares)):
    print(i)
    get_test(i)

