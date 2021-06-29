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
from astropy.visualization import hist
import seaborn as sns

sns.set_context("paper")


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



def get_pandas_df_of_results():
    errors = []
    results = []
    save=True
    for i in range(0, len(vlf_flares)):
        new_ts = pd.to_datetime(vlf_flares["event_starttime"].iloc[i])-datetime.timedelta(minutes=5)
        new_te = pd.to_datetime(vlf_flares["event_endtime"].iloc[i])+datetime.timedelta(minutes=5)

        # SID data
        sid_file = glob.glob(vlf_flares.iloc[i]["event_starttime"].strftime(sid_file_dir))[0]
        sid_data = sid_to_series(sid_file).truncate(new_ts, new_te)
        sid_data_db = sid_to_series(sid_file, amp=True).truncate(new_ts, new_te)

        # smoothing window defined in terms of cadence
        window_sec =  (sid_data.index[1] - sid_data.index[0]).total_seconds()
        window = int((3*60)/window_sec)
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


        peak_vlf = np.max(sid_resample_flare)
        peak_vlf2 = np.abs(np.max(sid_resample_flare) - sid_resample_flare[0])

        peak_vlf_db = np.max(sid_resample_flare_db)
        peak_vlf2_db = np.abs(np.max(sid_resample_flare_db) - sid_resample_flare_db[0])

        dt_value_gl = (sid_resample_flare.index[np.argmax(sid_resample_flare)] - gl_flare.index[np.argmax(gl_flare)]).total_seconds()
        dt_value_gs = (sid_resample_flare.index[np.argmax(sid_resample_flare)] - gs_flare.index[np.argmax(gs_flare)]).total_seconds()


        ## plots
        fig, ax = plt.subplots(2, sharex=True)
        ax[0].plot(gl, color="r", label="1-8$\mathrm{\AA}$")
        ax[0].plot(gs, color="b", label="0.5-4$\mathrm{\AA}$")
        ax[0].set_ylabel("Flux (Wm$^{-2}$)")
        ax[0].legend(loc="upper left")
        ax[0].set_yscale("log")

        ax[1].plot(sid_data_db - sid_data_db[0], label="raw data", color="grey")
        ax[1].plot(sid_resample_flare_db - sid_resample_flare_db[0], label="2min resample", color="k")
        ax[1].legend(loc="upper left")      

        for a in ax:
            a.axvline(gl_flare.index[np.argmax(gl_flare)], color="r")
            a.axvline(gs_flare.index[np.argmax(gs_flare)], color="b")
            a.axvline(pd.to_datetime(vlf_flares["event_starttime"].iloc[i]), ls="dashed", color="grey")
            a.axvline(pd.to_datetime(vlf_flares["event_endtime"].iloc[i]), ls="dashed", color="grey")
            a.axvline(sid_resample_flare.index[np.argmax(sid_resample_flare)], color="k")
        
        tstart_str = pd.to_datetime(vlf_flares["event_starttime"].iloc[i]).strftime("%Y-%m-%dT%H:%M")
        ax[1].set_xlabel("Time {:s}".format(pd.to_datetime(vlf_flares["event_starttime"].iloc[i]).strftime("%Y-%m-%d %H:%M")))
        ax[1].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
        plt.tight_layout()

        plt.subplots_adjust(hspace=0.01)
        plt.savefig("./test_plots/flare_amp_3min_{:d}_{:s}.png".format(i, tstart_str))
        plt.close()
