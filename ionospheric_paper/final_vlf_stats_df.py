import pandas as pd 
import numpy as np 
from sunpy.time import parse_time
import datetime
import glob
from scipy.signal import savgol_filter
from sunpy import timeseries as ts
import matplotlib.pyplot as plt
from matplotlib import dates
# Util functions used here

def calc_amp(x):
    """
    Convert volts measurments from SID to dB.

    """
    return 20*np.log10(x + 5) - 61 + 107

def sid_to_series(file, amp=False):
    """
    Function to read the VLF txt files and return a pandas Series.

    Parameters
    ----------
    file : ~`str`, or path
        file to be read
    amp : Boolean, optional, default False
        if True then will convert volts to amp 

    Returns 
    -------
    pd.Series

    """
    sid = pd.read_csv(file, comment="#", names=["times", "data"])
    tt = parse_time(sid["times"]).datetime
    if amp:
        ser = pd.Series(calc_amp(sid["data"].values), index=tt)
    else:
        ser = pd.Series(sid["data"].values, index=tt)       
    ser.sort_index(inplace=True)
    return ser


def final_runthrough():
    """
    Get the list of final flares to analyse - those with pngs left in the 
    final_flares dir.

    Saves list of flares to final_paper_vlf_flares2.csv

    """
    flare_list = pd.read_csv("flare_list_for_vlf.csv")
    flare_list["event_starttime"] = pd.to_datetime(flare_list["event_starttime"])

    gg = glob.glob("/Users/laurahayes/ionospheric_work/ionospheric-analysis/stats_study/final_flares/*.png")
    gg.sort()
    flare_starttimes = [x.split("/")[-1][6:22] for x in gg]
    matchtimes = [pd.to_datetime(x).strftime("%Y-%m-%d %H:%M:%S") for x in flare_starttimes]

    final_df = flare_list[flare_list["event_starttime"].isin([matchtimes[0]])]
    for i in range(1, len(matchtimes)):
        final_df = final_df.append(flare_list[flare_list["event_starttime"].isin([matchtimes[i]])])

    final_df.reset_index(inplace=True, drop=True)
    final_df.to_csv("final_paper_vlf_flares2.csv", index_label=False)



# File patterns
sid_file_dir = "/Users/laurahayes/ionospheric_work/vlf_data_all_birr/sid_alll/*%Y%m%d*NAA*"
goes_file_dir = "/Users/laurahayes/QPP/stats_study/TEBBS/goes_rawdata/*%Y%m%d*.fits"

# read final flares to analyse 
vlf_flares = pd.read_csv("final_paper_vlf_flares2.csv")
vlf_flares["event_starttime"] = pd.to_datetime(vlf_flares["event_starttime"])

##--------------Fix some backgrounds---------------------------##

# Need to fix time to choose background for several flares
# for remainder use the GOES start time as the background time.
vlf_flares["background_times"] = vlf_flares["event_starttime"].astype(str)

# updates to some bgtimes
update_bgtime = {'2013-04-25 17:24:00': '2013-04-25 17:13:00',
                 '2013-05-22 13:08:00': '2013-05-22 12:30:00',
                 '2013-10-25 14:51:00': '2013-10-25 14:30:00',
                 '2015-03-05 17:06:00': '2015-03-05 18:00:00',
                 '2015-03-12 12:09:00': '2015-03-12 11:38:00',
                 '2015-08-07 14:01:00': '2015-08-07 14:15:00',
                 '2015-10-02 17:55:00': '2015-10-02 17:56:00',
                 '2016-07-24 12:18:00': '2016-07-24 11:33:00',
                 '2016-07-24 12:56:00': '2016-07-24 11:33:00',
                 '2016-07-24 13:15:00': '2016-07-24 11:33:00',
                 '2016-07-24 13:34:00': '2016-07-24 11:33:00',
                 '2016-08-07 15:03:00': '2016-08-07 14:40:00'}

vlf_flares["background_times"] = vlf_flares["background_times"].apply(lambda x: x if x not in update_bgtime.keys() else update_bgtime[x])

##--------------remove bad events---------------------------##

events_to_remove = ["2016-08-07 15:03:00",
                    "2015-10-17 14:22:00",
                    "2015-04-22 08:30:00",
                    "2013-10-17 11:47:00", 
                    "2016-02-15 10:41:00",
                    "2017-04-02 07:50:00"]
vlf_flares.drop(np.where(vlf_flares["event_starttime"].isin(events_to_remove))[0], inplace=True)
vlf_flares.reset_index(drop=True, inplace=True)

##--------------fix_dropouts---------------------------##

# This is a shameful way to fix dropouts - manually doing it
# arg - just not arsed thinking of a smarter way for 5 events.
dropout_events = ["2013-05-03 16:39:00",
                  "2015-09-20 17:32:00",
                  "2016-07-19 10:00:00", 
                  "2015-09-30 10:49:00",
                  "2016-07-24 13:34:00"]

timerange_fix_start = ["2013-05-03 17:17:55", 
                       "2015-09-20 17:52:00",
                       "2016-07-19 12:16:50", 
                       "2015-09-30 11:09:55",
                       "2016-07-24 14:10:00"]

timerange_fix_end =   ["2013-05-03 17:18:55", 
                       "2015-09-20 17:53:00",
                       "2016-07-19 12:18:20", 
                       "2015-09-30 11:10:50",
                       "2016-07-24 14:10:30"]


fixes =pd.DataFrame({"event_starttime": dropout_events, 
                     "timerange_start": timerange_fix_start, 
                     "timerange_end": timerange_fix_end})



def get_pandas_df_of_results():
    errors = []
    results = []
    save=True
    for i in range(0, len(vlf_flares)):
        print("doing {:d}".format(i))
        try:
        # if save:
            new_ts = pd.to_datetime(vlf_flares["event_starttime"].iloc[i])-datetime.timedelta(minutes=5)
            new_te = pd.to_datetime(vlf_flares["event_endtime"].iloc[i])+datetime.timedelta(minutes=5)

            # SID data
            sid_file = glob.glob(vlf_flares.iloc[i]["event_starttime"].strftime(sid_file_dir))[0]
            sid_data_volts = sid_to_series(sid_file)
            sid_data_db = sid_to_series(sid_file, amp=True)

            # fix the dropout events
            if str(vlf_flares["event_starttime"].iloc[i]) in list(fixes["event_starttime"]):
                fix = fixes[fixes["event_starttime"].isin([str(vlf_flares["event_starttime"].iloc[i])])]
                sid_data_volts[(sid_data_volts.index>=fix["timerange_start"].values[0])&(sid_data_volts.index<=fix["timerange_end"].values[0])] = np.nan
                sid_data_db[(sid_data_db.index>=fix["timerange_start"].values[0])&(sid_data_db.index<=fix["timerange_end"].values[0])] = np.nan
                sid_data_volts.fillna(method="ffill", inplace=True)
                sid_data_db.fillna(method="ffill", inplace=True)

            # smoothing window defined in terms of cadence
            window_sec =  (sid_data_volts.index[1] - sid_data_volts.index[0]).total_seconds()
            window = int((2*60)/window_sec)
            if window%2 == 0:
                window = window+1

            # resample the SID data to get smoother timeseries
            # volts
            sid_resample_volts = pd.Series(savgol_filter(sid_data_volts, int(window), 3), index=sid_data_volts.index)
            sid_resample_volts_flare = sid_resample_volts.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])
            # dB
            sid_resample_db = pd.Series(savgol_filter(sid_data_db, int(window), 3), index=sid_data_db.index)
            sid_resample_flare_db = sid_resample_db.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])                
            
            # GOES data
            goes_file = glob.glob(pd.to_datetime(vlf_flares["event_starttime"].iloc[i]).strftime(goes_file_dir))[0]
            goes = ts.TimeSeries(goes_file)#.truncate(new_ts, new_te)
            gl = goes.to_dataframe()["xrsb"]
            gs = goes.to_dataframe()["xrsa"]
            gl_flare = gl.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])
            gs_flare = gs.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])


            # max amp of VLF in volts
            background_volts = sid_resample_volts[sid_resample_volts.index.get_loc(vlf_flares["background_times"].iloc[i], method='nearest')]
            peak_vlf_volts = np.max(sid_resample_volts_flare)
            peak_vlf2_volts = np.abs(np.max(sid_resample_volts_flare) - background_volts)
            
            # max amp of VLF in volts
            background_db = sid_resample_db[sid_resample_volts.index.get_loc(vlf_flares["background_times"].iloc[i], method='nearest')]
            peak_vlf_db = np.max(sid_resample_flare_db)
            peak_vlf2_db = np.abs(np.max(sid_resample_flare_db) - background_db)

            # time delays
            dt_value_gl = (sid_resample_flare_db.index[np.argmax(sid_resample_flare_db)] - gl_flare.index[np.argmax(gl_flare)]).total_seconds()
            dt_value_gs = (sid_resample_flare_db.index[np.argmax(sid_resample_flare_db)] - gs_flare.index[np.argmax(gs_flare)]).total_seconds()


            # write data to row
            event = {}
            event["start_time_goes"] = vlf_flares.iloc[i]["event_starttime"]
            event["peak_flare_gl"] = np.max(gl_flare)
            event["peak_flare_gs"] = np.max(gs_flare)
            event["peak_flare_gl-bg"] = np.max(gl_flare - gl_flare[0])
            event["peak_flare_gs-bg"] = np.max(gs_flare - gs_flare[0])
            event["max_vlf_volts"] = peak_vlf_volts
            event["abs_vlf_volts"] = peak_vlf2_volts
            event["background_sid_db"] = background_db
            event["background_sid_volts"] = background_volts
            event["background_gl"] = gl_flare[0]
            event["background_gs"] = gs_flare[0]
            event["max_vlf_db"] = peak_vlf_db
            event["abs_vlf_db"] = peak_vlf2_db
            event["dt_value_gl"] = dt_value_gl
            event["dt_value_gs"] = dt_value_gs


            results.append(event)


            # plot in amps
            fig, ax = plt.subplots(4, figsize=(8, 8))
            ax[0].plot(goes.to_dataframe()["xrsb"], color="r")
            ax[0].plot(goes.to_dataframe()["xrsa"], color="b")
            ax[0].set_yscale("log")
            ax[1].plot(sid_data_db, color="k")

            for a in [ax[0], ax[1]]:
                a.axvline(new_ts, color="grey")
                a.axvline(new_te, color="grey")
                a.axvline(vlf_flares["event_starttime"].iloc[i], ls='dashed', color="grey")
                a.axvline(vlf_flares["event_endtime"].iloc[i], ls='dashed', color="grey")
                a.set_xlim(pd.to_datetime(sid_data_db.index[0].strftime("%Y-%m-%d 06:00")), pd.to_datetime(sid_data_db.index[0].strftime("%Y-%m-%d 21:00")))
                a.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))

            ax[2].plot(gl_flare - gl_flare[0], color="r")
            ax[2].plot(gs_flare - gs_flare[0], color="b")
            ax[2].set_yscale("log")
            
            ax[3].plot(sid_data_db.truncate(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])-sid_resample_flare_db[0], color="grey")
            ax[3].plot(sid_resample_flare_db - background_db, color="k")
            

            for a in (ax[2], ax[3]):
                a.set_xlim(vlf_flares["event_starttime"].iloc[i], vlf_flares["event_endtime"].iloc[i])
                a.axvline(sid_resample_flare_db.index[np.argmax(sid_resample_flare_db)], color="k")
                a.axvline(gl_flare.index[np.argmax(gl_flare)], color="r")
                a.axvline(gs_flare.index[np.argmax(gs_flare)], color="b")
                a.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))

            ax[3].set_xlabel("Time {:s}".format(pd.to_datetime(vlf_flares["event_starttime"].iloc[i]).strftime("%Y-%m-%d %H:%M")))
            plt.tight_layout()
            tstart_str = pd.to_datetime(vlf_flares["event_starttime"].iloc[i]).strftime("%Y-%m-%dT%H:%M")
            plt.savefig("./final_checks2/flare_db_{:s}.png".format(tstart_str))
            plt.close()

        except:
            print("error", i)
    

        





    results = pd.DataFrame(results)





    merged_results = pd.merge(results, vlf_flares, left_on="start_time_goes", right_on="event_starttime")
    merged_results.to_csv("vlf_stats_results_new.csv", index_label=False)

