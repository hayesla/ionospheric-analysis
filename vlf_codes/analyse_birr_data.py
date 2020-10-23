import pandas as pd 
import matplotlib.pyplot as plt 
from pathlib import Path
from sunpy import timeseries as ts
from sunpy.time import parse_time
from matplotlib import dates
import glob
import numpy as np 
import sys
sys.path.append("..")
from goes_event_list import get_goes_event_list 


def get_flarelist(goes_class_filter, filename):
    """
    getting flare list with flares greater than specifed flare class.
    Saved as a CSV file filename
    """ 
    t_start = "2012-08-22 00:00"
    t_end = "2018-04-20 00:00"
    get_goes_event_list(t_start, t_end, filename=Path.cwd().joinpath(filename), goes_class_filter=goes_class_filter)

# get_flarelist('C1', filename='goes_c_flares_birr_dates.csv')

goes_data_dir = "/Users/laurahayes/QPP/stats_study/TEBBS/goes_rawdata/"
vlf_data_dir = "/Users/laurahayes/ionospheric_work/vlf_data_all_birr/sid_alll/"
save_dir = "/Users/laurahayes/ionospheric_work/ionospheric-analysis/vlf_codes/vlf_plots_birr/"


# goes_flares = pd.read_csv("goes_sc_flares_xmc5.csv")
goes_flares = pd.read_csv("goes_c_flares_birr_dates.csv")
goes_flares = goes_flares.drop_duplicates(subset="start_time") 

goes_flares["peak_times_hours"] = [x.hour for x in pd.to_datetime(goes_flares["peak_time"])]
daytime_flares = goes_flares[(goes_flares["peak_times_hours"]>8) & (goes_flares["peak_times_hours"]<20)]

days_to_plot = daytime_flares["event_date"].unique()

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


def make_vlf_flare_list():
    vlf_days = []
    for i in range(len(days_to_plot)):
        tt = parse_time(days_to_plot[i]).strftime("%Y%m%d")
        files_vlf = glob.glob(vlf_data_dir + tt + '*.csv')
        if len(files_vlf) != 0:
            vlf_days.append(days_to_plot[i])


def plot(i):

    tt = parse_time(days_to_plot[i]).strftime("%Y%m%d")

    files_vlf = glob.glob(vlf_data_dir + tt + '*.csv')
    if len(files_vlf) == 0:
        print("No VLF data")
        return

    goes_file = goes_data_dir + "go15" + tt + ".fits"
    if not Path(goes_file).exists():
        print("No goes data")
        return

    data_vlf = read_files(files_vlf)
    goes_data = ts.TimeSeries(goes_file).to_dataframe()


    flares_ind = np.where(daytime_flares["event_date"].isin([days_to_plot[i]])==True)[0]
    flares = daytime_flares.iloc[flares_ind]


    fig, ax = plt.subplots(2, figsize=(8,6), sharex=True)

    ax[0].plot(goes_data['xrsb'], color='b', label='1-8$\mathrm{\AA}$')
    ax[0].plot(goes_data['xrsa'], color='r', label='0.5-4$\mathrm{\AA}$')
    ax[0].set_yscale('log')
    ax[0].set_xlim(days_to_plot[i] + " 00:00", days_to_plot[i] + " 23:59")

    ax[1].plot(pd.to_datetime(data_vlf['time']), data_vlf['volts'], color='grey')
    for f in flares["peak_time"]:
        ax[0].axvline(parse_time(f).datetime, color="k", ls="dashed")
        ax[1].axvline(parse_time(f).datetime, color="k", ls="dashed")

    ax[1].xaxis.set_major_locator(dates.HourLocator(interval=3))
    ax[1].xaxis.set_minor_locator(dates.HourLocator(interval=1))
    ax[1].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
    ax[0].set_ylim(1e-9, 1e-3)
    #ax[1].set_ylim(-5, 5)

    for a in ax:
        a.tick_params(which='both', direction='in')

    ax[0].set_ylabel('Flux Wm$^{-2}$')
    ax[1].set_ylabel('Volts')
    ax[1].set_xlabel('Time ' + days_to_plot[i] + ' UT')
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.05)
    plt.savefig(save_dir + 'birr_vlf_' +  days_to_plot[i] + '.png', dpi=100)
    plt.close()


    