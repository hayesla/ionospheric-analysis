from sunpy.time import parse_time
from sunpy import timeseries as ts
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib import dates
from pathlib import Path
import glob
import datetime
import urllib
import os
import numpy as np 
import sys
sys.path.append("..")
from goes_event_list import get_goes_event_list 

goes_data_dir = "/Users/laurahayes/QPP/stats_study/TEBBS/goes_rawdata/"

def get_flarelist(goes_class_filter, filename):
    """
    getting flare list with flares greater than specifed flare class.
    Saved as a CSV file filename
    """ 
    t_start = "2011-02-01 00:00"
    t_end = "2019-01-01 00:00"
    get_goes_event_list(t_start, t_end, filename=Path.cwd().joinpath(filename), goes_class_filter=goes_class_filter)




goes_flares = pd.read_csv("goes_sc_flares_xmc5.csv")
goes_flares = goes_flares.drop_duplicates(subset="start_time") 
goes_flares["peak_times_hours"] = [x.hour for x in pd.to_datetime(goes_flares["peak_time"])]
daytime_flares = goes_flares[(goes_flares["peak_times_hours"]>8) & (goes_flares["peak_times_hours"]<20)]

events_to_download = daytime_flares["event_date"].unique()


def check_url(url):
    try: 
        urllib.request.urlopen(url)
        return True
    except:
        print("nah")
        return False


def check_urls_bas(transmitter='NAA'):
    url_list = []
    failed = []
    for i in range(len(events_to_download)):
        date = parse_time(events_to_download[i])
        url_base = "http://psddb.nerc-bas.ac.uk/data/psddata/atmos/space/vlf/ultra/eskdalemuir//"
        url_extent = "{:s}/data/{:s}{:s}.txt".format(date.strftime("%Y"), transmitter, date.strftime("%Y%m%d"))
        url1 = url_base + url_extent
        url2 = url_base + url_extent + ".gz"
        aa = check_url(url1)
        bb = check_url(url2)
        if aa:
            print("url_exists!")
            url_list.append(url1)
            print(aa)
        elif bb:
            print("url exists!")
            url_list.append(url2)
            print(bb)
        else:
            print("nay")
            failed.append(events_to_download[i])

    return url_list, failed


# list_of_files = check_urls_bas(transmitter='DHO')
# worked, failed = list_of_files


from parfive import Downloader
def download_files_parfive(list_of_files, pathy="/Users/laurahayes/ionospheric_work/ionospheric-analysis/vlf_codes/vlf_bas_files/"):
    dl = Downloader()
    for f in list_of_files:
        filename = f.split('/')[-1]
        dl.enqueue_file(f, path=pathy)


    files = dl.download()

def download_files(list_of_files, pathy="/Users/laurahayes/ionospheric_work/ionospheric-analysis/vlf_codes/vlf_bas_files/"):
    for f in list_of_files:
        print(f)
        if not os.path.exists(pathy + f.split('/')[-1]):

            try:
                urllib.request.urlretrieve(f, pathy + f.split('/')[-1])
            except:
                print("problemo!")



def read_vlf_data(file, tt):

    data = pd.read_csv(file, comment="%", names=["tsec", "amp", "phase"], delim_whitespace=True)
    date = parse_time(tt[0:4]+'-'+tt[4:6] + '-'+tt[6:8]).datetime
    times = [date + datetime.timedelta(seconds=t) for t in data['tsec']]
    amp = pd.Series(data['amp'].values, index=times)
    phase = pd.Series(data['phase'].values, index=times)
    



    return amp, phase


goes_data_dir = "/Users/laurahayes/QPP/stats_study/TEBBS/goes_rawdata/"
save_dir = "/Users/laurahayes/ionospheric_work/ionospheric-analysis/vlf_codes/vlf_plots/"

def get_flare_amps(i):
    """
    Function to return the GOES amp and the VLF amp

    """

    all_data = []
    for i in range(62, len(events_to_download)):
        tt = parse_time(events_to_download[i]).strftime("%Y%m%d")
        files_vlf = glob.glob("./vlf_bas_files/{:s}{:s}*".format("NAA", tt))
        if len(files_vlf)>0:
            print(i, 'out of ', len(events_to_download))
            # print("No VLF data")



        # goes_file = goes_data_dir + "go15" + tt + ".fits"
        # if not Path(goes_file).exists():
        #     print("No goes data")

        # goes_data = ts.TimeSeries(goes_file)
        # gl = goes_data.data["xrsb"]
        # gs = goes_data.data["xrsa"]
            flares_ind = np.where(daytime_flares["event_date"].isin([events_to_download[i]])==True)[0]
            flares = daytime_flares.iloc[flares_ind]

            vlf_amp, vlf_phase = read_vlf_data(files_vlf[0], tt)

            day_date = parse_time(vlf_amp.index[0]).strftime('%Y-%m-%d')
            vlf_amp.sort_index(inplace=True)
            average_vlf_amp = vlf_amp.truncate(day_date + ' 09:00', day_date + ' 17:00')
            average_vlf_amp_mean = average_vlf_amp.mean()

            df_test = {}
            for f in range(len(flares)):
                amp_flare = vlf_amp.truncate(flares.iloc[f]['start_time'], flares.iloc[f]['peak_time'])
                amp_max = amp_flare.max()
                if len(amp_flare)>0:
                    amp_diff = amp_flare[-1] - amp_flare[0]
                else:
                    amp_diff = np.nan
                df_test['amp_flare'] = amp_max
                df_test['amp_diff'] = amp_diff
                df_test['amp_average'] = average_vlf_amp_mean
                df_test['flare'] = flares.iloc[f]['goes_class'] 
                df_test['start_time'] = flares.iloc[f]['start_time']
                df_test['peak_time'] = flares.iloc[f]['peak_time']



        all_data.append(df_test)
        all_data_vlf = pd.DataFrame(all_data)

def get_goes_flux(x):

        if x[0] == 'C':
            gc = float(x[1:])*(1e-6)
        elif x[0] == 'M':
            gc = float(x[1:])*(1e-5)    
        elif x[0] == 'X':
            gc = float(x[1:])*(1e-4)    
        else:
            print ('somethings wrong')
        return gc


def plot_flares(i, transmitter='NRK'):

    tt = parse_time(events_to_download[i]).strftime("%Y%m%d")
    files_vlf = glob.glob("./vlf_bas_files/{:s}{:s}*".format(transmitter, tt))
    if len(files_vlf)==0:
        print("No VLF data")
        return

    goes_file = goes_data_dir + "go15" + tt + ".fits"
    if not Path(goes_file).exists():
        print("No goes data")
        return

    goes_data = ts.TimeSeries(goes_file)
    gl = goes_data.data["xrsb"]
    gs = goes_data.data["xrsa"]
    flares_ind = np.where(daytime_flares["event_date"].isin([events_to_download[i]])==True)[0]
    flares = daytime_flares.iloc[flares_ind]

    vlf_amp, vlf_phase = read_vlf_data(files_vlf[0], tt)


    fig, ax = plt.subplots(3, sharex=True, figsize=(8, 10))

    ax1 = ax[0]
    ax2 = ax[1]
    ax3 = ax[2]

    ax1.plot(gl, color="r", label="1-8 $\mathrm{\AA}$")
    ax1.plot(gs, color="b", label="0.5-4 $\mathrm{\AA}$")   
    ax1.set_ylim(1e-9, 1e-3)
    ax1.set_yscale("log")
    ax1.tick_params(which="both", direction="in", right=True, top=True)
    ax1.set_ylabel("Flux (Wm$^{-2}$)")
    ax1.legend(loc="upper right")

    ax2.plot(vlf_amp, label='NAA', color='grey')
    ax2.set_ylabel('VLF Amplitude (dB)')

    ax3.plot(vlf_phase, label='NAA', color='k')
    ax3.set_ylabel('Phase (degrees)')


    ax3.set_xlabel("Time {:s} UT".format(events_to_download[i]))
    ax3.set_xlim(events_to_download[i] + " 00:00", events_to_download[i] + " 23:59")
    ax3.xaxis.set_major_locator(dates.HourLocator(interval=3))
    ax3.xaxis.set_minor_locator(dates.HourLocator(interval=1))
    ax3.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
    ax3.tick_params(which="both", direction="in", right=True, top=True)
    for f in flares["peak_time"]:
        ax1.axvline(parse_time(f).datetime, color="k", ls="dashed")
        ax2.axvline(parse_time(f).datetime, color="k", ls="dashed")
        ax3.axvline(parse_time(f).datetime, color="k", ls="dashed")
    ax1.grid()
    ax2.grid()
    ax3.grid()
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.05)
    plt.savefig(save_dir + transmitter + parse_time(events_to_download[i]).strftime("%Y%m%d.png"), dpi=200)
    plt.close()


def plot_all():
    for i in range(len(events_to_download)):
        print("doing flare ", i)
        try:
            plot_flares(i, transmitter='DHO')
        except:
            print(i, "failed")