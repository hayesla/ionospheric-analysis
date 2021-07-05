import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import scipy.stats
from sunpy.time import parse_time
from matplotlib import dates
from scipy.signal import savgol_filter
from sunpy import timeseries as ts 
from scipy.io import readsav 
from astropy.io import fits
import datetime

sid_file = '/Users/laurahayes/ionospheric_work/vlf_data_all_birr/sid_alll/20130522_000000_NAA_S-0055.csv'
goes_file = '/Users/laurahayes/QPP/stats_study/TEBBS/goes_rawdata/go1520130522.fits'

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

flare_start = "2013-05-22 12:00:00"
flare_end = "2013-05-22 15:00:00"


sid_data = sid_to_series(sid_file).truncate(flare_start, flare_end)
goes_data = ts.TimeSeries(goes_file).truncate(flare_start, flare_end)

gl = goes_data.to_dataframe()["xrsb"]
gs = goes_data.to_dataframe()["xrsa"]

def make_rhessi_lc(file):
    #file = 'hsi_spectrum_20130515_012024.fits'
    a = fits.open(file)
    start_time = a[0].header['DATE_OBS']
    t_start = datetime.datetime.strptime(start_time[0:10] + ' '+start_time[11:], '%Y-%m-%d %H:%M:%S.%f')
    start_time_day = datetime.datetime.strptime(str(t_start)[0:10]+' 00:00:00', '%Y-%m-%d %H:%M:%S')
    #print a[1].data.columns

    time = a[1].data['TIME']
    time_array = []
    for i in range(len(time)):
        time_array.append(start_time_day + datetime.timedelta(seconds = time[i]))


    e_min = a[2].data['E_MIN']
    e_max = a[2].data['E_MAX']
    channel = zip(list(e_min), list(e_max))
    dif_channels = a[1].data['RATE'].T

    r36 = dif_channels[0]
    r612 = dif_channels[1]
    r1225 = dif_channels[2]
    r2550 = dif_channels[3]
    r50100 = dif_channels[4]
    r100300 = dif_channels[5]


    df = pd.DataFrame({'3-6 keV': r36, '6-12 keV': r612, '12-25 keV': r1225, '25-50 keV': r2550,'50-100 keV': r50100, '100-300 keV': r100300},index=time_array)
    

    return df

rhe = make_rhessi_lc('hsi_spectrum_20130522_121644.fits').truncate(flare_start, flare_end)



atten = readsav('rhessi_atten_states_20130522.sav')['atten']
att_ts = [parse_time(x, format="utime").datetime for x in atten['start_times'][0]]
att_te = [parse_time(x, format="utime").datetime for x in atten['end_times'][0]]
state = atten['state'][0]

atten_info = np.array(list(zip(att_ts, att_te, state)))
aa_0 = atten_info[atten_info[:,2] == 0]
aa_1 = atten_info[atten_info[:,2] == 1]
aa_3 = atten_info[atten_info[:,2] == 3]
aa_4 = atten_info[atten_info[:,2] == 4]





rhe_data = readsav("rhessi_lc_corrected_20130522.sav")
rhe_times = parse_time(rhe_data["hsi_times"], format="utime").datetime
rhe_2550 = np.squeeze(rhe_data["hsi_counts_2550"])
rhe_50100 = np.squeeze(rhe_data["hsi_counts_50100"])
rhe_1225 = np.squeeze(rhe_data["hsi_counts_1225"])
rhe_612 = np.squeeze(rhe_data["hsi_counts_0612"])


fig, ax = plt.subplots(3, sharex=True)

ax[0].plot(gl, color="r")
ax[0].plot(gs, color="b")
ax[0].set_xlim(gl.index[0], gl.index[-1])
ax[0].set_yscale("log")

ax[1].plot(rhe_times, rhe_612, label="6-12keV")
ax[1].plot(rhe_times, rhe_1225, label="12-25keV")
ax[1].plot(rhe_times, rhe_2550, label="25-50keV")
ax[1].plot(rhe_times, rhe_50100, label="50-100keV")
ax[1].set_yscale("log")


ax[2].plot(calc_amp(sid_data), color="k")


def plot_atten_states(aa, y_val, color = 'r'):
	for i in range(len(aa)):
		plt.plot([aa[i][0], aa[i][1]], [y_val, y_val], color = color)

mask = np.logical_and(rhe.index < aa_0[0][1], rhe.index> aa_0[0][0])
for i in range(1, len(aa_0)):
	mask = mask + np.logical_and(rhe.index < aa_0[i][1], rhe.index > aa_0[i][0])
	
for i in range(0, 4):
	mask = mask + np.logical_and(rhe.index < aa_1[i][1], rhe.index > aa_1[i][0])



