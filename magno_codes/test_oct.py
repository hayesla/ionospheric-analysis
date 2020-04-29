from sunpy.time import parse_time
from sunpy import timeseries as ts
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib import dates
from pathlib import Path
import glob
from download_BIRR_data import search_data_and_download, find_files
import datetime
import sys
sys.path.append('..')
from goes_event_list import get_goes_event_list 
from scipy.io import readsav
import numpy as np

goes_data_dir = '/Users/laurahayes/QPP/stats_study/TEBBS/goes_rawdata/'

t_s = '2013-10-28 09:00'
t_e = '2013-10-28 20:00'

goes_file = goes_data_dir + 'go1520131028.fits'
goes_data = ts.TimeSeries(goes_file)
gl = goes_data.data['xrsb'].truncate(t_s, t_e)
gs = goes_data.data['xrsa'].truncate(t_s, t_e)


filey = '/Users/laurahayes/ionospheric_work/ionospheric-analysis/magno_codes/magno_files/birr_mag_20131028_000001.txt'
magno = pd.read_csv(filey, delim_whitespace=True, skiprows=1,
					names=['Date','Time', 'Index', 'Bx', 'By', 'Bz', 'E1', 'E2', 'E3', 'E4', 'T(FG)', 'T(E)', 'volts'])

magno_time = [datetime.datetime.strptime(magno.iloc[i]['Date'] + ' ' +magno.iloc[i]['Time'], '%d/%m/%Y %H:%M:%S') for i in range(len(magno))]

bx = pd.Series(np.array(magno['Bx']), index=magno_time).truncate(t_s, t_e)
by = pd.Series(np.array(magno['By']), index=magno_time).truncate(t_s, t_e)
bz = pd.Series(np.array(magno['Bz']), index=magno_time).truncate(t_s, t_e)

h = np.sqrt(np.array(bx)**2 + np.array(by)**2)
H = pd.Series(h, index = magno_time)

euve = readsav('idlsave_goes_28oct.sav')
euve_times = [parse_time(euve['UTBASE']+euve['tarray'][i], format='utime').datetime for i in range(len(euve['tarray']))]
eve_gt = pd.Series(euve['yclean'][0], index=euve_times)

evve = eve_gt.replace(-99999.0, np.nan).truncate(t_s, t_e)

def normalize(x):
	return (x - x.min())/(x.max() - x.min())

fig, ax = plt.subplots(3, sharex=True, figsize=(6, 10))

ax1 = ax[0]
ax2 = ax[1]
ax3 = ax[2]

ax1.plot(gl, color='r', label='1-8 $\mathrm{\AA}$')
ax1.plot(gs, color='b', label='0.5-4 $\mathrm{\AA}$')	
ax1.set_ylim(1e-9, 1e-3)
ax1.set_yscale('log')
ax1.tick_params(which='both', direction='in', right=True, top=True)
ax1.set_ylabel('Flux (Wm$^{-2}$)')
ax1.legend(loc='upper right')
ax1.set_ylabel('Flux')

ax2.plot(normalize(bx), label='Bx', color='k')
ax2.plot(normalize(by), label='By', color='grey')
ax2.plot(normalize(bz), label='Bz', color='green')
ax2.set_ylabel('Normalized')
ax2.set_ylim(0.6, -0.1)
ax2.legend()
ax3.plot(evve, label='EUVE')
ax3.legend()
ax3.set_xlim(t_s, t_e)
ax3.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
ax3.set_ylabel('Flux')
ax3.set_xlabel('Time 2013-10-28 UT')
plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.savefig('oct_test2.png', dpi=200)
plt.close()




fig, ax = plt.subplots(2, sharex=True, figsize=(10, 8))

ax1 = ax[0]
ax2 = ax[1]

ax1.plot(gl, color='r', label='1-8 $\mathrm{\AA}$')
ax1.plot(gs, color='b', label='0.5-4 $\mathrm{\AA}$')	
ax1.set_ylim(1e-9, 1e-3)
ax1.set_yscale('log')
ax1.tick_params(which='both', direction='in', right=True, top=True)
ax1.set_ylabel('Flux (Wm$^{-2}$)')
ax1.legend(loc='upper right')


ax2.plot(magno_time, magno['Bx'], label='Bx', color='k')
ax2.plot(np.nan, color='grey', label='By')
ax2.plot(np.nan, color='green', label='Bx')
ax2.legend(loc='lower right')
ax3 = ax2.twinx()
ax3.plot(magno_time, magno['By'], label='By', color='grey')
ax4 = ax2.twinx()
ax4.plot(magno_time, magno['Bz'], label='Bz', color='green')

ax2.set_xlabel('Time 2013-10-28 UT')
ax2.set_xlim(gl.index[0], gl.index[-1])
ax2.xaxis.set_major_locator(dates.HourLocator(interval=3))
ax2.xaxis.set_minor_locator(dates.HourLocator(interval=1))
ax2.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
ax2.tick_params(which='both', direction='in', right=True, top=True)

ax1.grid()
ax2.grid()
plt.tight_layout()
plt.subplots_adjust(hspace=0.05)
plt.savefig('oct_test.png', dpi=200)
plt.close()

