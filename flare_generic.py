import matplotlib.pyplot as plt 
from matplotlib import dates
import numpy as np 
from sunpy import timeseries as ts 
from sunpy.time import parse_time
from read_files import euve_to_series, mag_to_series, sid_to_series

# tstart = '2017-09-10 15:00'
# tend = '2017-09-10 22:00'

tstart = '2015-11-04 08:00'
tend = '2015-11-04 20:00'

euve_data = euve_to_series("./magno_codes/euve_data/g15_euve_{:s}.txt".format(parse_time(tstart).strftime('%Y%m%d')))
goes_data = ts.TimeSeries("/Users/laurahayes/QPP/stats_study/TEBBS/goes_rawdata/go15{:s}.fits".format(parse_time(tstart).strftime('%Y%m%d')))
magno_data = mag_to_series("./magno_codes/magno_files/birr_mag_{:s}_000001.txt".format(parse_time(tstart).strftime('%Y%m%d')))
sid_data = sid_to_series("./vlf_codes/vlf_files/BIR_sid_{:s}_000000.txt".format(parse_time(tstart).strftime('%Y%m%d')))

#euve_flare = euve_data.truncate(tstart, tend)
bx, by, bz = magno_data[0].truncate(tstart, tend), magno_data[1].truncate(tstart, tend), magno_data[2].truncate(tstart, tend)
gl, gs = goes_data.to_dataframe().truncate(tstart, tend)['xrsb'], goes_data.to_dataframe().truncate(tstart, tend)['xrsa']
euve_data = euve_data.truncate(tstart, tend)

h = np.sqrt(np.array(bx)**2 + np.array(by)**2)
H = pd.Series(h, index=bx.index)

dd = np.arctan(by/bx)*(180/np.pi)
D = pd.Series(dd, index=bx.index)

def norm(x):
	return (x - np.min(x))/(np.max(x) - np.min(x))


def overall_plot():
	fig, ax = plt.subplots(3, sharex=True, figsize=(8, 10))

	ax[0].plot(gl, color='r', label='GOES 1-8$\mathrm{\AA}$')
	ax[0].plot(np.nan, color='k', label=r'Ly$\alpha$')
	ax[0].set_yscale('log')
	ax[0].spines['bottom'].set_color('red')
	ax[0].spines['top'].set_color('red')
	ax[0].xaxis.label.set_color('red')



	ax02 = ax[0].twinx()
	ax02.plot(sid_data, color='k', label=r'Ly$\alpha$', drawstyle='steps-mid')
	ax02.set_ylabel('Flux (Wm$^{-2}$)')
	ax[0].set_title('a. Flare emission', loc='left')
	ax[0].set_ylabel("Flux (Wm$^{-2}$)", color='r')





	ax[1].plot(norm(bx), label='Bx')
	ax[1].plot(norm(by), label='By')
	ax[1].plot(norm(bz), label='Bz')
	ax[1].set_title('b. Magnetometer data', loc='left')
	ax[1].set_ylabel("Normalized components")

	ax[2].plot(euve_data, label='VLF', color='grey')
	ax[2].set_ylabel('VLF amplitude (Volts)')


	ax[2].set_xlim(tstart, tend)
	ax[2].set_xlabel('Time (UT) {:s}'.format(tstart))
	ax[2].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
	ax[2].xaxis.set_minor_locator(dates.MinuteLocator(interval=5))
	ax[2].xaxis.set_major_locator(dates.MinuteLocator(interval=15))

	ax[2].set_title('c. VLF data', loc='left')

	for a in ax:
		a.tick_params(which='both', direction='in', top=True)
		a.legend(loc='upper left')
		a.xaxis.grid()


	plt.tight_layout()


	plt.savefig('flare_{:s}.png'.format(tstart[0:10]), dpi=200)
	plt.close()

