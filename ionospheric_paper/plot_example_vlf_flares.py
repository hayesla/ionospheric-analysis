import matplotlib.pyplot as plt 
from matplotlib import dates
import numpy as np 
import pandas as pd 
from sunpy import timeseries as ts 
import datetime
import sys
sys.path.append("..")
from read_files import sid_to_series
from matplotlib import gridspec
import matplotlib
import seaborn as sns
sns.set_context("paper", font_scale=1.2)


"""------------------------------
Script for figure 2 of VLF paper


An example from 2015-10-02 with
lots of flares.

------------------------------"""

matplotlib.rcParams['xtick.direction'] = "in"
matplotlib.rcParams['ytick.direction'] = "in"
matplotlib.rcParams['xtick.minor.visible'] = True
matplotlib.rcParams['ytick.minor.visible'] = True
plt.rcParams['font.family'] = 'Helvetica'

tstart = "2015-10-02 07:00"
tend = "2015-10-02 22:00"

goes_file = "/Users/laurahayes/QPP/stats_study/TEBBS/goes_rawdata/go1520151002.fits"
goes_file = "sci_gxrs-l2-irrad_g13_d20151002_v0-0-0.nc"
vlf_file = "/Users/laurahayes/ionospheric_work/vlf_data_all_birr/sid_alll/20151002_000000_NAA_S-0055.csv"
vlf_quiet = "/Users/laurahayes/ionospheric_work/vlf_data_all_birr/sid_alll/20151006_000000_NAA_S-0055.csv"

goes_data = ts.TimeSeries(goes_file)
gl, gs = goes_data.to_dataframe()['xrsb'], goes_data.to_dataframe()['xrsa']

sid_data = sid_to_series(vlf_file)
sid_data_quiet = sid_to_series(vlf_quiet)
sid_data_quiet.index = sid_data_quiet.index - datetime.timedelta(days=4)


def calc_amp(x):
	return 20*np.log10(x + 5) - 61 + 107

def make_paper_plot(sid_data, sid_data_quiet, amp=False):
	"""
	Make plot of example flare day for paper
	if amp is True, the vlf amplitude is plotted in db rather than volts
	"""
	if amp:
		sid_data = calc_amp(sid_data)
		sid_data_quiet = calc_amp(sid_data_quiet)

	# fig, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(8, 8))
	fig = plt.figure(figsize=(8, 6))
	grids = gridspec.GridSpec(2, 1, height_ratios=[1, 2]) 
	ax1 = plt.subplot(grids[0])
	ax2 = plt.subplot(grids[1], sharex=ax1)
	ax1.plot(gl, color="tab:red", label="1-8 $\mathrm{\AA}$")
	ax1.plot(gs, color="tab:blue", label="0.5-4 $\mathrm{\AA}$")
	ax1.set_yscale("log")
	ax1.set_ylim(1e-8, 4e-4)
	ax1.set_ylabel("Flux Wm$^{-2}$")
	ax1.tick_params(which="both", labelbottom=False)
	ax1_rhs = ax1.twinx()

	ax1_rhs.set_yscale("log")
	ax1_rhs.set_ylim(1e-8, 4e-4)
	ax1_rhs.set_yticks((1e-8, 1e-7, 1e-6, 1e-5, 1e-4))
	ax1_rhs.set_yticklabels(('A', 'B', 'C', 'M', 'X'))

	ax1.yaxis.grid(True, 'major')
	ax1.xaxis.grid(False, 'major')
	ax1.legend(loc="upper right")

	ax2.plot(sid_data_quiet, color="grey", lw=0.5, label="Quiet day")
	ax2.plot(sid_data, label="VLF amplitude", color="k")
	ax2.legend(loc="upper right")

	ax2.set_xlim(tstart, tend)
	ax2.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
	ax2.yaxis.grid(True, "major")

	if amp:
		ax2.set_ylim(35, 68)
		ax2.set_ylabel("VLF Amplitude (dB)")
	else:
		ax2.set_ylim(-5, 5)
		ax2.set_ylabel("VLF Amplitude (volts)")

	ax2.set_xlabel("Time {:s} (UT)".format(gl.index[100].strftime("%Y-%m-%d")))
	for a in (ax1, ax2):
		a.axvspan("2015-10-02 07:00", "2015-10-02 09:30", color="grey", alpha=0.5, zorder=3)
		a.axvspan("2015-10-02 19:00", "2015-10-02 22:00", color="grey", alpha=0.5, zorder=3)
	ax1.text(0.01, 0.91, "a. Solar X-ray Observations", transform=ax1.transAxes)
	ax2.text(0.01, 0.95, "b. Ionospheric VLF Observations", transform=ax2.transAxes)


	plt.tight_layout()
	plt.subplots_adjust(hspace=0.01)
	if amp:
		plt.savefig("./final_paper_plots/example_active_day.png", dpi=300, bbox_inches="tight")
	else:
		plt.savefig("./final_paper_plots/example_active_day2.png", dpi=300, bbox_inches="tight")
	plt.close()

# make plots 
make_paper_plot(sid_data, sid_data_quiet)
make_paper_plot(sid_data, sid_data_quiet, amp=True)
