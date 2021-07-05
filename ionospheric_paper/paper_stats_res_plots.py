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
plt.ion()
from astropy.visualization import hist
import seaborn as sns

# setting plotting defaults
sns.set_context("paper")
plt.rcParams['xtick.direction'] = "in"
plt.rcParams['ytick.direction'] = "in"
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['font.family'] = 'Helvetica'
cmap_paper = "viridis"


""" ------------------------------------------------------------------
This is the script to make the plots for the paper. 
The data read in here is created in the script `final_vlf_stats_df.py`



-----------------------------------------------------------------------"""



### Read in data ###
vlf_stats = pd.read_csv("vlf_stats_results_new.csv")
vlf_stats["event_starttime"] = pd.to_datetime(vlf_stats["event_starttime"])


def flux_vs_amp():
	
	fig, ax = plt.subplots()
	im = ax.scatter(vlf_stats["peak_flare_gl-bg"], vlf_stats["abs_vlf_db"], 
	                c=vlf_stats["peak_flare_gl"], norm=colors.LogNorm(), 
	                cmap=cmap_paper)
	ax.set_xscale("log")
	ax.set_yscale("log")

	ax.set_xlabel("GOES 1-8$\mathrm{\AA}$ X-ray peak (Wm$^{-2}$)")
	ax.set_ylabel("VLF amplitude excess (dB)")

	ax.yaxis.set_major_formatter(ScalarFormatter())
	# cbar = plt.colorbar(im)
	# cbar.set_label("GOES 1-8$\mathrm{\AA}$ X-ray peak (Wm$^{-2}$)")
	cc = scipy.stats.spearmanr(vlf_stats["peak_flare_gl-bg"], vlf_stats["abs_vlf_db"])
	print(cc)
	ax.text(0.85, 0.06, "C.C. {:.2f}".format(round(cc.correlation, 2)), transform=ax.transAxes)
