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
	ax.plot(vlf_stats["peak_flare_gl"], (vlf_stats["abs_vlf_db"]), marker='.', ls='')
	for axis in [ax.xaxis, ax.yaxis]:
	    axis.set_major_formatter(ScalarFormatter())
	ax.set_xscale("log")
	ax.set_yscale("log")

