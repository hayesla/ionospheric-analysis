import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
import glob
import datetime
from matplotlib import dates
from astropy import units as u

# all the files.
files = glob.glob("/Users/laurahayes/ionospheric_work/vlf_data_all_birr/sid_alll/*.csv")
files.sort()

# get all the dates from the files.
dates_files = [datetime.datetime.strptime(f.split("/")[-1][0:8], "%Y%m%d") for f in files]
# get unique dates - some days have multiple files.
unique_dates = list(set(dates_files))
unique_dates.sort()

# find all the days between the start and end date of available days to 
# allow us to populate nans for days with no data.
all_days = []
day1 = unique_dates[0]
all_days.append(day1)
while day1 <= unique_dates[-1]:
	day1 = day1 + datetime.timedelta(days=1)
	all_days.append(day1)


def find_files_by_date(date):
	"""
	Find the files available for a particular date.
	Some days have multiple files that need to be 
	concatenated - this returns a list of files
	for that date. 

	Parameters
	----------
	date : ~ datetime.datetime

	Returns
	-------
	file_list : `list`

	"""
	date_str = date.strftime('%Y%m%d')
	file_list = []
	for f in files:
		if date_str in f:
			file_list.append(f)
	return file_list

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
		new_df.reset_index(drop=True)
		return new_df

def make_empty_df(date):
	"""
	Returns an dataframe of all times of the day for a particular 
	date with nans for the 'volts'.

	Parameters
	----------
	date : `datetime.datetime`

	Returns
	-------
	pd.DataFrame

	"""
	dti = pd.date_range(date.strftime('%Y-%m-%d'), periods=17280, freq='5S')
	df = pd.DataFrame({'time':np.array(dti.astype('str')), 'volts':empty_arr})
	return df



