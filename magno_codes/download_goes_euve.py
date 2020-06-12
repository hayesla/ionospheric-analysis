import urllib
import pandas as pd 
from sunpy.time import parse_time


save_dir = "/Users/laurahayes/ionospheric_work/ionospheric-analysis/magno_codes/euve_data/"
x_flares = pd.read_csv("x_class_w_files.csv")


def download_euvs():
	base_url = "https://hesperia.gsfc.nasa.gov/goes_euv/euve/"
	url_list = parse_time(x_flares["event_date"]).strftime(base_url+"%Y/g15_euve_%Y%m%d.txt") 

	for url in url_list:
		try:
			urllib.request.urlretrieve(url, save_dir + url.split('/')[-1])
		except:
			print(url)

def download_euve_date(date):
	base_url = "https://hesperia.gsfc.nasa.gov/goes_euv/euve/"

	url = parse_time(date).strftime(base_url+"%Y/g15_euve_%Y%m%d.txt")

	try:
		urllib.request.urlretrieve(url, save_dir + url.split('/')[-1])
	except:
		print("didnt work", url)