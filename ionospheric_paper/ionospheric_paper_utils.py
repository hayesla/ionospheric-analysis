import pandas as pd 
import datetime



def get_flarelist_for_vlf():
# flare list
	flare_list = pd.read_csv("/Users/laurahayes/ml_project_flares/flare_analysis/goes_flare_list/final_flare_list.csv")
	flare_list["event_starttime"] = pd.to_datetime(flare_list["event_starttime"])
	flare_list["event_end"] = pd.to_datetime(flare_list["event_endtime"])
	flare_list["event_peaktime"] = pd.to_datetime(flare_list["event_peaktime"])
	flare_list.loc[(pd.to_datetime(flare_list.event_endtime)<flare_list.event_starttime),'event_endtime']=pd.to_datetime(flare_list["event_endtime"]) + datetime.timedelta(days=1)
	flare_list.loc[(pd.to_datetime(flare_list.event_peaktime)<flare_list.event_starttime),'event_peaktime']=pd.to_datetime(flare_list["event_peaktime"]) + datetime.timedelta(days=1)
	flare_list["peak_hour"] = flare_list.event_peaktime.dt.hour
	flare_list["unique_day"] = flare_list.event_peaktime.dt.strftime("%Y-%m-%d")

	# trim between times for which data is available
	flare_list_vlf = flare_list[(flare_list.event_starttime>="2012-08-22")& \
								(flare_list.event_starttime<="2018-04-19")]
	# trim from 8am to 8pm
	flare_list_vlf = flare_list[(flare_list.peak_hour>=8)&(flare_list.peak_hour<=20)]
	flare_list_vlf.reset_index(inplace=True, drop=True)
	flare_list_vlf.to_csv("flare_list_for_vlf.csv", index_label=False)