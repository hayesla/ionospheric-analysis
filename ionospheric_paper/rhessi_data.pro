pro get_data_hsi

tstart = "2013-05-22 12:00"
tend = "2013-05-22 14:00"

search_network, /enable

o = hsi_obs_summary()
o-> set, obs_time_interval = [tstart, tend]
d = o->getdata(/corrected, class = 'full_rate')


changes = o->changes()
atten = changes.attenuator_state
help,atten,/st


save, filename = 'rhessi_atten_states_20130522.sav', atten


hsi_times = o->getaxis(/ut)

hsi_counts_0612 = d.countrate[1, 0, *]
hsi_counts_1225 = d.countrate[2, 0, *]
hsi_counts_2550 = d.countrate[3, 0, *]
hsi_counts_50100 = d.countrate[4, 0, *]
hsi_counts_100300 = d.countrate[5, 0, *]


save, filename = 'rhessi_lc_corrected_20130522.sav', hsi_times, hsi_counts_0612, hsi_counts_1225, hsi_counts_2550, hsi_counts_50100, hsi_counts_100300

end