from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from sunpy.time import parse_time
from sunpy.coordinates import frames

#plt.rcParams["font.family"] = "Times New Roman"

fig, ax = plt.subplots(figsize=(6, 4))
# m = Basemap(projection='cyl', llcrnrlat=-20, urcrnrlat=89, 
# 			llcrnrlon=-260, urcrnrlon=-20, resolution='c')


m = Basemap(projection='cyl',llcrnrlat=-30,urcrnrlat=90,\
             llcrnrlon=-120,urcrnrlon=80,resolution='c')
# date = parse_time('2015-11-04 13:30').datetime
date = parse_time('2020-11-29 12:00').datetime

sun = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, 
				obstime=date, frame=frames.HeliographicStonyhurst )    

coords = sun.transform_to('itrs').earth_location.to_geodetic()  
lon, lat = coords.lon.value, coords.lat.value

# m = Basemap(projection='cyl',lat_0=lat, lon_0=lon, resolution='c')
# m = Basemap(projection='cyl',lat_0=0, lon_0=0, resolution='c')


#m.drawlsmask()
m.fillcontinents()
# m.drawcoastlines()
m.drawcountries()
#m.drawstates()

mag_color='tab:green'
vlf_color='tab:blue'


naa_lat, naa_lon = 44.644506, -67.284565
birr_lon, birr_lat = -7.9, 53

esk_lon, esk_lat = -3.10, 55.16
# rht_lat, rht_lon = 41.720, -111.822
# nml_lat, nml_lon = 46.365987, -98.335667
# nlk_lat, nlk_lon = 48.203633,  -121.916828

# kak_lat, kak_lon = 36.1356, (140.1111-180)-180

date = parse_time('2015-11-04 13:30').datetime

m.plot(lon, lat, color='yellow', marker='o', markeredgecolor='k', ms=10)
m.plot(naa_lon, naa_lat, color=vlf_color, marker='o',ms=5)
# m.plot(birr_lon, birr_lat, color=mag_color, marker='o', ms=5)
m.plot(esk_lon, esk_lat, color=mag_color, marker='o', ms=5)
# m.plot(nlk_lon, nlk_lat, color='b', marker='x')
# m.plot(birr_lon, birr_lat, color=vlf_color, marker='o', ms=5)
m.drawgreatcircle(naa_lon,naa_lat,esk_lon,esk_lat,linewidth=2,color=vlf_color)
# m.drawgreatcircle(nlk_lon,nlk_lat,rht_lon,rht_lat,linewidth=2,color='b')
# m.drawgreatcircle(nml_lon,nml_lat,rht_lon,rht_lat,linewidth=2,color='b')

m.nightshade(date, alpha=0.6)
m.drawparallels(np.arange(-30.,91.,30.), labels=[True,True,False,False])
m.drawmeridians(np.arange(-180.,181.,60.), labels=[False,False,False,True])

plt.title(date.strftime('%Y-%m-%d %H:%M'))

# plt.text(birr_lon+5, birr_lat+3, 'Magnetometer', color=mag_color)
# plt.text(rht_lon+5, rht_lat-5, 'VLF', color=vlf_color)


plt.tight_layout()

plt.savefig('test_map2.png', dpi=200)
plt.close()



# import ephem
# import math
# from datetime import datetime
# from sunpy import parse_time 



# sun = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, 
# 				obstime='2020-02-12 19:43', frame='heliographic_stonyhurst')    

# sun.transform_to('itrs').earth_location.to_geodetic()  



# def test_plot():
# 	date_lat, date_lon = [], []
# 	date_new = parse_time('2020-02-12 00:00').datetime
# 	for i in range(1440):
# 		date_new = date+datetime.timedelta(minutes=i)		

# 		sun = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU, 
# 						obstime=date_new, frame='heliographic_stonyhurst')    

# 		coords = sun.transform_to('itrs').earth_location.to_geodetic()  
# 		lon, lat = coords.lon.value, coords.lat.value
# 		date_lat.append(lat)
# 		date_lon.append(lon)

# 	m = Basemap(projection='cyl',lat_0=lat, lon_0=lon, resolution='c')

# 	m.drawlsmask()
# 	m.fillcontinents()
# 	m.drawcoastlines()
# 	m.drawcountries()
# 	for i in range(len(date_lon)):
# 		m.plot(date_lon[i], date_lat[i], color='b', marker='o',ms=5)

# 	plt.show()

# def find_subsolar(date):
# 	greenwich = ephem.Observer()
# 	greenwich.lat = "0"
# 	greenwich.lon = "0"
# 	#greenwich.date = datetime.utcnow()
# 	greenwich.date = date.datetime
# 	sun = ephem.Sun(greenwich)
# 	sun.compute(greenwich.date)
# 	sun_lon = math.degrees(sun.ra - greenwich.sidereal_time() )
# 	if sun_lon < -180.0 :
# 	  sun_lon = 360.0 + sun_lon 
# 	elif sun_lon > 180.0 :
# 	  sun_lon = sun_lon - 360.0
# 	sun_lat = math.degrees(sun.dec)

# 	print("Subsolar Point Sun Lon:",sun_lon, "Lat:",sun_lat)