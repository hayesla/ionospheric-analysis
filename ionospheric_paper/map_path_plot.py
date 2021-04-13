import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature.nightshade import Nightshade
import cartopy.feature as cfeature
import datetime
import numpy as np 
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import seaborn as sns

sns.set_context("paper", font_scale=1.3)
# projection=ccrs.PlateCarree()
# #hacked_proj = projection
# projection._threshold /= 20.


point_colors = "blue"
img_colors = False

# there is an issue with the resolution for
# plotting great circle, this is a hack to overwrite
# the threshold.
projection = ccrs.Miller()
class hacked_miller(ccrs.Miller):
    @property
    def threshold(self):
        return 0.05

projection = hacked_miller()

fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=projection)

# setting up tick labels etc
ax.set_xticks(np.arange(0, 360,20) , crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(-90, 90, 10), crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.gridlines()

# add in features of interest on map
ax.set_extent((-90, 10, 30, 65), crs=ccrs.PlateCarree())
ax.coastlines(resolution="50m")
ax.add_feature(cfeature.BORDERS, lw=0.5, color="grey")
ax.add_feature(cfeature.LAND, facecolor="tan")
ax.add_feature(cfeature.OCEAN, facecolor="lightblue")

# hack to allow LAKE feature to have both edgecolor and facecolor
hacked_lake = cfeature.NaturalEarthFeature(cfeature.LAKES.category, cfeature.LAKES.name, cfeature.LAKES.scale, 
									edgecolor="k", lw=0.5)
ax.add_feature(hacked_lake, facecolor="lightblue")

# if True then overplot colors of land/sea from image
if img_colors:
	ax.stock_img()

# plot the locations of transmitter/reciever
naa_lat, naa_lon = 44.644506, -67.284565
birr_lon, birr_lat = -7.9, 53

ax.text(naa_lon+1, naa_lat-1, "NAA", color=point_colors, weight="bold")
ax.text(birr_lon-3, birr_lat, "Birr", color=point_colors, weight="bold")

ax.plot([naa_lon, birr_lon], [naa_lat, birr_lat],
         color=point_colors, linewidth=2, marker='o',
         transform=ccrs.Geodetic())


ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
plt.tight_layout()
plt.savefig("map_of_path.png", dpi=200)
plt.close()