import warnings
import math
import numpy as np

# Plotting Imports
import shapely
import cartopy.crs as ccrs
import cartopy.geodesic
import matplotlib.pyplot as plt

# Brahe Imports
import brahe as bh
import brahe.data_models as bdm
import brahe.access.access as ba
import functions

# sat 1: height = 555000m, inclination = 30 deg
# sat 2: height = 600000m, inclination = 60 deg
sats_initial = [555000, 600000, 30, 60]
coords_initial = functions.get_coords_mult_sats(sats_initial, 0)
alt_initial = functions.get_alt_mult_sats(sats_initial, 0)
# print(coords_initial)
# print(alt_initial)

# propagate by 10 minutes, extract coordinate data every 60 seconds 
end_time = 10 * 60
dt = 60

height_incl = sats_initial
coords = []
coords_sat1 = []
coords_sat2 = []
alt = []
for time in range(0, end_time + 60, dt):
    coords_current = functions.get_coords_mult_sats(height_incl, time)
    coords_sat1.append(coords_current[0])
    coords_sat2.append(coords_current[1])
    coords.append(coords_current)
    alt_current = functions.get_alt_mult_sats(height_incl, time)
    alt.append(alt_current)
    # update height and inclination 
    # height_incl = [alt_current[0], alt_current[1], 30, 60]

# print(coords_sat1)
# print(coords_sat2)
# print(alt)

# Create the figure
fig = plt.figure(figsize=(10,5))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_global()
ax.stock_img()
c = 'b' # Set the plot color

it1 = len(coords_sat1)
for i in range(it1):
    sat1_long = coords_sat1[i][0]
    sat1_lat = coords_sat1[i][1]
    ax.plot(sat1_long, sat1_lat, color=c, marker='o', markersize=3, transform=ccrs.Geodetic())

it2 = len(coords_sat2)
for i in range(it2):
    sat2_long = coords_sat2[i][0]
    sat2_lat = coords_sat2[i][1]
    ax.plot(sat2_long, sat2_lat, color=c, marker='o', markersize=3, transform=ccrs.Geodetic())
plt.show()

elevation_min = 20.0
# for i in len(alt):
    # sat 1
    # lam1 = functions.compute_earth_interior_angle(ele=elevation_min, alt=alt[i][1])
    # lam2 = functions.compute_earth_interior_angle(ele=elevation_min, alt=alt[i][2])


   

