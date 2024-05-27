import warnings
import math
import numpy as np
import csv

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
# sat 3: height = 500000m, inclination = 45 deg
# sat 4: height = 530000m, inclination = 35 deg
sats_initial = [555000, 600000, 500000, 530000, 30, 60, 45, 35] 
coords_initial = functions.get_coords_mult_sats(sats_initial, 0) # find intial coordinates of each satellite 
alt_initial = functions.get_alt_mult_sats(sats_initial, 0)       # find initial altitude of each satellite 


num_sats = int(len(sats_initial) / 2) # define number of satellites in problem 
coord_info = {}                       # define emtpy dictionary to carry through coordinate information 
alt_info = {}                         # define emtpy dictionary to carry through altitude information 
for i in range(num_sats):
    coord_info[f'coords_sat_{i}'] = []
    alt_info[f'alt_sat_{i}'] = []

end_time = 60 * 60      # propagate by 60 minutes 
dt = 60                 # extract coordinate data every 60 seconds 

height_incl = sats_initial
coords = []             # empty list for coordinate information 
alt = []                # empty list for altitude information 
for time in range(0, end_time + 60, dt):
    coords_current = functions.get_coords_mult_sats(height_incl, time)  # get coord info for current time step
    alt_current = functions.get_alt_mult_sats(height_incl, time)        # get alt info for current time step
    
    coords_row = []   # define empty array for csv writing 

    for i in range(num_sats):
        coord_info[f'coords_sat_{i}'].append(coords_current[i])
        alt_info[f'alt_sat_{i}'].append(alt_current[i])
        coords_row.extend(coords_current[i])            # "flatten" lat/long for each satellite into the same row  
    
    coords.append(coords_row)                           # keep track of each satellites lat/long at specific time step 

# Write the collected coordinates to a CSV file
with open('satellite_coordinates.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    
    # Write header
    header = []
    for i in range(num_sats):
        header.append(f'lat_sat_{i}')
        header.append(f'lon_sat_{i}')
    csvwriter.writerow(header)
    
    # Write data rows
    for row in coords:
        csvwriter.writerow(row)
    

'''
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


   

'''