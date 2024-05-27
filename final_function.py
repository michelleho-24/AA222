# This is where the final function for optimization will live, as well as the initial conditions 
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

from shapely.geometry import Polygon, MultiPolygon
import cartopy.feature as cfeature
import csv
from matplotlib.patches import Circle
from shapely.geometry import Point

n = 6
# An n-satellite constellation loses m satellites and must reconfigure

# Initialize heights and altitudes for satellites 
alt_sats = [550 for _ in range(n)]
# print(len(alt_sats))
inclination = np.arange(0, 360, 360/n).tolist()
sats_initial = []
sats_initial.extend(alt_sats)
sats_initial.extend(inclination)
original_coverage = MultiPolygon()

# alt_sats = sats_initial[:n]
# inclination = sats_initial[n:]
# print(alt_sats) 

coords = functions.calculate_lat_lon(sats_initial)
original_area = functions.compute_area(coords, alt_sats=alt_sats, inclination = inclination) #finds the original area of the satellite 

missing = 0 # the first satellite goes missing
del alt_sats[missing]
del inclination[missing]


# print(alt_sats)


input_vec = []
input_vec.extend(alt_sats) 
input_vec.extend(inclination)

# new_coords = functions.calculate_lat_lon(input_vec)

# OPTIMIZE find_difference(new_coords) 
# Find the difference between the original area and the new area




# optimize the area_coverage objective function 
def find_difference(input_vec):
    coords = functions.calculate_lat_lon(input_vec)
    alt_sats = input_vec[:n-1]
    inclination = input_vec[n-1:]
    constellation = functions.compute_area(coords, alt_sats, inclination)
    area_overlap = Polygon()
    # for sat in constellation: 
    #     
    area_overlap = area_overlap.union(constellation)
    # for x in original_area: 
    area_overlap = area_overlap.intersection(original_area)
    
    return area_overlap



area = find_difference(input_vec)
# print(area)
# functions.plot_coverage(area, "test_coverage.png")
# functions.plot_coverage(original_area, "original_coverage.png")
print(original_area.area)
print(area.area)