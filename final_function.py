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
alt_sats = [550000 for _ in range(n)]
inclination = [i*(n/360) for i in range(361)] 
sats_initial = []
sats_initial.extend(alt_sats)
sats_initial.extend(inclination)
original_coverage = MultiPolygon()
coords = functions.calculate_lat_lon(sats_initial)
original_area = functions.compute_area(coords) #finds the original area of the satellite 

missing = 3 # the third satellite goes missing
del alt_sats[missing]
del inclination[missing]
input_vec = []
input_vec.extend(alt_sats) 
input_vec.extend(inclination)

new_coords = functions.calculate_lat_lon(input_vec)

# OPTIMIZE find_difference(new_coords) 



# optimize the area_coverage objective function 
def find_difference(coords):
    constellation = functions.compute_area(coords)
    area_overlap = Polygon()
    for sat in constellation: 
        area_overlap = area_overlap.union(sat)
    for x in original_area: 
        area_overlap = area_overlap.union(x)
    return area_overlap

#plot_coverage(constellation, "test_coverage.png")