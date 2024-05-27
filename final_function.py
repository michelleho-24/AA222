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

from scipy.optimize import minimize, Bounds 
import pyswarms as ps 

n = 6
# An n-satellite constellation loses m satellites and must reconfigure

# Initialize heights and altitudes for satellites 
alt_sats = [550 for _ in range(n)]
# print(len(alt_sats))
inclination = np.arange(0, 180, 180/n).tolist()
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

    functions.write_polygon(area_overlap, "area_overlap.csv")
    
    return -area_overlap.area

x0 = input_vec

# Define lower and upper bounds for alt_sats and inclination
lower_bounds = [500]*(n-1) + [0]*(len(x0)-(n-1))

upper_bounds = [650]*(n-1) + [180]*(len(x0)-(n-1))

# # Create bounds object
bounds = Bounds(lower_bounds, upper_bounds)

# for method in ['SLSQP']: # 'COBYLA', 'SLSQP'
for method in ['L-BFGS-B']: 
    for i in range(10):
        res = minimize(find_difference, x0, bounds=bounds, method=method,
                options={'xatol': 1e-8, 'disp': True})
        if res.success:
            print(method, res.x)
            break

# Set up hyperparameters
options = {'c1': 0.5, 'c2': 0.3, 'w':0.9}

# Create bounds
bounds = [(lower, upper) for lower, upper in zip(lower_bounds, upper_bounds)]
max_bound = np.array(upper_bounds)
min_bound = np.array(lower_bounds)
bounds = (min_bound, max_bound)

# Initialize swarm
optimizer = ps.single.GlobalBestPSO(n_particles=1, dimensions=len(x0), options=options, bounds=bounds)

# Perform optimization
cost, pos = optimizer.optimize(find_difference, iters=1000)
print(pos)
# Evaluate the function at the optimal position
optimal_function_value = find_difference(pos)
print(optimal_function_value)

# optimize the area_coverage objective function 
def find_difference_2(input_vec):
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

    functions.write_polygon(area_overlap, "area_swarm.geojson")
    functions.write_polygon(original_area, "original.geojson")

area = find_difference_2(pos)


# print(area)
# functions.plot_coverage(area, "test_coverage.png")
# functions.plot_coverage(original_area, "original_coverage.png")
# print(original_area.area)
# print(area)