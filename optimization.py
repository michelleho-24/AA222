# This is where the final function for optimization will live, as well as the initial conditions 
import warnings
import math
import numpy as np
import csv
import random

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

with open('output.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Type","Area", "Original Satellites", "Original Coordinates"])

# def optimize(num_sats, deleted)
# number of original satellites
n = 20

# Initialize heights and altitudes for satellites 
alt_sats_orig = [650 for _ in range(n)]
inclination_orig = np.arange(0, 180, 180/n).tolist()
sats_initial = []
sats_initial.extend(alt_sats_orig)
sats_initial.extend(inclination_orig)

coords = functions.calculate_lat_lon(sats_initial)
original_constellation = functions.compute_area(coords, alt_sats=alt_sats_orig, inclination = inclination_orig) 

print("Original Area: ", original_constellation.area)
print("Original Satellites: ", sats_initial)
# print("Original Coordinates", coords)
functions.write_polygon(original_constellation, "original_constellation")

# lose a satellite 
missing = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] # the first satellite goes missing
alt_sats_remaining = [alt_sats_orig[i] for i in range(len(alt_sats_orig)) if i not in missing]
inclination_remaining = [inclination_orig[i] for i in range(len(inclination_orig)) if i not in missing]

sats_remaining = []
sats_remaining.extend(alt_sats_remaining) 
sats_remaining.extend(inclination_remaining)

coords_remaining = functions.calculate_lat_lon(sats_remaining)
remaining_constellation = functions.compute_area(coords_remaining, alt_sats=alt_sats_remaining, inclination = inclination_remaining)

print("Remaining Area: ", remaining_constellation.area)
print("Percent Difference: ", 100*(original_constellation.area - remaining_constellation.area)/original_constellation.area)

# optimize the area_coverage objective function 
def obj_func(input_vec):
    coords = functions.calculate_lat_lon(input_vec)
    alt_sats = input_vec[:n-1]
    inclination = input_vec[n-1:]
    constellation = functions.compute_area(coords, alt_sats, inclination)
    area_overlap = Polygon()
    # for sat in constellation: 
    #     
    area_overlap = area_overlap.union(constellation)
    # for x in original_area: 
    area_overlap = area_overlap.intersection(original_constellation)
    
    return -area_overlap.area

# initial guess 
x0 = sats_remaining


n = n - len(missing)
# Define lower and upper bounds for alt_sats and inclination
lower_bounds = [500]*(n-1) + [0]*(len(x0)-(n-1))
upper_bounds = [650]*(n-1) + [180]*(len(x0)-(n-1))
bounds = Bounds(lower_bounds, upper_bounds)

# LBFGS-B

for method in ['L-BFGS-B']: 
    for i in range(10):
        lbfgs_sats = minimize(obj_func, x0, bounds=bounds, method=method,
                options={'xatol': 1e-8, 'disp': True})
        if lbfgs_sats.success:
            break

coords_lbfgs = functions.calculate_lat_lon(lbfgs_sats.x)
alt_sats_lbfgs = lbfgs_sats.x[:n-1]
inclination_lbfgs = lbfgs_sats.x[n-1:]
lbfgs_constellation = functions.compute_area(coords_lbfgs, alt_sats=alt_sats_lbfgs, inclination = inclination_lbfgs) 

print("LBFGS-B Satellites", lbfgs_sats.x)
print("LBFGS Area: ", lbfgs_constellation.area)
print("Percent Difference from Original: ", 100*(original_constellation.area - lbfgs_constellation.area)/original_constellation.area)
# print("LBFGS-B Coordinates", coords_lbfgs)
functions.write_polygon(lbfgs_constellation, "lbfgs_constellation")

# particle swarm optimization
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
cost, pso_sats = optimizer.optimize(obj_func, iters=1000)
print("PSO Satellites: ", pso_sats)
coords_pso = functions.calculate_lat_lon(pso_sats)
alt_sats_pso = pso_sats[:n-1]
inclination_pso = pso_sats[n-1:]
pso_constellation = functions.compute_area(coords_pso, alt_sats=alt_sats_pso, inclination = inclination_pso)

print("PSO Area: ", pso_constellation.area)
print("Percent Difference from Original: ", 100*(original_constellation.area - pso_constellation.area)/original_constellation.area)
# print("PSO Coordinates", coords_pso)
functions.write_polygon(pso_constellation, "pso_constellation")

random_sats = []
for i in range(n-1):
    random_sats.append(random.uniform(500, 650))
for i in range(n-1):
    random_sats.append(random.uniform(0, 180))

coords_random = functions.calculate_lat_lon(random_sats)
alt_sats_random = random_sats[:n-1]
inclination_random = random_sats[n-1:]
random_constellation = functions.compute_area(coords_random, alt_sats=alt_sats_random, inclination = inclination_random)

print("Random Satellites: ", random_sats)
print("Random Area: ", random_constellation.area)
print("Percent Difference from Original: ", 100*(original_constellation.area - random_constellation.area)/original_constellation.area)
# print("Random Coordinates", coords_random)
functions.write_polygon(random_constellation, "random_constellation")


with open('coords_random.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(coords_random)