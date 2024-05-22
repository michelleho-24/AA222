from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pyplot as plt
import cartopy.geodesic
import numpy as np
import math
import shapely
import cartopy.crs as ccrs
import cartopy.geodesic
import matplotlib.pyplot as plt

def compute_earth_interior_angle(ele=0.0, alt=525):
    '''This function computes the earth interior angle for a given elevation angle and altitude.
    This is the angle between the satellite and the observer on the Earth's surface viewing the satellite
    at the given elevation angle. This is more useful for plotting than the look angle, since it
    can used alonside Earth's Radius to draw a circle around the subsatellite point to get
    the view cone of observers on the Earth's surface that would be able to see the satellite.

    Args:
    - ele (float): Elevation angle of the satellite [deg]
    - alt (float): Altitude of the satellite [km]

    '''
    ele = ele * math.pi / 180.0

    rho = math.asin(bh.R_EARTH/(bh.R_EARTH + alt * 1e3))

    eta = math.asin(math.cos(ele)*math.sin(rho))
    lam = math.pi/2.0 - eta - ele

    return lam

# Input latitude and longitude matrix (dummy for now)
lat_lon = np.array([[1,2,3,4], [5,6,7,8],[9,10,11,12], [13, 14, 15, 16]])

# Assume an input of latitude and longitude input as a matrix where columnes correspond to latitude and longitude for each satellite 
num_points = lat_lon.shape[0]
num_sat = lat_lon.shape[1]/2
heights = np.array([500,500,500,500])
ele = 0.0

coverage = MultiPolygon()
original_coverage = MultiPolygon()

# could calculate the original surface area coverage 
# and then do a union of that with the new area coverage 
# and that would give the best matching to the previous area coverage
# better way to return to the original area coverage because you can 
# enforce area coverage of the perious constellation 

# While running through every time step
for i in range(0, num_sat): 
    x = i*2
    latitude = lat_lon[:,x] #index the i-th column to get current latitude
    longitude = lat_lon[:,x+1] #index the i+1-th column to get current longitude 
    alt = heights[i]
    d_horizon_deg = compute_earth_interior_angle(ele,alt)
    circle_points = cartopy.geodesic.Geodesic().circle(lon=longitude, lat=latitude, radius=d_horizon_deg, n_samples=100, endpoint = False)
    coverage_area = Polygon(circle_points) 
    coverage = coverage.union(coverage_area)

fig, ax = plt.subplots(figsize=(10,5))

# Plot coverage area for satellite constellation
for polygon in coverage: 
    ax.add_patch(Circle((polygon.centroid.x, polygon.centroid.y), d_horizon_deg, color = 'b', alpha=0.1))
    ax.set_xlabel('Longitude (degrees)')
    ax.set_ylabel('Latitude (degrees)')
    ax.set_title('Area visible from satellites over one hour') 
    ax.grid(True) 
    plt.show()

total_area = original_coverage.union(coverage) # find area coverage relative to original area