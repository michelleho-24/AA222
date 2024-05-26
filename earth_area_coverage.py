from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pyplot as plt
import cartopy.geodesic
import numpy as np
import math
import shapely
import cartopy.crs as ccrs
import cartopy.geodesic
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import csv
from matplotlib.patches import Circle
from shapely.geometry import Point

def read_csv_file(file_path):
    lat_lon = []
    with open(file_path, 'r') as csvfile: 
        csvreader = csv.reader(csvfile)
        headers_skipped = False
        for row in csvreader: 
            if not headers_skipped: 
                headers_skipped = True
                continue
            lat_lon.append(row)
    return lat_lon

R_EARTH = 6278*1e3 #in meters
def compute_earth_interior_angle(ele, alt):
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

    rho = math.asin(R_EARTH/(R_EARTH + alt * 1e3))

    eta = math.asin(math.cos(ele)*math.sin(rho))
    lam = math.pi/2.0 - eta - ele

    return lam #returns in radians 

ele = 0.0

def calculate_coverage_area(lat, lon, alt): 
    angle = compute_earth_interior_angle(ele,alt)
    radius = math.tan(angle)
    coverage_area = Point(lon, lat).buffer(radius)
    return coverage_area

# Input latitude and longitude matrix (dummy for now)
file_path = 'satellite_coordinates.csv'
lat_lon_ar = read_csv_file(file_path)
lat_lon = np.array(lat_lon_ar)
#lat_lon = np.array([[1,2,3,4], [5,6,7,8],[9,10,11,12], [13, 14, 15, 16]])

# Assume an input of latitude and longitude input as a matrix where columnes correspond to latitude and longitude for each satellite 
rows, columns = lat_lon.shape
num_points  = rows
num_sat = columns/2
heights = np.array([500,500,500])
ele = 0.0

constellation_coverage = MultiPolygon()
original_coverage = MultiPolygon()

# While running through every time step
constellation_coverage = []
for i in range(0,int(num_sat)-1):
    x = i*2
    satellite_coverage = [] # Initialize satellite coverage
    lat = lat_lon[:,x] # extract latitude 
    n = x+1
    lon = lat_lon[:,n] # extract longitude 
    alt = heights[i] # extract height
    for t in range(0,num_points): # iterate through all of the data points 
        coverage = calculate_coverage_area(lat[t], lon[t],alt)
        satellite_coverage.append(coverage)

    satellite_overall = Polygon()
    for area in satellite_coverage:
        satellite_overall = satellite_overall.union(area)
    constellation_coverage.append(satellite_overall)

constellation = Polygon()
for sat in constellation_coverage: 
    constellation = constellation.union(sat)

area_covered = constellation.area

def plot_coverage(coverage, filename): 
    fig, ax = plt.subplots(figsize=(10,8))
    ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor = 'lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor = 'lightblue')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, facecolor = 'lightblue', alpha=0.5)
    ax.add_feature(cfeature.RIVERS, edgecolor = 'blue')
    #ax.set_aspect('equal', adjustable = 'box')
    ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    for polygon in coverage.geoms: 
        ax.plot(*polygon.exterior.xy, color = 'red', transform=ccrs.PlateCarree())
    
    xmin, ymin, xmax, ymax = coverage.bounds
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Coverage Area')
    # plt.grid(True)
    # plt.show()
    plt.savefig(filename, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close()

plot_coverage(constellation, "test_coverage.png")

'''
for i in range(0, int(num_sat)): 
    sat_coverage = MultiPolygon()
    x = i*2
    latitude = lat_lon[:,x] #index the i-th column to get current latitude
    latitude = latitude.astype(float)
    longitude = lat_lon[:,x+1] #index the i+1-th column to get current longitude 
    longitude = longitude.astype(float)
    alt = heights[i]
    d_horizon_deg = compute_earth_interior_angle(ele,alt)
    lat_rad = np.radians(latitude)
    lon_rad = np.radians(longitude)
    for k in range(0,num_points):
        circle_points = np.empty((100,2))
        for j in range(100):
            angle = 360*j/100 #degrees
            circle_points[j, 0] = longitude[k] + d_horizon_deg*np.cos(angle)
            circle_points[j, 1] = latitude[k] + d_horizon_deg*np.sin(angle)
        #circle_points = cartopy.geodesic.Geodesic().circle(lon=longitude, lat=latitude, radius=d_horizon_deg, n_samples=100, endpoint = False)
        coverage_area = Polygon(circle_points) 
        sat_coverage = sat_coverage.union(coverage_area)
    constellation_coverage = constellation_coverage.union(sat_coverage)

'''

'''
# find area coverage relative to original area
total_area = original_coverage.union(constellation_coverage)

# Plot coverage area for satellite constellation
fig, ax = plt.subplots(figsize=(10,5))
for polygon in coverage.geoms: 
    ax.add_patch(Circle((polygon.centroid.x, polygon.centroid.y), d_horizon_deg, color = 'b', alpha=0.1))
    ax.set_xlabel('Longitude (degrees)')
    ax.set_ylabel('Latitude (degrees)')
    ax.set_title('Area visible from satellites over one hour') 
    ax.grid(True) 
    plt.show()
'''
# could calculate the original surface area coverage 
# and then do a union of that with the new area coverage 
# and that would give the best matching to the previous area coverage
# better way to return to the original area coverage because you can 
# enforce area coverage of the perious constellation 