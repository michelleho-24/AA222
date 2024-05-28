# from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pyplot as plt
import numpy as np
import math
import shapely
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import csv
from matplotlib.patches import Circle
from shapely.geometry import Point
import json
from shapely.geometry import shape
import functions 
import brahe as bh

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

R_EARTH = 6278*1e3 # radius of earth in meters
def compute_earth_interior_angle(ele, alt):
    ele = ele * math.pi / 180.0
    rho = math.asin(R_EARTH/(R_EARTH + alt * 1e3))
    eta = math.asin(math.cos(ele)*math.sin(rho))
    lam = math.pi/2.0 - eta - ele
    return lam #returns in radians 

ele = 10.0

def calculate_coverage_area(lat, lon, alt,ele): 
    angle = compute_earth_interior_angle(ele,alt)
    radius = angle*R_EARTH
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
        coverage = calculate_coverage_area(lat[t], lon[t],alt, ele)
        satellite_coverage.append(coverage)

    satellite_overall = Polygon()
    for area in satellite_coverage:
        satellite_overall = satellite_overall.union(area)
    constellation_coverage.append(satellite_overall)

constellation = Polygon()
for sat in constellation_coverage: 
    constellation = constellation.union(sat)

def write_polygon(polygon, filename):
    geojson ={
        "type":"Feature",
        "properties":{},
        "geometry": polygon.__geo_interface__
    }
    with open(filename, "w") as f: 
        json.dump(geojson,f)

def read_polygon(filename):
    with open(filename, "r") as f: 
        geojson = json.load(f)

    polygon = shape(geojson['geometry'])
    return polygon

# area_covered = constellation.area
write_polygon(constellation, "test.geojson")

const = read_polygon("test.geojson")

def plot_coverage(coverage, filename): 
    fig, ax = plt.subplots(figsize=(10,8))
    ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
    ax.stock_img()
    # for polygon in coverage.geoms: 
    ax.plot(*coverage.exterior.xy, color = 'red', transform=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()
    ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Coverage Area')
    plt.grid(True)
    # plt.show()
    plt.savefig(filename, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close()

plot_coverage(const, "test_coverage.png")

def plot_eddy(sat_coords):
    # og_coords = [550, 550, 550, 550, 550, 550, 0.0, 30.0, 60.0, 90.0, 120.0, 150.0]
    # LGBFS_coords = [553.69411486,550,548.15294257,553.69411486,550,30,60,90,120,150]
    # PSO_coords = [605.33456691,624.50644854,593.98916237,591.4515896 ,625.9709508,60.53297263 ,25.3114063 ,120.00667778 ,26.98097552 ,91.31959186]
    # random_coords = [633, 508, 584, 640, 566, 179, 77, 110, 93, 8]

    elevation_min = 10.0
    # Create the figure
    fig = plt.figure(figsize=(10,5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    ax.stock_img()
    c = 'b' # Set the plot color
    # plot the optimized sats 
    for i in range(len(sat_coords)/2):
        lam = compute_earth_interior_angle(ele=elevation_min, alt=sat_coords[i])

        # Get the satellite position from above
        x_ecef = functions.tle.state_ecef(t) # Get the ECEF state one day later
        x_geod = bh.sECEFtoGEOD(x_ecef[0:3], use_degrees=True) # Need the array slice to get only the position
        lon, lat = x_geod[0], x_geod[1]

        # Plot Groundstation Location
        ax.plot(lon, lat, color=c, marker='o', markersize=3, transform=ccrs.Geodetic())

        # Get a bunch of points in a circle space at the the right angle offset from the sub-satellite point to draw the view cone
        circle_points = cartopy.geodesic.Geodesic().circle(lon=lon, lat=lat, radius=lam*bh.R_EARTH, n_samples=100, endpoint=False)
        geom = shapely.geometry.Polygon(circle_points)
        ax.add_geometries((geom,), crs=ccrs.Geodetic(), facecolor=c, alpha=0.5, edgecolor='none', linewidth=0)

        plt.show()

plot_eddy([550, 550, 550, 550, 550, 550, 0.0, 30.0, 60.0, 90.0, 120.0, 150.0])
# could calculate the original surface area coverage 
# and then do a union of that with the new area coverage 
# and that would give the best matching to the previous area coverage
# better way to return to the original area coverage because you can 
# enforce area coverage of the perious constellation 