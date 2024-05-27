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
        coverage = calculate_coverage_area(lat[t], lon[t],alt, ele)
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

plot_coverage(constellation, "test_coverage.png")

# could calculate the original surface area coverage 
# and then do a union of that with the new area coverage 
# and that would give the best matching to the previous area coverage
# better way to return to the original area coverage because you can 
# enforce area coverage of the perious constellation 