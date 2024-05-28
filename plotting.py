from shapely.geometry import Polygon, MultiPolygon
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
<<<<<<< Updated upstream
import cartopy
=======
>>>>>>> Stashed changes

def read_polygon(filename):
    with open(filename, "r") as f: 
        geojson = json.load(f)

    polygon = shape(geojson['geometry'])
    return polygon

<<<<<<< Updated upstream
def plot_coverage(coverage_original, coverage_optimized, filename): 
    fig, ax = plt.subplots(figsize=(10,8), dpi = 1000)
    ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
    ax.stock_img()
    # for polygon in coverage.geoms: 
    ax.plot(*coverage_original.exterior.xy, color = 'red', transform=ccrs.PlateCarree())
    #ax.plot(*coverage_optimized.exterior.xy, color = 'blue', transform=ccrs.PlateCarree())
=======
def plot_coverage(coverage, filename): 
    fig, ax = plt.subplots(figsize=(10,8))
    ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
    ax.stock_img()
    # for polygon in coverage.geoms: 
    ax.plot(*coverage.exterior.xy, color = 'red', transform=ccrs.PlateCarree())
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
file1 = "area_swarm.geojson"
file2 = "original.geojson"

coverage1 = read_polygon(file1)
print("Coverage 1:", coverage1.area)
# get overage before optimization
coverage2 = read_polygon(file2)
print("Coverage 2:", coverage2.area)

#circle_points = cartopy.geodesic.Geodesic().circle(lon=lon, lat=lat, radius=lam*bh.R_EARTH, n_samples=100, endpoint=False)
#geom = shapely.geometry.Polygon(circle_points)

def plotting(polygon):
    x,y = polygon.exterior.xy
    plt.figure()
    plt.plot(x,y)
    plt.fill(x,y,alpha=0.5)
    plt.grid(True)
    plt.axis('equal')
    plt.show()

plotting(coverage2)

#plot_coverage(coverage1, coverage2,"overlapping.png")
# plot_coverage(coverage2, "original.png")
=======
file1 = "area_overlap.geojson"
file2 = "original.geojson"

coverage1 = read_polygon(file1)
coverage2 = read_polygon(file2)

plot_coverage(coverage1, "overlapping.png")
plot_coverage(coverage2, "original.png")
>>>>>>> Stashed changes
