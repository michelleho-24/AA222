from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pyplot as plt
import cartopy.geodesic

coverage = MultiPolygon()

# While running through every time step
for min in minutes: 
    # do emmies calculation thing to get latitude and longitude 
    circle_points = cartopy.geodesic.Geodesic().circle(lon=longitude, lat=latitude, radius=d_horizon_deg, n_samples=100, endpoint = False)
    coverage_area = Polygon(circle_points) 
    coverage = coverage.union(coverage_area)

fig, ax = plt.subplots(figsize=(10,5))
for polygon in coverage: 
    ax.add_patch(Circle((polygon.centroid.x, polygon.centroid.y), d_horizon_deg, color = 'b', alpha=0.1))
    ax.set_xlabel('Longitude (degrees)')
    ax.set_ylabel('Latitude (degrees)')
    ax.set_title('Area visible from satellites over one hour') 
    ax.grid(True) 
    plt.show()





