from shapely.geometry import Polygon

# make a square using the Polygon class
polygon_1 = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
polygon_2 = Polygon([(0.5, 0), (0.5, 1), (1.5, 1), (1.5, 0)])

# find the intersection of the two squares
intersection = polygon_1.intersection(polygon_2)
print(intersection.area) # 1.0

union = polygon_1.union(polygon_2)
print(union.area) # 1.0
