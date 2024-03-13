import matplotlib.pyplot as plt 
from shapely import MultiPolygon, Polygon, area
import math
from shapely.ops import transform
import pyproj
from functools import partial


def calculate_wgs84_polygon_area_in_acres(coords):
    """
    Calculate the area of a polygon specified in WGS84 coordinates and returns the area in square acres.

    Parameters:
    - coords: A list of tuples, where each tuple represents a point (longitude, latitude).

    Returns:
    - The area of the polygon in square acres.
    """
    # Create a Polygon from the coordinates
    polygon = Polygon(coords)

    # Define the projection transformation: from WGS84 to an equal-area projection
    # Using a commonly used equal area projection (EPSG:6933)
    proj = partial(
        pyproj.transform,
        pyproj.Proj(init='epsg:4326'), # Source projection: WGS84
        pyproj.Proj(init='epsg:6933')  # Target projection: World Equidistant Cylindrical (equal area)
    )

    # Transform the polygon to the new projection and calculate the area
    transformed_polygon = transform(proj, polygon)
    area_in_sq_meters = transformed_polygon.area
    
    # Convert square meters to square acres
    area_in_sq_acres = round( area_in_sq_meters * 0.000247105, 3)
    
    return area_in_sq_acres


def get_scaled_polygon(polygon, scale):

    poly = Polygon(polygon)
    minx, miny, maxx, maxy = poly.bounds

    width = maxx - minx
    height = maxy - miny

    # if width > height:
    #     # Width is larger, scale height
    #     height = width * (2/3) 
    # elif height > width:
    #     # Height is larger, scale width
    #     width = height * (2/3)

    centerx = (minx + maxx) / 2 
    centery = (miny + maxy) / 2

    scaled_width = width * scale
    scaled_height = height * scale

    minx = centerx - scaled_width/2
    maxx = centerx + scaled_width/2
    miny = centery - scaled_height/2
    maxy = centery + scaled_height/2

    return Polygon([(minx, miny), (maxx, miny), (maxx, maxy), (minx, maxy)])


def plot_polygons(original, scaled):
    fig, ax = plt.subplots()
    x, y = original.exterior.xy

    ax.plot(x, y, color='blue', label=f'Original {calculate_wgs84_polygon_area_in_acres(original)} acres')

    x, y = scaled.exterior.xy
    ax.plot(x, y, color='red', label=f'Scaled ({calculate_wgs84_polygon_area_in_acres(scaled)} acres)')

    ax.set_title("Original and Scaled Polygons")
    ax.legend()

    ax.set_aspect('equal')
    plt.show()


polygon = Polygon([(0,0), (-3,-4), (-1,0), (-5,3), (0,5),(3,1)])
polygon2 = Polygon([
    (-103.64803185718695,48.0246192603926),
    (-103.5838470682935, 48.02446091270679),
    (-103.58379619616841, 47.99326020338409),
    (-103.64816832746658,  47.99327398116597),
    (-103.64803185718695, 48.0246192603926)
    ])
scaled = get_scaled_polygon(polygon, 1.2)
plot_polygons(polygon, scaled)