import numpy as np
import csv
import math
import random
import matplotlib.pyplot as plt
from scipy.spatial import QhullError
from scipy.interpolate import griddata
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import transform
import pyproj
from functools import partial
from google.cloud import bigquery

# Replace 'your_project_id' with your actual Google Cloud project ID
# Replace 'your_dataset' and 'your_table' with your dataset and table names
project_id = 'siros-tech'
dataset_id = 'develop'
table_id = 'lg_directionals_24_02'
client = bigquery.Client(project=project_id)
api_10=3305306144



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



def get_scaled_polygon(poly : Polygon, scale=1):
  minx, miny, maxx, maxy = poly.bounds
  
  width = maxx - minx
  height = maxy - miny

#   if width > height:
#     # Width is larger, scale height
#     height = width * (2/3) 
#   elif height > width:
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


def plot_polygons_areas(original, scaled):
    fig, ax = plt.subplots()
    x, y = original.exterior.xy

    ax.plot(x, y, color='blue', label=f'Original {calculate_wgs84_polygon_area_in_acres(original)} acres')

    x, y = scaled.exterior.xy
    ax.plot(x, y, color='red', label=f'Scaled ({calculate_wgs84_polygon_area_in_acres(scaled)} acres)')

    ax.set_title("Original and Scaled Polygons")
    ax.legend()

    ax.set_aspect('equal')
    plt.show()


def get_elevation(api_10):
    # print(api_10)
    # Open the CSV file for reading
    with open('elevations_38k.csv', 'r') as file:
        csv_reader = csv.reader(file)
        
        # Skip the header row
        next(csv_reader)
        
        # Iterate through each row in the CSV
        for row in csv_reader:
            # print(row)
            # Check if the API number in the row starts with the given api_10
            # The API number is expected to be in the 5th column (index 4)

            #row[0] - API10
            if row[2] and row[0].startswith(str(api_10)):
                # If a match is found, return the value from the Latitude column (index 52)
                print(f"Elevation: {row[1]}, KB elev: {row[2]}")
                return int(row[2])  # Adjust the index if necessary based on the actual CSV structure
            # else:
                # print("asdasd")
        print(f"No KB elevation found for: {api_10}")
    # Return None if no matching API number is found
    return False 


def  get_formation(api10):
    query = f"""
                SELECT siros_formation FROM `siros-tech.lg_well_data.lg_well_data_12_06_23_origin_for_calc` 
                where uwi like ("{api10}%") limit 1
            """
    query_job = client.query(query)
    results = query_job.result()
    for row in results:
        print(f"Formation: {row[0]}")
        return row[0] # formation code example: TF1
    return False
    

def get_wellbores(api10:str):
    query = f"""
        SELECT api_10, uwi, wellbore, md, latitude, longitude, tvd 
        FROM `{project_id}.{dataset_id}.{table_id}`
        WHERE uwi BETWEEN {api10}0000 AND {api10}9999 
        AND tvd IS NOT NULL 
        AND ( wellbore LIKE "LAT%" OR wellbore LIKE "STK%" )
        ORDER BY md ASC
    """
    kb_elevation =get_elevation(str(api10))
    if kb_elevation:
        # print(query)
        elevation = int(get_elevation(str(api10)))
        query_job = client.query(query)
        results = list(query_job.result())

        longitude = [row.longitude for row in results]
        latitude = [row.latitude for row in results]
        tvd = [(row.tvd - kb_elevation) for row in results] 
        # print(latitude, longitude, tvd)
        return latitude, longitude, tvd
    return False


def find_min_max_expanded(latitude, longitude, expansion_distance=13500):
    """
    Expand longitude and latitude min/max values by a specified distance in meters.
    
    Parameters:
    - longitude: List of longitude values.
    - latitude: List of latitude values.
    - expansion_distance: Distance by which to expand the min and max values, in meters.
    
    Returns:
    Tuple of expanded min and max values for both longitude and latitude: 
    (min_longitude, max_longitude, min_latitude, max_latitude)
    """
    # Original min and max
    min_longitude = min(longitude)
    max_longitude = max(longitude)
    min_latitude = min(latitude)
    max_latitude = max(latitude)
    
    # Approximate conversion from distance to degrees
    lat_expansion_deg = expansion_distance / 111000  # 111,000 meters per degree of latitude
    avg_lat = (min_latitude + max_latitude) / 2
    lon_expansion_deg = expansion_distance / (111000 * math.cos(math.radians(avg_lat)))  # Adjust for longitude
    
    # Apply the expansion
    expanded_min_longitude = min_longitude - lon_expansion_deg
    expanded_max_longitude = max_longitude + lon_expansion_deg
    expanded_min_latitude = min_latitude - lat_expansion_deg
    expanded_max_latitude = max_latitude + lat_expansion_deg

    # Return the expanded values directly
    return expanded_min_longitude, expanded_max_longitude, expanded_min_latitude, expanded_max_latitude


def get_geo(api10 : int, formation : str, map_type : str ="STRUCTURE"):
    latitude, longitude, md = get_wellbores(api10)
    min_long, max_long, min_lat, max_lat = find_min_max_expanded(latitude, longitude)
    
    if not formation:
        formation = get_formation(api10)
        print("Formation not found for the given API10.")
        return [], [], []  # Return empty lists for lats, longs, vals

    # Updated query to filter by latitude and longitude boundaries
    query = f"""
        SELECT *
        FROM `siros-tech.develop.das_geology_10_24_23`
        WHERE formation = '{formation}' AND map_type = '{map_type}'
        AND lat BETWEEN {min_lat} AND {max_lat}
        AND long BETWEEN {min_long} AND {max_long}
    """
    # print(query)
    
    query_job = client.query(query)
    results = list(query_job.result())

    lats = [row.lat for row in results]
    longs = [row.long for row in results]
    vals = [row.value for row in results]

    return lats, longs, vals


def plot_surf(lats, longs, vals):
    """
    Plot a surface plot of the given latitude, longitude, and values.
    
    Parameters:
    - lats: Array of latitude values.
    - longs: Array of longitude values.
    - vals: Array of values corresponding to each lat/long pair.
    """
    # Create a grid to interpolate onto
    grid_x, grid_y = np.mgrid[min(longs):max(longs):100j, min(lats):max(lats):100j]
    
    # Interpolate values onto the grid
    grid_z = griddata((longs, lats), vals, (grid_x, grid_y), method='cubic')
    
    # Create the plot
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the surface
    surf = ax.plot_surface(grid_x, grid_y, grid_z, cmap='viridis', edgecolor='none')
    
    # Add labels and show
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Value')
    ax.set_title('Surface Plot')
    plt.colorbar(surf, ax=ax, shrink=0.5, aspect=5)  # Add a color bar to a plot
    plt.show()


def plot_wellbores(api10):
    query = f"""
        select api_10, uwi, wellbore,md, latitude, longitude , tvd
        from `{project_id}.{dataset_id}.{table_id}`
        where uwi BETWEEN {api_10}0000 and {api_10}9999 
        AND tvd IS NOT NULL 
        AND ( wellbore LIKE "LAT%" OR wellbore LIKE "STK%" )
        order by md ASC
    """
    # print(query)

    query_job = client.query(query)

    # Fetch results
    results = query_job.result()

    # Prepare data for plotting
    longitude = []
    latitude = []
    tvd = []

    for row in results:
        longitude.append(row.longitude)
        latitude.append(row.latitude)
        tvd.append(row.tvd)

    # Plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(longitude, latitude, tvd, c='r', marker='o')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('TVD')
    ax.set_zlim(max(tvd), min(tvd))

    plt.show()

#doesn't support area expansion
def surf_and_wellbore(api10):
    # Get wellbore data
    wellbore_lats, wellbore_longs, wellbore_md = get_wellbores(api10)
    wellbore_tvd = tuple(-1 * x for x in wellbore_md)
    # Ensure there is wellbore data to proceed
    if not wellbore_lats or not wellbore_longs:
        print("No wellbore data found.")
        return

    # Get isopach data based on the wellbore data
    formation=get_formation(api10)
    isopach_lats, isopach_longs, isopach_vals = get_geo(api10,formation=formation)

    if not isopach_lats or not isopach_longs:
        print("No isopach data found.")
        return

    try:
        # # Create a grid to interpolate onto
        grid_x, grid_y = np.mgrid[min(isopach_longs):max(isopach_longs):100j, min(isopach_lats):max(isopach_lats):100j]

        # # Interpolate isopach values onto the grid using 'cubic' if possible, 'linear' as fallback
        grid_z = griddata((isopach_longs, isopach_lats), isopach_vals, (grid_x, grid_y), method='cubic')

        # # Create the plot
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')

        # # Plot the isopach surface
        surf = ax.plot_surface(grid_x, grid_y, grid_z, cmap='viridis', edgecolor='none', alpha=0.5, label=formation)

        # Overlay wellbore data
        ax.scatter(wellbore_longs, wellbore_lats, wellbore_tvd, c='red', marker='o', label='Wellbore')

        # Add labels and legend
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_zlabel('Depth/Value')
        ax.legend()

        plt.show()

    except Exception as e:
        print(f"Error plotting data: {e}")


def plot_surfaces_and_wellbores(api10s, surfaces):
    """
    Plot multiple surfaces and wellbore data in a 3D plot.

    Parameters:
    - api10: API number for the wellbores.
    - surfaces: A list of dictionaries, where each dictionary represents a surface dataset
      with keys 'lats', 'longs', 'vals', and 'cmap' for the colormap.
    """
    

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot each surface
    
    for surface in surfaces:
        try:
            # Create a grid for the surface
            grid_x, grid_y = np.mgrid[min(surface['longs']):max(surface['longs']):5j, min(surface['lats']):max(surface['lats']):5j]
            
            # Try cubic interpolation
            grid_z = griddata((surface['longs'], surface['lats']), surface['vals'], (grid_x, grid_y), method='cubic')
        except QhullError:
            # Fallback to linear interpolation if cubic fails
            print("Cubic interpolation failed due to insufficient points, falling back to linear.")
            grid_z = griddata((surface['longs'], surface['lats']), surface['vals'], (grid_x, grid_y), method='linear')

        # Plot the surface
        ax.plot_surface(grid_x, grid_y, grid_z, cmap=surface['cmap'], edgecolor='none', alpha=0.5,  label=surface['formation'])
    # Overlay wellbore data
        
    for api10 in api10s:
        # Get wellbore data
        wellbore_lats, wellbore_longs, wellbore_md = get_wellbores(api10)
        ax.scatter(wellbore_longs, wellbore_lats, wellbore_md, c='brown', s=3, marker='o', label='Wellbore', linewidths=0, alpha=0.8)
        ax.plot(wellbore_longs, wellbore_lats, wellbore_md, c='black', label='Wellbore', alpha=0.7)

    # Add labels and legend
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Depth/Value')
    plt.legend()

    plt.show()

def get_surface_data(api_key, surface_codes):
    surfaces = []
    # colors = ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Grays', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_grey', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gist_yerg', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'grey', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r']
    colors = [  'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
    # colors = ['viridis', 'plasma', 'inferno', 'magma', 'cividis']
    # colors = ['Pastel1', 'Pastel2', 'Paired', 'Accent',
    #         'Dark2', 'Set1', 'Set2', 'Set3',
    #         'tab10', 'tab20', 'tab20b', 'tab20c']
    for code in surface_codes:
        clr=colors[random.randint(0,len(colors)-1)]
        lats, longs, vals = get_geo(api_key, code)
        surface = {
            'lats': lats,
            'longs': longs,
            'vals': vals,
            'cmap': clr,
            'formation': code
        }
        surfaces.append(surface)
    return surfaces

def get_wellbores_from_polygon(polygon) :
    query = f"""
        SELECT dir.api_10, dir.longitude, dir.latitude, dir.md, dir.tvd, dir.wellbore
        FROM `siros-tech.develop.lg_directionals_24_02` dir
        WHERE ST_CONTAINS(
            ST_GEOGFROMTEXT(
                'POLYGON(({polygon}))'), 
            ST_GEOGPOINT(longitude, latitude)) 
        AND tvd IS NOT NULL 
        AND ( wellbore LIKE "LAT%" OR wellbore LIKE "STK%" )
        # OR wellbore LIKE "DIR" )    
        order by api_10, md
    """
    query_job = client.query(query)
    results = list(query_job.result())


    query2 = f"""
        SELECT distinct dir.api_10
        FROM `siros-tech.develop.lg_directionals_24_02` dir
        WHERE ST_CONTAINS(
            ST_GEOGFROMTEXT(
                'POLYGON(({polygon}))'), 
            ST_GEOGPOINT(longitude, latitude)) 
        AND tvd IS NOT NULL 
        AND ( wellbore LIKE "LAT%" OR wellbore LIKE "STK%" )    
        order by api_10
    """
    query_job2 = client.query(query2)
    results2 = list(query_job2.result())


    elevations = {'3305301403' : 2323}

    for row in results2:
        api_10 = str(row.api_10)  # Ensure api_10 is a string, which will be used as the dictionary key
        kb_elevation = get_elevation(api_10)  # Call your function to get the elevation based on api_10

        # Store the elevation in the dictionary with api_10 as the key
        elevations[api_10] = kb_elevation


    
    lats = [row.latitude for row in results]
    longs = [row.longitude for row in results]
    md = [row.md for row in results]
    tvd = [row.tvd-elevations[str(row.api_10)] for row in results]

    return  longs, lats, tvd

def plot_wellbores_from_polygon(polygon,surfaces) :
    # Plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(lats, longs, tvd, c='r', marker='o', s=1 )

    formtions_list = ', '.join([f"'{elem}'" for elem in surfaces])
    formations_query = f"""
        SELECT  formation, map_type, long, lat, value as depth
        FROM `siros-tech.develop.das_geology_10_24_23` 
        WHERE map_type = 'STRUCTURE' and formation IN ({formtions_list}) and 
        ST_CONTAINS(ST_GEOGFROMTEXT('POLYGON((-103.64803185718695 48.0246192603926,-103.5838470682935 48.02446091270679,-103.58379619616841 47.99326020338409,-103.64816832746658  47.99327398116597,-103.64803185718695 48.0246192603926))'), 
  ST_GEOGPOINT(long, lat))"""
    query_job = client.query(formations_query)
    results = list(query_job.result())

    lats = [row.lat for row in results]
    longs = [row.long for row in results]
    depth = [row.depth for row in results]


    # Assuming you have your depth and corresponding coordinates stored in arrays or lists
    depths = depth  # Your 1D array or list of depth values
    x_coords = lats  # Your 1D array or list of x coordinates
    y_coords = longs  # Your 1D array or list of y coordinates

    # Generate a regular grid for interpolation
    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)
    x_grid, y_grid = np.meshgrid(np.linspace(x_min, x_max, 100), np.linspace(y_min, y_max, 100))

    # Interpolate depth values onto the regular grid
    z_grid = griddata((x_coords, y_coords), depths, (x_grid, y_grid), method='linear')

    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface
    surf = ax.plot_surface(x_grid, y_grid, z_grid, cmap='viridis')


    # ax.plot_surface(lats, longs, np.array(depth), edgecolor='none', alpha=0.5)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('MD')
    ax.set_zlim(max(md), min(md))

    plt.show()

    return lats, longs, tvd


def random_color():
    colors = [  'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
    return colors[random.randint(0,len(colors)-1)]


def dsu(): 
    surfaces = """
        SELECT  formation, map_type, long, lat, value as depth
        FROM `siros-tech.develop.das_geology_10_24_23` 
        WHERE map_type = 'STRUCTURE' and formation IN ("TF1","TF3") and 
        ST_CONTAINS(ST_GEOGFROMTEXT('POLYGON((-103.64803185718695 48.0246192603926,-103.5838470682935 48.02446091270679,-103.58379619616841 47.99326020338409,-103.64816832746658  47.99327398116597,-103.64803185718695 48.0246192603926))'), 
  ST_GEOGPOINT(long, lat))"""

    wellbores = """
        SELECT dir.api_10, dir.longitude, dir.latitude, dir.md, dir.tvd
        FROM `siros-tech.develop.lg_directionals_24_02` dir
        WHERE ST_CONTAINS(
            ST_GEOGFROMTEXT(
                'POLYGON((-103.64803185718695 48.0246192603926,-103.5838470682935 48.02446091270679,-103.58379619616841 47.99326020338409,-103.64816832746658  47.99327398116597,-103.64803185718695 48.0246192603926))'), 
            ST_GEOGPOINT(longitude, latitude)) order by api_10, md
    """

    elevation = 2223
    query_job = client.query(wellbores)
    results = list(query_job.result())
    longitude = [row.longitude for row in results]
    latitude = [row.latitude for row in results]
    tvd = [(row.md - elevation) for row in results] 

    # Formation section
    query_surf = client.query(surfaces)
    results_surf = list(query_surf.result())

    # Initialize lists to hold surface data
    surf_long = []
    surf_lat = []
    surf_depth = []

    # Iterate through the results and store surface data
    for row in results_surf:
        surf_long.append(float(row.long))
        surf_lat.append(float(row.lat))
        surf_depth.append(float(row.depth))

    # Plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    try:
        # Create a grid for the surface
        grid_x, grid_y = np.mgrid[min(surf_long):max(surf_long):25j, min(surf_lat):max(surf_lat):25j]

        # Perform cubic interpolation
        grid_z = griddata((surf_long, surf_lat), surf_depth, (grid_x, grid_y), method='cubic')
    except QhullError:
        # Fallback to linear interpolation if cubic fails
        print("Cubic interpolation failed due to insufficient points, falling back to linear.")
        grid_z = griddata((surf_long, surf_lat), surf_depth, (grid_x, grid_y), method='linear')

    # Plot the surface
    ax.plot_surface(grid_x, grid_y, grid_z, cmap=random_color(), edgecolor='none', alpha=0.5)

    # Overlay wellbore data
    # ax.scatter(longitude, latitude, tvd, c='r', marker='o')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('TVD')
    ax.set_zlim(max(tvd), min(tvd))
    plt.legend()
    plt.show()

    return True

polygon_bq = "-103.64803185718695 48.0246192603926,-103.5838470682935 48.02446091270679,-103.58379619616841 47.99326020338409,-103.64816832746658  47.99327398116597,-103.64803185718695 48.0246192603926"
polygon_shaply = Polygon([
    (-103.64803185718695,48.0246192603926),
    (-103.5838470682935, 48.02446091270679),
    (-103.58379619616841, 47.99326020338409),
    (-103.64816832746658,  47.99327398116597),
    (-103.64803185718695, 48.0246192603926)
    ])
surface_codes = ["UB","PRNG"]
polygon_extended = get_scaled_polygon(polygon_shaply,3)
P_wellbores = get_wellbores_from_polygon(polygon=polygon_bq)
plot_polygons_areas(polygon_shaply, polygon_extended)



def get_structure_data(polygon_extended, surface_codes, wellbores_data):

    client = bigquery.Client()
    polygon_wkt = polygon_extended.wkt
    surface_codes_str = ",".join([f"'{code}'" for code in surface_codes])

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    # Prepare data for plotting wellbores
    wellbore_longs, wellbore_lats, wellbore_tvd = wellbores_data
    wellbore_tvd_negative = [-x for x in wellbore_tvd]

    # Overlay wellbore data as 3D scatter plot
    ax.scatter(wellbore_longs, wellbore_lats, wellbore_tvd_negative, c='red', marker='o', s=2, label='Wellbores')

    proxy_artists = []  # List to hold proxy artists for the legend

    # Loop through each surface code to fetch and plot data separately
    for code in surface_codes:
        # Construct the query for each surface code
        query = f"""
            SELECT formation, map_type, long , lat , value
            FROM `siros-tech.develop.das_geology_10_24_23`
            WHERE formation = '{code}' 
            AND map_type = 'STRUCTURE'
            AND ST_CONTAINS(ST_GEOGFROMTEXT('{polygon_wkt}'), ST_GEOGPOINT(long, lat))
        """
        print(query)
        # Execute the query
        query_job = client.query(query)
        results = query_job.result()
        # Prepare data for plotting structure data for the current formation
        longs, lats, values = zip(*[(row.long, row.lat, row.value) for row in results])

        # Check if data is present for the formation
        if longs and lats and values:
            # Create grid and interpolate values for the current formation's structure data
            grid_x, grid_y = np.meshgrid(np.linspace(min(longs), max(longs), 15),
                                         np.linspace(min(lats), max(lats), 15))
            grid_z = griddata((longs, lats), values, (grid_x, grid_y), method='cubic')

            # Plot 3D surface for the current formation
            color = random_color()
            ax.plot_surface(grid_x, grid_y, grid_z, cmap=color, edgecolor='none', alpha=0.5, label=code)
            # proxy_artists.append(Patch(color=color, label=code))
    

    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Value')
    ax.set_title('Structure Data and Wellbores by Formation')

    # Handling legend for multiple formations might be tricky due to how matplotlib handles legends for 3D surfaces.
    # Consider adding a custom legend or annotating each surface plot if necessary.
    plt.legend()

    plt.show()
# Example call to the function (adjust parameters as needed)
get_structure_data(polygon_extended, surface_codes, P_wellbores)


# get_structure_data(polygon_extended, surface_codes,P_wellbores)


# plot_wellbores(polygon_extended)
# surfaces = get_surface_data(api_10, surface_codes)
# get_wellbores_from_polygon(polygon=polygon_extended,surfaces=surfaces)



#test section
#get_elevation(3306101354)
#get_formation(3306101354)
#get_wellbores(3306101354)

# api_10 = 3306101354
# plot_wellbores(api_10)
# surf_and_wellbore(api_10)
# latitude, longitude, tvd = get_wellbores(api_10)
# min_long, max_long, min_lat, max_lat = find_min_max_expanded(latitude, longitude)
# polygon_original = Polygon([(min_long, min_lat), (min_long, max_lat), (max_long, max_lat), (max_long, min_lat)])
# polygon_scaled = get_scaled_polygon(polygon_original,5)
# plot_polygons_areas(polygon_original, polygon_scaled)


# api_10 = 3305306144
# api_10 = 3305304974
# api_10 = 3305304840
# api_10 = 3310501765
# api_10 = 3302500638 #surface_codes = ["TF2","TF3",]
# api_10 = 3305307416

api_10 = 3310503936
# api_10 = 3300701195
# api_10 = 3306104075
# api_10 = 3310503936
# api_10 = 



# surface_codes = ["UB","MB","LB","PRNG","TF1","TF2","TF3"]
# surface_codes = ["UB","MB","LB","TF1","TF2","TF3"]
# surface_codes = []
# surface_codes = ["UB","MB","LB"]

# api_10 = 3305306389 #TF2
# api_10 = 3305306390 #
surface_codes = ["MB","LB","PRNG"]

# 3305309995	TF1
# 3305309982	MB
# 3305309983	TF1
# 3305309984	MB
# 3305309985	TF1


# plot_surfaces_and_wellbores([api_10],surfaces=surfaces)
