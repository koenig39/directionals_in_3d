import math
import numpy as np
from scipy.spatial import QhullError
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from google.cloud import bigquery
from mpl_toolkits.mplot3d import Axes3D



# Replace 'your_project_id' with your actual Google Cloud project ID
# Replace 'your_dataset' and 'your_table' with your dataset and table names
project_id = 'siros-tech'
dataset_id = 'develop'
table_id = 'lg_directionals_24_02'
client = bigquery.Client(project=project_id)
api_10=3305308133


def  get_formation(api10):
    query = f"""
                SELECT siros_formation FROM `siros-tech.lg_well_data.lg_well_data_12_06_23_origin_for_calc` 
                where uwi like ("{api10}%") limit 1
            """
    query_job = client.query(query)
    results = query_job.result()
    for row in results:
        return row[0] # formation code example: TF1
    return False
    

def get_wellbores(api10):
    query = f"""
        SELECT api_10, uwi, wellbore, md, latitude, longitude  
        FROM `{project_id}.{dataset_id}.{table_id}`
        WHERE uwi BETWEEN {api10}0000 AND {api10}9999 
        ORDER BY md ASC
    """
    query_job = client.query(query)
    results = list(query_job.result())

    longitude = [row.longitude for row in results]
    latitude = [row.latitude for row in results]
    md = [-1 * row.md for row in results]

    return latitude, longitude, md

def find_min_max_expanded(latitude, longitude, expansion_distance=2000):
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


def get_geo(api10, formation : str, map_type="STRUCTURE"):
    # formation = get_formation(api10)
    
    # Correct variable name from `api_10` to `api10`
    latitude, longitude, md = get_wellbores(api10)
    min_long, max_long, min_lat, max_lat = find_min_max_expanded(latitude, longitude)
    
    if not formation:
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
    print(query)
    
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
    fig = plt.figure()
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
        select api_10, uwi, wellbore,md, latitude, longitude  
        from `{project_id}.{dataset_id}.{table_id}`
        where uwi BETWEEN {api10}0000 and {api10}9999 
        order by md ASC
    """

    query_job = client.query(query)

    # Fetch results
    results = query_job.result()

    # Prepare data for plotting
    longitude = []
    latitude = []
    md = []

    for row in results:
        longitude.append(row.longitude)
        latitude.append(row.latitude)
        md.append(row.md)

    # Plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(longitude, latitude, md, c='r', marker='o')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('MD')
    ax.set_zlim(max(md), min(md))

    plt.show()

def surf_and_wellbore(api10):
    # Get wellbore data
    wellbore_lats, wellbore_longs, wellbore_md = get_wellbores(api10)

    # Ensure there is wellbore data to proceed
    if not wellbore_lats or not wellbore_longs:
        print("No wellbore data found.")
        return

    # Get isopach data based on the wellbore data
    formation="UB"
    isopach_lats, isopach_longs, isopach_vals = get_geo(api10,formation=formation)

    if not isopach_lats or not isopach_longs:
        print("No isopach data found.")
        return

    try:
        # # Create a grid to interpolate onto
        grid_x, grid_y = np.mgrid[min(isopach_longs):max(isopach_longs):10j, min(isopach_lats):max(isopach_lats):10j]

        # # Interpolate isopach values onto the grid using 'cubic' if possible, 'linear' as fallback
        grid_z = griddata((isopach_longs, isopach_lats), isopach_vals, (grid_x, grid_y), method='cubic')

        # # Create the plot
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')

        # # Plot the isopach surface
        surf = ax.plot_surface(grid_x, grid_y, grid_z, cmap='viridis', edgecolor='none', alpha=0.5, label=formation)

        # Overlay wellbore data
        ax.scatter(wellbore_longs, wellbore_lats, wellbore_md, c='red', marker='o', label='Wellbore')

        # Add labels and legend
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_zlabel('Depth/Value')
        ax.legend()

        plt.show()

    except Exception as e:
        print(f"Error plotting data: {e}")

def plot_surfaces_and_wellbores(api10, surfaces):
    """
    Plot multiple surfaces and wellbore data in a 3D plot.

    Parameters:
    - api10: API number for the wellbores.
    - surfaces: A list of dictionaries, where each dictionary represents a surface dataset
      with keys 'lats', 'longs', 'vals', and 'cmap' for the colormap.
    """
    # Get wellbore data
    wellbore_lats, wellbore_longs, wellbore_md = get_wellbores(api10)

    fig = plt.figure(figsize=(20, 14))
    ax = fig.add_subplot(111, projection='3d')

    # Plot each surface
    for surface in surfaces:
        try:
            # Create a grid for the surface
            grid_x, grid_y = np.mgrid[min(surface['longs']):max(surface['longs']):100j, min(surface['lats']):max(surface['lats']):100j]
            
            # Try cubic interpolation
            grid_z = griddata((surface['longs'], surface['lats']), surface['vals'], (grid_x, grid_y), method='cubic')
        except QhullError:
            # Fallback to linear interpolation if cubic fails
            print("Cubic interpolation failed due to insufficient points, falling back to linear.")
            grid_z = griddata((surface['longs'], surface['lats']), surface['vals'], (grid_x, grid_y), method='linear')

        # Plot the surface
        # ax.plot_surface(grid_x, grid_y, grid_z, cmap=surface['cmap'], edgecolor='none', alpha=0.7)
        ax.scatter(surface['longs'], surface['lats'], surface['vals'], c='red', marker='o', label='Wellbore')

    # Overlay wellbore data
    ax.scatter(wellbore_longs, wellbore_lats, wellbore_md, c='red', marker='o', label='Wellbore')

    # Add labels and legend
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Depth/Value')
    plt.legend()

    plt.show()


# tf_1_lats, tf_1_longs, tf_1_vals = get_geo(api_10,'UB')
# tf_3_lats, tf_3_longs, tf_3_vals = get_geo(api_10,'TF1')
# tf_2_lats, tf_2_longs, tf_2_vals = get_geo(api_10,'TF3')


# surfaces = [
#         {
#         'lats': tf_1_lats,
#         'longs': tf_1_longs,
#         'vals': tf_1_vals,
#         'cmap': 'viridis'
#     },
#             {
#         'lats': tf_2_lats,
#         'longs': tf_2_longs,
#         'vals': tf_2_vals,
#         'cmap': 'viridis'
#     },
#                 {
#         'lats': tf_3_lats,
#         'longs': tf_3_longs,
#         'vals': tf_3_vals,
#         'cmap': 'inferno'
#     },
# ]


# surf_and_wellbore(3310505584)

plot_wellbores(3310505384)
# plot_surfaces_and_wellbores(api_10,surfaces)
# plot_surf(tf_2_lats, tf_2_longs, tf_2_vals)

# plot_surf(tf_1_lats, tf_1_longs, tf_1_vals)
# surf_and_wellbore(api_10)