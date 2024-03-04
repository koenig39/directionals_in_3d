import math
import random
import csv
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
api_10=3305306144


def get_elevation(api_10):
    # Open the CSV file for reading
    with open('elevations.csv', 'r') as file:
        csv_reader = csv.reader(file)
        
        # Skip the header row
        next(csv_reader)
        
        # Iterate through each row in the CSV
        for row in csv_reader:
            # Check if the API number in the row starts with the given api_10
            # The API number is expected to be in the 5th column (index 4)
            if row[4].startswith(api_10):
                # If a match is found, return the value from the Latitude column (index 52)
                print(row[52])
                return row[52]  # Adjust the index if necessary based on the actual CSV structure
                
    # Return None if no matching API number is found
    return None


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
        SELECT api_10, uwi, wellbore, md, latitude, longitude, tvd 
        FROM `{project_id}.{dataset_id}.{table_id}`
        WHERE uwi BETWEEN {api10}0000 AND {api10}9999 
        AND tvd IS NOT NULL 
        AND ( wellbore LIKE "LAT%" OR wellbore LIKE "STK%" )
        ORDER BY md ASC
    """

    elevation = int(get_elevation(str(api10)))
    print(query)
    query_job = client.query(query)
    results = list(query_job.result())

    longitude = [row.longitude for row in results]
    latitude = [row.latitude for row in results]
    tvd = [(-1*row.tvd + elevation) for row in results] 
    return latitude, longitude, tvd

def find_min_max_expanded(latitude, longitude, expansion_distance=2500):
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
        where uwi BETWEEN {api10}0000 and {api10}9999 
        AND tvd is not null 
        order by md ASC
    """
    print(query)

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
        md.append(row.tvd)

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
        ax.scatter(wellbore_longs, wellbore_lats, wellbore_md, c='red', marker='o', label='Wellbore')

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
            grid_x, grid_y = np.mgrid[min(surface['longs']):max(surface['longs']):100j, min(surface['lats']):max(surface['lats']):100j]
            
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
        # ax.plot(wellbore_longs, wellbore_lats, wellbore_md, c='black', label='Wellbore', alpha=0.7)

    # Add labels and legend
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Depth/Value')
    plt.legend()

    plt.show()

def get_surface_data(api_key, surface_codes):
    surfaces = []
    # colors = ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Grays', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_grey', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gist_yerg', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'grey', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r']
    # colors = [  'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            # 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            # 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
    # colors = ['viridis', 'plasma', 'inferno', 'magma', 'cividis']
    colors = ['Pastel1', 'Pastel2', 'Paired', 'Accent',
            'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c']
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

# api_10 = 3305306144
# api_10 = 3306101354
# api_10 = 3305304974
# api_10 = 3305304840
# api_10 = 3310501765
# api_10 = 3302500638 #surface_codes = ["TF2","TF3",]
# api_10 = 3305307416

# api_10 = 3310503936
# api_10 = 3300701195
api_10 = 3306104075
# plot_wellbores(3305303948)


# surf_and_wellbore(api_10)

# surface_codes = ["UB","MB","LB","PRNG","TF1","TF2","TF3"]
# surface_codes = ["LB","PRNG","TF1", "TF2"]
# surface_codes = ["UB","MB","LB","TF1","TF2","TF3"]
surface_codes = ["UB","TF1"]
# surface_codes = ["UB","MB","LB"]
# surface_codes = []

# 3305309995	TF1
# 3305309982	MB
# 3305309983	TF1
# 3305309984	MB
# 3305309985	TF1


surfaces = get_surface_data(api_10, surface_codes)
plot_surfaces_and_wellbores([api_10],surfaces=surfaces)
