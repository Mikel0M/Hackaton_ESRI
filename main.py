"""
Created on Tue Feb  4 10:12:31 2025

@author: mikel
"""
from fastapi import FastAPI, File, UploadFile, HTTPException
import shutil
import ifcopenshell
import ifcopenshell.geom
import numpy as np
import json
import os
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import MultiPoint, Polygon
from pyproj import Transformer

# directory to save uploaded files
UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

app = FastAPI()


"""ifcopenshell functions"""

def get_building_height(model):
    """
    To get simple height from OverallHeight
    """
    heights = []
    for building in model.by_type("IfcBuilding"):
        if hasattr(building, "OverallHeight"):
            heights.append(building.OverallHeight)
    return max(heights) if heights else None


def get_building_height_complex(model):
    """
    To get height from either OverallHeight or estimate it from the storey elevations
    """
    heights = []
    floors = model.by_type("IfcBuildingStorey")

    #first try to gert the height of the IfcBuilding
    for building in model.by_type("IfcBuilding"):
        if hasattr(building, "OverallHeight"):
            heights.append(building.OverallHeight)

    #If no OverallHeight, estimate based on storey elevations
    if not heights:
        storey_heights = [
            storey.Elevation for storey in floors if hasattr(storey, "Elevation")
        ]
        if storey_heights:
            estimated_height = max(storey_heights)
        else:
            estimated_height = "Unknown height"
    else:
        estimated_height = max(heights)

    return estimated_height, len(floors)

def get_floor_area_simple(model):
    """
    Use get_floor_area_first_test if your file structure is simple and uses properties directly on elements 
    (like GrossFloorArea or NetArea on storeys, slabs, or spaces)
    """
    total_area = 0

    # Try getting the area from IfcBuildingStorey (preferred)
    for storey in model.by_type("IfcBuildingStorey"):
        if hasattr(storey, "GrossFloorArea") and storey.GrossFloorArea is not None:
            total_area += storey.GrossFloorArea
    
    # If no GrossFloorArea, sum up areas from IfcSlab
    if total_area == 0:
        for slab in model.by_type("IfcSlab"):
            for relDefinesByProperties in slab.IsDefinedBy:
                if relDefinesByProperties.is_a("IfcRelDefinesByProperties"):
                    propSet = relDefinesByProperties.RelatingPropertyDefinition
                    if propSet.is_a("IfcPropertySet"):
                        for prop in propSet.HasProperties:
                            if prop.Name == "NetArea": # Some models use NetArea
                                total_area += prop.NominalValue.wrappedValue

    # If still no area, try summing up IfcSpace areas
    if total_area == 0:
        for space in model.by_type("IfcSpace"):
            for relDefinesByProperties in space.IsDefinedBy:
                if relDefinesByProperties.is_a("IfcRelDefinesByProperties"):
                    propSet = relDefinesByProperties.RelatingPropertyDefinition
                    if propSet.is_a("IfcpropertySet"):
                        for prop in propSet.HasProperties:
                            if prop.Name == "GrossFloorArea":
                                total_area += prop.NominalValue.wrappedValue

    return total_area if total_area > 0 else "Unknown area"


def save_to_json(data, file_path):
    """
    Save data to a JSON file, ensuring all values are JSON serializable.
    """
    def convert_numpy_types(obj):
        """ Recursively convert NumPy types to native Python types """
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()  # Convert NumPy arrays to lists
        elif isinstance(obj, dict):
            return {key: convert_numpy_types(value) for key, value in obj.items()}
        elif isinstance(obj, list):
            return [convert_numpy_types(item) for item in obj]
        return obj  # Return as is if it's already a native Python type
    
    # Convert NumPy types before saving
    clean_data = convert_numpy_types(data)

    # Save to JSON
    with open(file_path, "w") as json_file:
        json.dump(clean_data, json_file, indent=4)

    print(f"✅ JSON saved successfully to {file_path}")

def get_floor_area(model):
    """Use get_floor_area for more complex files where area might be stored in related properties or element quantities, 
    and if you want more flexibility in handling different ways area data can be structured
    """
    total_area = 0

    # Try getting area from IfcBuildingStorey
    for storey in model.by_type("IfcBuildingStorey"):
        for relDefinesByProperties in storey.IsDefinedBy:
            if relDefinesByProperties.is_a("IfcRelDefinesByProperties"):
                propSet = relDefinesByProperties.RelatingPropertyDefinition
                
                # If area is in IfcPropertySet
                if propSet.is_a("IfcPropertySet"):
                    for prop in propSet.HasProperties:
                        if prop.Name in ["GrossFloorArea", "NetFloorArea"]:  # Update based on debug output
                            total_area += prop.NominalValue.wrappedValue

                # If area is in IfcElementQuantity
                elif propSet.is_a("IfcElementQuantity"):
                    for quantity in propSet.Quantities:
                        if quantity.is_a("IfcQuantityArea"):
                            total_area += quantity.AreaValue

    return total_area if total_area > 0 else "Unknown area"

"""Functions based on Geopandas"""

def generate_building_coords(center_x, center_y, width=12, length=24):
    """
    Generate building coordinates for a rectangular building centered at (center_x, center_y).
    """
    half_width = width / 2
    half_length = length / 2
    
    return [
        (center_x + half_width, center_y + half_length),
        (center_x - half_width, center_y + half_length),
        (center_x - half_width, center_y - half_length),
        (center_x + half_width, center_y - half_length)
    ]

def check_zoning(building_polygon, zoning_map):
    """
    Check if a building is within any zoning area.
    """
    within_zoning = zoning_map.contains(building_polygon)
    return zoning_map[within_zoning]

import matplotlib.pyplot as plt
from shapely.geometry import Point

def plot_building_and_zoning(building_polygon, zoning_map, lv95_coords, vertices, zoom_factor=1):
    """
    Plot the zoning map, the building polygon, and an additional polygon from vertices.
    Also zooms into the area of both polygons and adds a red point at the building's center.
    
    :param building_polygon: The building polygon to plot.
    :param zoning_map: The zoning map to plot as the background.
    :param lv95_coords: Coordinates to add a red point (typically the centroid or center of the building).
    :param vertices: A list of vertices for an additional polygon to plot.
    :param zoom_factor: Factor by which to zoom in on the building and vertices (default 1 means no zoom).
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    zoning_map.plot(ax=ax, color='lightgray', edgecolor='black')

    # Plot the red point at the lv95_coords (center of the building)
    ax.plot(lv95_coords[0], lv95_coords[1], 'ro', markersize=8, label="Building Center")

    # Create the building polygon and plot it (if provided)
    if building_polygon:
        x_building, y_building = building_polygon.exterior.xy
        ax.fill(x_building, y_building, color='red', alpha=0.5, label="Building Polygon")
    
    # If vertices are provided, plot the additional polygon
    if vertices:
        # Create the polygon from the vertices
        poly = Polygon(vertices)
        x_poly, y_poly = poly.exterior.xy
        ax.fill(x_poly, y_poly, color='blue', alpha=0.5, label="Additional Polygon")

    # Calculate the bounding box of both the building polygon and the additional polygon
    minx_building, miny_building, maxx_building, maxy_building = building_polygon.bounds
    if vertices:
        # Get the bounding box of the additional polygon
        poly = Polygon(vertices)
        minx_vertices, miny_vertices, maxx_vertices, maxy_vertices = poly.bounds
    else:
        minx_vertices, miny_vertices, maxx_vertices, maxy_vertices = minx_building, miny_building, maxx_building, maxy_building

    # Combine the bounding boxes of both the building and additional polygon
    minx = min(minx_building, minx_vertices)
    miny = min(miny_building, miny_vertices)
    maxx = max(maxx_building, maxx_vertices)
    maxy = max(maxy_building, maxy_vertices)

    # Apply zoom based on the combined bounding box of both polygons
    x_margin = (maxx - minx) * zoom_factor
    y_margin = (maxy - miny) * zoom_factor

    # Set the limits of the plot to zoom into both the building and additional polygon
    ax.set_xlim(minx - x_margin, maxx + x_margin)
    ax.set_ylim(miny - y_margin, maxy + y_margin)

    plt.legend()
    plt.show()


def extract_zoning_restrictions(zoning_with_building):
    """
    Extract zoning restrictions for the building and store them in a dictionary.
    """
    zoning_restrictions = {}

    if not zoning_with_building.empty:
        zoning_restrictions = {
            "Max Building Height": zoning_with_building.get('GEBAEUDEHO', -99).values[0],
            "Total Allowed Height": zoning_with_building.get('GESAMTHOEH', -99).values[0],
            "Max Full Floors": zoning_with_building.get('VOLLGESCHO', 0).values[0],
            "Allowed Attic Floors": zoning_with_building.get('DACHGESCHO', 0).values[0],
            "Floor Area Ratio (FAR)": zoning_with_building.get('AUSNUETZUN', -99).values[0]
        }
    else:
        print("🚨 Building is not within any zoning area!")

    return zoning_restrictions


"""Here is for compdealing with coordinates"""

def dms_to_decimal(degrees, minutes, seconds, microseconds=0):
    """
    Convert Degrees, Minutes, Seconds (DMS) with microseconds to Decimal Degrees.
    """
    # Convert microseconds to milliseconds first
    milliseconds = microseconds / 1000  # Convert microseconds to milliseconds

    # Now convert DMS to decimal degrees
    decimal_value = degrees + (minutes / 60) + ((seconds + (milliseconds / 1000)) / 3600)

    print(f"DMS to Decimal: {degrees}° {minutes}' {seconds}\" {microseconds}µs → {decimal_value}")  # Debug print
    return decimal_value



def get_world_coordinates(model):
    """Extracts WGS84 latitude and longitude from IfcSite."""
    for site in model.by_type("IfcSite"):
        if hasattr(site, "RefLatitude") and hasattr(site, "RefLongitude"):
            #print(f"Extracted DMS from IFC: {site.RefLatitude}, {site.RefLongitude}")  # Debug
            return site.RefLatitude, site.RefLongitude
    return None  # No georeferencing found

def wgs84_to_lv95(latitude, longitude):
    """Convert WGS84 (EPSG:4326) to Swiss LV95 (EPSG:2056)."""
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:2056", always_xy=True)
    easting, northing = transformer.transform(longitude, latitude)
    print(f"WGS84 to LV95: {latitude}, {longitude} → {easting}, {northing}")  # Debug
    return easting, northing

"""Here is for comparing"""

def check_Max_Building_Height():
    return None

""" Here is for getting IFC floor coordinates"""

def get_floor_vertices(model):
    """
    Extract vertices of floor areas (e.g., from IfcSlab or IfcSpace).
    """
    settings = ifcopenshell.geom.settings()
    vertices = []

    # Loop through all IfcSlab objects to get the geometry
    for slab in model.by_type("IfcSlab"):
        shape = ifcopenshell.geom.create_shape(settings, slab)
        if shape:
            geometry = shape.geometry.verts  # Extract vertex list
            for i in range(0, len(geometry), 3):  # Extract (x, y, z) triplets
                vertices.append((geometry[i], geometry[i + 1], geometry[i + 2]))

    # Loop through IfcSpace objects (if no slabs found)
    if not vertices:
        for space in model.by_type("IfcSpace"):
            shape = ifcopenshell.geom.create_shape(settings, space)
            if shape:
                geometry = shape.geometry.verts  # Extract vertex list
                for i in range(0, len(geometry), 3):  # Extract (x, y, z) triplets
                    vertices.append((geometry[i], geometry[i + 1], geometry[i + 2]))

    return vertices


def get_boundary_polygon(vertices):
    """
    Project the 3D points onto a 2D plane (X, Y) and compute the boundary polygon.
    """
    # Project to 2D by removing Z values
    points_2d = [(x, y) for x, y, _ in vertices]

    # Compute the convex hull (outer boundary)
    multipoint = MultiPoint(points_2d)
    boundary_polygon = multipoint.convex_hull  # This creates a polygon

    return boundary_polygon

def get_boundary_polygon_lv95(vertices, lv95_coords):
    """
    Move the vertices by adding LV95 coordinates to their x and y values.
    Project the 3D points onto a 2D plane (X, Y) and compute the boundary polygon.
    
    :param vertices: List of tuples with (x, y, z) coordinates.
    :param lv95_coords: Tuple of LV95 coordinates (lv95_x, lv95_y).
    :return: A polygon representing the boundary.
    """
    # Shift the x and y values of each vertex by the LV95 coordinates
    shifted_points = [(x + lv95_coords[0], y + lv95_coords[1]) for x, y, _ in vertices]
    
    # Compute the convex hull (outer boundary)
    multipoint = MultiPoint(shifted_points)
    boundary_polygon = multipoint.convex_hull  # This creates the polygon

    return boundary_polygon


def get_centroid(polygon):
    """
    Compute the centroid of the polygon.
    """
    return polygon.centroid  # Shapely provides centroid directly

def move_vertices(vertices, x_translation, y_translation):
    """
    Move each vertex by a given translation vector (x_translation, y_translation).

    :param vertices: List of vertices in the format [(x, y, z), (x2, y2, z2), ...]
    :param x_translation: The amount to move in the x-direction.
    :param y_translation: The amount to move in the y-direction.
    :return: List of translated vertices.
    """
    # Apply the translation vector to each vertex
    translated_vertices = [(x + x_translation, y + y_translation, z) for x, y, z in vertices]
    
    return translated_vertices

"""Here is for testing"""


""" Implementing functions to simplify """

def get_lv95_coords(path):
    """
    """
    model = ifcopenshell.open(path)
    world_coords = get_world_coordinates(model)
    if world_coords:
        lat_dms, lon_dms = world_coords  # Unpack tuple of tuples

        # Convert DMS to Decimal
        latitude = dms_to_decimal(*lat_dms)  
        longitude = dms_to_decimal(*lon_dms)

        print(f"Converted WGS84: Lat = {latitude}, Lon = {longitude}")

        # Convert to Swiss LV95
        lv95_coords = wgs84_to_lv95(latitude, longitude)
        #print(f"Swiss LV95 Coordinates: {lv95_coords}")  # (Easting, Northing)
    else:
        print("No georeferencing data found in the IFC model.")
        # Get floor vertices

    return lv95_coords

def get_Building_data(path, zoning_map_path):
    """
    Function to extract building data and zoning information, including restrictions.
    """

    model = ifcopenshell.open(path)
    zoning_map = gpd.read_file(zoning_map_path)

    lv95_coords = get_lv95_coords(path)
    vertices = get_floor_vertices(model)
    boundary_polygon = get_boundary_polygon(vertices)
    vertices_95 = move_vertices(vertices, lv95_coords[0], lv95_coords[1])
    boundary_polygon_lv95 = get_boundary_polygon(vertices_95)
    centroid = get_centroid(boundary_polygon)

    # We are going to use the centroid of the polygonal projection of the slabs to get the position of the building.
    centroid_lv95 = [centroid.x + lv95_coords[0], centroid.y + lv95_coords[1]]

    # Generate building coordinates
    center_x, center_y = centroid_lv95[0], centroid_lv95[1]  
    building_coords = generate_building_coords(center_x, center_y)
    building_polygon = Polygon(building_coords)

    # Check zoning
    zoning_with_building = check_zoning(building_polygon, zoning_map)

    # Initialize building data dictionary
    building_data = {
        "height": get_building_height_complex(model)[0],
        "floors": len(model.by_type("IfcBuildingStorey")),
        "floor_area": get_floor_area(model)
    }

    # Store zoning information if available (only first matching zone)
    if not zoning_with_building.empty:
        first_zone = zoning_with_building.iloc[0]  # Take the first matching zoning entry
        building_data.update({
            "OBJID": first_zone['OBJID'],
            "R1_CODE": first_zone['R1_CODE'],
            "R1_BEZEICH": first_zone['R1_BEZEICH'],
            "R1_TYP_KAN": first_zone['R1_TYP_KAN']
        })

    # Extract zoning restrictions and add them to building_data
    zoning_restrictions = extract_zoning_restrictions(zoning_with_building)
    building_data.update(zoning_restrictions)  # Merging dictionaries

    # Save data to JSON
    save_to_json(building_data, "building_data.json")
    print("JSON saved successfully:", building_data)

    # Visualization
    plot_building_and_zoning(boundary_polygon_lv95, zoning_map, lv95_coords, vertices_95, zoom_factor=2)

    return building_data


# IFC Test File
path = r"tests\Test_Esri.ifc"

# SHP Test File
zoning_map_path = r"data\Zonenplan.shp"


# actual test
if __name__ == "__main__":
    if os.path.exists(path):
        if os.path.exists(zoning_map_path):
            get_Building_data(path, zoning_map_path)

        else:
            print("SHP file missing")

    else:
        print("IFC file missing")    

UPLOAD_DIR = "uploads"

@app.get("/")
def read_root():
    return {"message": "FastAPI is running!"}

@app.post("/upload/")
async def upload_ifc(file: UploadFile = File(...)):
    os.makedirs(UPLOAD_DIR, exist_ok=True)

    # Save the uploaded IFC file
    ifc_path = os.path.join(UPLOAD_DIR, file.filename)
    with open(ifc_path, "wb") as buffer:
        shutil.copyfileobj(file.file, buffer)

    # Load the existing SHP file
    if os.path.exists(zoning_map_path):
        try:
            gdf = gpd.read_file(zoning_map_path)
            shp_info = gdf.head().to_json()  # Convert first few rows to JSON
        except Exception as e:
            return {"error": f"Failed to read SHP file: {str(e)}"}
    else:
        return {"error": "SHP file is missing on the server"}

    # Call the function to process the IFC file after uploading it
    try:
        building_data = get_Building_data(ifc_path, zoning_map_path)  # Make sure the correct path is passed
    except Exception as e:
        return {"error": f"Failed to process IFC file: {str(e)}"}

    return {
        "message": "IFC file uploaded and processed successfully",
        "data": {
            "data": building_data

        }
    }


    
    