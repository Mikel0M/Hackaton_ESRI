"""
Created on Tue Feb  4 10:12:31 2025

@author: mikel
"""

import ifcopenshell
import ifcopenshell.geom
import json
import os
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import MultiPoint, Polygon
from pyproj import Transformer


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



def get_building_data(model):
    """
    getting the data gathered
    """
    
    building_data = {}

    # get building height (if available)
    building_data['height'] = get_building_height_complex(model)[0]

    # get number of floors
    building_data['floors'] = len(model.by_type("IfcBuildingStorey"))

    # get total floor area
    building_data['floor_area'] = get_floor_area(model)



    return building_data


def save_to_json(data, file_path):
    """
    to extract the data to a file given a path
    """
    with open(file_path, "w") as json_file:
        json.dump(data, json_file, indent=4) # pretty print with indent



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
    Extract and print zoning restrictions for the building.
    """
    if not zoning_with_building.empty:
        max_height = zoning_with_building.get('GEBAEUDEHO', None).values[0]
        total_height = zoning_with_building.get('GESAMTHOEH', None).values[0]
        max_floors = zoning_with_building.get('VOLLGESCHO', None).values[0]
        attic_floors = zoning_with_building.get('DACHGESCHO', None).values[0]
        floor_area_ratio = zoning_with_building.get('AUSNUETZUN', None).values[0]

        print("🏗 **Zoning Restrictions for the Building:**")
        print(f"- Max Building Height: {max_height} meters")
        print(f"- Total Allowed Height: {total_height} meters")
        print(f"- Max Full Floors: {max_floors} floors")
        print(f"- Allowed Attic Floors: {attic_floors} floors")
        print(f"- Floor Area Ratio (FAR): {floor_area_ratio}")
    else:
        print("🚨 Building is not within any zoning area!")


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

# IFC Test File
path = r"tests\Hackaton_test.ifc"

# SHP Test File
zoning_map_path = r"data\Zonenplan.shp"
zoning_map = gpd.read_file(zoning_map_path)

# actual test
if __name__ == "__main__":
    if os.path.exists(path):
        model = ifcopenshell.open(path)

        if model:
            world_coords = get_world_coordinates(model)
            if world_coords:
                lat_dms, lon_dms = world_coords  # Unpack tuple of tuples

                # Debug Print DMS
                #print(f"Raw DMS Latitude: {lat_dms}")
                #print(f"Raw DMS Longitude: {lon_dms}")

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

        vertices = get_floor_vertices(model)

        if vertices:
       

            boundary_polygon = get_boundary_polygon(vertices)
            vertices_95 = move_vertices(vertices, lv95_coords[0], lv95_coords[1])
            boundary_polygon_lv95 = get_boundary_polygon(vertices_95)
            centroid = get_centroid(boundary_polygon)


        # We are going to use the centroid of the polygonal projection of the slabs to get the position of the building.
        # This way, if the 0 point is outside of the building it will not be a problem.
        centroid_lv95 = [centroid.x + lv95_coords[0],centroid.y + lv95_coords[1]]

        #now, in order to represent the building

        # this part reads IFC and makes a couple of tests
        height, num_floors = get_building_height_complex(model)
        floor_area = get_floor_area(model)
        building_data = get_building_data(model)
        save_to_json(building_data, "building_data.json")
        print("json saved")
        # this part is intended to test GIS interactions using Geopandas
        # Example usage
        center_x, center_y = centroid_lv95[0], centroid_lv95[1]  # Example center coordinates
        building_coords = generate_building_coords(center_x, center_y)
        building_polygon = Polygon(building_coords)

        # Assuming zoning_map is a GeoDataFrame containing zoning polygons
        zoning_with_building = check_zoning(building_polygon, zoning_map)

        if not zoning_with_building.empty:
            print("Building is within the following zoning areas:")
            for idx, row in zoning_with_building.iterrows():
                print(f"Zoning Area {idx + 1}:")
                print(f"  OBJID: {row['OBJID']}")
                print(f"  R1_CODE: {row['R1_CODE']}")
                print(f"  R1_BEZEICH: {row['R1_BEZEICH']}")
                print(f"  R1_TYP_KAN: {row['R1_TYP_KAN']}")
                #print(f"  Geometry: {row['geometry']}")
                print()

        plot_building_and_zoning(boundary_polygon_lv95, zoning_map, lv95_coords,vertices_95, zoom_factor = 2)
        extract_zoning_restrictions(zoning_with_building)

    else:
        print("Error: IFC file not found!")





    
    