import ifcopenshell
import json
import os
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon

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

def plot_building_and_zoning(building_polygon, zoning_map, buffer_distance=100):
    """
    Plot the zoning map, the building, and its buffer.
    """
    building_buffer = building_polygon.buffer(buffer_distance)

    fig, ax = plt.subplots(figsize=(10, 10))
    zoning_map.plot(ax=ax, color='lightgray', edgecolor='black')

    minx, miny, maxx, maxy = zoning_map.total_bounds
    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)

    # Plot building
    x, y = building_polygon.exterior.xy
    ax.fill(x, y, color='red', alpha=0.5)

    # Plot buffer
    x_buffer, y_buffer = building_buffer.exterior.xy
    ax.fill(x_buffer, y_buffer, color='blue', alpha=0.3)

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


"""Here is for testing"""

# IFC Test File
path = r"tests\AC20-FZK-Haus.ifc"

# SHP Test File
zoning_map_path = r"data\Zonenplan.shp"
zoning_map = gpd.read_file(zoning_map_path)

# actual test
if __name__ == "__main__":
    if os.path.exists(path):
        model = ifcopenshell.open(path)

        # this part reads IFC and makes a couple of tests
        """height, num_floors = get_building_height_complex(model)
        floor_area = get_floor_area(model)
        building_data = get_building_data(model)
        save_to_json(building_data, "building_data.json")
        print("json saved")
    else:
        print("Error: IFC file not found!")"""

        # this part is intended to test GIS interactions using Geopandas
        # Example usage
        center_x, center_y = 2696411, 1262181  # Example center coordinates
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
                print(f"  Geometry: {row['geometry']}")
                print()

        plot_building_and_zoning(building_polygon, zoning_map)
        extract_zoning_restrictions(zoning_with_building)



    
    