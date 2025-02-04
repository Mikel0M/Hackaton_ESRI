import ifcopenshell
import json
import os

"""To get simple height from OverallHeight"""
def get_building_height(model):
    heights = []
    for building in model.by_type("IfcBuilding"):
        if hasattr(building, "OverallHeight"):
            heights.append(building.OverallHeight)
    return max(heights) if heights else None

"""To get height from either OverallHeight or estimate it from the storey elevations"""
def get_building_height_complex(model):
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

"""Use get_floor_area_first_test if your file structure is simple and uses properties directly on elements 
(like GrossFloorArea or NetArea on storeys, slabs, or spaces)"""
def get_floor_area_simple(model):
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

#getting the data gathered
def get_building_data(model):
    
    building_data = {}

    # get building height (if available)
    building_data['height'] = get_building_height_complex(model)[0]

    # get number of floors
    building_data['floors'] = len(model.by_type("IfcBuildingStorey"))

    # get total floor area
    building_data['floor_area'] = get_floor_area(model)

    return building_data


""" to extract the data to a file given a path """
def save_to_json(data, file_path):
    with open(file_path, "w") as json_file:
        json.dump(data, json_file, indent=4) # pretty print with indent


""" Use get_floor_area for more complex files where area might be stored in related properties or element quantities, 
and if you want more flexibility in handling different ways area data can be structured"""
def get_floor_area(model):
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


path = r"tests\AC20-FZK-Haus.ifc"

if __name__ == "__main__":
    if os.path.exists(path):
        model = ifcopenshell.open(path)


        height, num_floors = get_building_height_complex(model)
        floor_area = get_floor_area(model)
        print(f"Building Height: {height} meters")
        print(f"Number of Floors: {num_floors}")
        print(f"Total Floor Area: {floor_area} square meters")
        print("Checking available elements in IFC file:")
        print(f"Number of IfcBuildingStorey: {len(model.by_type('IfcBuildingStorey'))}")
        print(f"Number of IfcSlab: {len(model.by_type('IfcSlab'))}")
        print(f"Number of IfcSpace: {len(model.by_type('IfcSpace'))}")
        building_data = get_building_data(model)
        save_to_json(building_data, "building_data.json")
    else:
        print("Error: IFC file not found!")

    
    