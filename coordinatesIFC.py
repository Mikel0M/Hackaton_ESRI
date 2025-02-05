import ifcopenshell
from shapely.geometry import Polygon
import os

import ifcopenshell
from shapely.geometry import Polygon

def get_building_polygon(model):
    """
    Extracts building polygon coordinates from IfcSlab, IfcSpace, or IfcPolygonalFaceSet objects in an IFC 4 file.
    Returns the coordinates as a list of tuples (x, y).
    """
    building_coords = []

    # Get all the vertices (IfcVertex) in the model first, indexed by their ID
    vertices = {}
    for vertex in model.by_type("IfcVertex"):
        if hasattr(vertex, "Coordinates"):
            coordinates = vertex.Coordinates
            print(f"Vertex ID: {vertex.id()}, Coordinates: {coordinates}")
            vertices[vertex.id()] = (coordinates[0], coordinates[1])  # Store only X, Y coordinates
        else:
            print(f"Warning: Vertex ID {vertex.id()} does not have coordinates.")

    # Log the number of vertices found
    print(f"Found {len(vertices)} vertices.")

    # Check for alternative types of geometry
    for slab in model.by_type("IfcSlab"):
        for geometry in slab.Representation.Representations:
            if geometry.is_a("IfcShapeRepresentation"):
                for item in geometry.Items:
                    if item.is_a("IfcPolygonalFaceSet"):
                        if hasattr(item, "Faces"):
                            faces = item.Faces
                            for face in faces:
                                if face.is_a("IfcIndexedPolygonalFace"):
                                    indices = face.CoordIndex
                                    print(f"Processing IfcIndexedPolygonalFace with indices: {indices}")
                                    coords = []
                                    for index in indices:
                                        print(f"Checking vertex index: {index}")
                                        if index > 0:  # Ensure the index is greater than 0
                                            index -= 1  # Convert to 0-based index
                                        if index in vertices:
                                            coords.append(vertices[index])
                                        else:
                                            print(f"Warning: Vertex index {index+1} not found in vertices.")
                                    if coords:
                                        building_coords.append(coords)
                                    else:
                                        print("Warning: No valid coordinates found for this face.")
                                elif face.is_a("IfcIndexedPolygonalFaceWithVoids"):
                                    indices = face.CoordIndex
                                    coords = []
                                    print(f"Processing IfcIndexedPolygonalFaceWithVoids with indices: {indices}")
                                    for index in indices:
                                        if index > 0:  # Ensure the index is greater than 0
                                            index -= 1  # Convert to 0-based index
                                        if index in vertices:
                                            coords.append(vertices[index])
                                        else:
                                            print(f"Warning: Vertex index {index+1} not found in vertices.")
                                    if coords:
                                        building_coords.append(coords)
                                    else:
                                        print("Warning: No valid coordinates found for this face.")

                    # Check if the item is an IfcFacetedBrep (which can contain geometry as well)
                    elif item.is_a("IfcFacetedBrep"):
                        print(f"Processing IfcFacetedBrep...")
                        for face in item.Faces:
                            if hasattr(face, "Bounds"):
                                # Extract the coordinates of the face's boundary
                                for bound in face.Bounds:
                                    if hasattr(bound, "Vertices"):
                                        for vertex in bound.Vertices:
                                            if hasattr(vertex, "Coordinates"):
                                                coords = (vertex.Coordinates[0], vertex.Coordinates[1])
                                                print(f"Found coordinates in IfcFacetedBrep: {coords}")
                                                building_coords.append(coords)

    # If no coordinates are found, return a fallback set
    if not building_coords:
        print("Warning: No polygonal face found for the building geometry. Using fallback coordinates.")
        building_coords = [(2696755, 1261971), (2696756, 1261971), (2696756, 1261972), (2696755, 1261972)]

    return building_coords


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
            print("IfcBuildingStorey")
    
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
                                print("IfcSlab")

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
                                print("IfcSpace")

    return total_area if total_area > 0 else "Unknown area"


# Testing with your IFC model
if __name__ == "__main__":
    # Path to your IFC model file
    ifc_file_path = "tests/Test_hackaton.ifc"

    # Open the IFC model
    if os.path.exists(ifc_file_path):
        model = ifcopenshell.open(ifc_file_path)

        # Call the function to get the building polygon
        building_coords = get_building_polygon(model)
        #building_area = get_floor_area_simple(model)

        # Output the coordinates
        print(f"Building Coordinates: {building_coords}")
        #print(f"Building Area: {building_area}")
    else:
        print("IFC file not found!")
