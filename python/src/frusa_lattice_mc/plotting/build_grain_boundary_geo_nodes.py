import bpy, mathutils

#initialize geometry_nodes node group
def geometry_nodes_node_group(material_name):
    geometry_nodes = bpy.data.node_groups.new(type = 'GeometryNodeTree', name = "Geometry Nodes")

    geometry_nodes.color_tag = 'NONE'
    geometry_nodes.description = ""
    geometry_nodes.default_group_node_width = 140
    

    geometry_nodes.is_modifier = True

    #geometry_nodes interface
    #Socket Geometry
    geometry_socket = geometry_nodes.interface.new_socket(name = "Geometry", in_out='OUTPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket.attribute_domain = 'POINT'

    #Socket Geometry
    geometry_socket_1 = geometry_nodes.interface.new_socket(name = "Geometry", in_out='INPUT', socket_type = 'NodeSocketGeometry')
    geometry_socket_1.attribute_domain = 'POINT'


    #initialize geometry_nodes nodes
    #node Group Input.001
    group_input_001 = geometry_nodes.nodes.new("NodeGroupInput")
    group_input_001.name = "Group Input.001"

    #node Group Output.001
    group_output_001 = geometry_nodes.nodes.new("NodeGroupOutput")
    group_output_001.name = "Group Output.001"
    group_output_001.is_active_output = True

    #node Mesh to Curve
    mesh_to_curve = geometry_nodes.nodes.new("GeometryNodeMeshToCurve")
    mesh_to_curve.name = "Mesh to Curve"
    #Selection
    mesh_to_curve.inputs[1].default_value = True

    #node Curve to Mesh
    curve_to_mesh = geometry_nodes.nodes.new("GeometryNodeCurveToMesh")
    curve_to_mesh.name = "Curve to Mesh"
    #Fill Caps
    curve_to_mesh.inputs[2].default_value = False

    #node Join Geometry
    join_geometry = geometry_nodes.nodes.new("GeometryNodeJoinGeometry")
    join_geometry.name = "Join Geometry"

    #node Set Material
    set_material = geometry_nodes.nodes.new("GeometryNodeSetMaterial")
    set_material.name = "Set Material"
    #Selection
    set_material.inputs[1].default_value = True
    set_material.inputs[2].default_value = bpy.data.materials[material_name]

    #node Quadrilateral
    quadrilateral = geometry_nodes.nodes.new("GeometryNodeCurvePrimitiveQuadrilateral")
    quadrilateral.name = "Quadrilateral"
    quadrilateral.mode = 'RECTANGLE'
    #Width
    quadrilateral.inputs[0].default_value = 0.10000000149011612
    #Height
    quadrilateral.inputs[1].default_value = 0.10000000149011612





    #Set locations
    group_input_001.location = (-693.641845703125, 90.26908874511719)
    group_output_001.location = (378.86260986328125, 70.57115173339844)
    mesh_to_curve.location = (-320.322509765625, 169.60801696777344)
    curve_to_mesh.location = (-177.20066833496094, 60.37898254394531)
    join_geometry.location = (230.75022888183594, 51.678619384765625)
    set_material.location = (2.3581085205078125, 23.295822143554688)
    quadrilateral.location = (-539.0096435546875, -105.19743347167969)

    #Set dimensions
    group_input_001.width, group_input_001.height = 140.0, 100.0
    group_output_001.width, group_output_001.height = 140.0, 100.0
    mesh_to_curve.width, mesh_to_curve.height = 140.0, 100.0
    curve_to_mesh.width, curve_to_mesh.height = 140.0, 100.0
    join_geometry.width, join_geometry.height = 140.0, 100.0
    set_material.width, set_material.height = 140.0, 100.0
    quadrilateral.width, quadrilateral.height = 140.0, 100.0

    #initialize geometry_nodes links
    #group_input_001.Geometry -> mesh_to_curve.Mesh
    geometry_nodes.links.new(group_input_001.outputs[0], mesh_to_curve.inputs[0])
    #mesh_to_curve.Curve -> curve_to_mesh.Curve
    geometry_nodes.links.new(mesh_to_curve.outputs[0], curve_to_mesh.inputs[0])
    #set_material.Geometry -> join_geometry.Geometry
    geometry_nodes.links.new(set_material.outputs[0], join_geometry.inputs[0])
    #join_geometry.Geometry -> group_output_001.Geometry
    geometry_nodes.links.new(join_geometry.outputs[0], group_output_001.inputs[0])
    #curve_to_mesh.Mesh -> set_material.Geometry
    geometry_nodes.links.new(curve_to_mesh.outputs[0], set_material.inputs[0])
    #quadrilateral.Curve -> curve_to_mesh.Profile Curve
    geometry_nodes.links.new(quadrilateral.outputs[0], curve_to_mesh.inputs[1])
    #group_input_001.Geometry -> join_geometry.Geometry
    geometry_nodes.links.new(group_input_001.outputs[0], join_geometry.inputs[0])
    return geometry_nodes
