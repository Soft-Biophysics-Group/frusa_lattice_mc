import bpy, mathutils

def build_composition_nodes(
    min_map_range_distance: float = 60.0,
    max_map_range_distance: float = 75.0,
    color_ramp_max_value: float = 65.0,
    erode_filter_distance:int = 2
):
# Generate unique scene name
    base_name = "Scene"
    end_name = base_name
    if bpy.data.scenes.get(end_name) != None:
        i = 1
        end_name = base_name + f".{i:03d}"
        while bpy.data.scenes.get(end_name) != None:
            end_name = base_name + f".{i:03d}"
            i += 1

    scene = bpy.context.window.scene

    scene.name = end_name
    scene.use_fake_user = True
    bpy.context.window.scene = scene
    scene.use_nodes = True
#initialize Scene node group
    scene_1 = scene.node_tree
    #start with a clean node tree
    for node in scene_1.nodes:
        scene_1.nodes.remove(node)
    scene_1.color_tag = "NONE"
    scene_1.description = ""
    scene_1.default_group_node_width = 140

    #scene_1 interface

    #initialize scene_1 nodes
    #node Render Layers
    render_layers = scene_1.nodes.new("CompositorNodeRLayers")
    render_layers.name = "Render Layers"
    render_layers.layer = 'WithoutBox'

    #node Composite
    composite = scene_1.nodes.new("CompositorNodeComposite")
    composite.name = "Composite"
    composite.use_alpha = True
    #Alpha
    composite.inputs[1].default_value = 1.0

    #node Mix
    mix = scene_1.nodes.new("CompositorNodeMixRGB")
    mix.name = "Mix"
    mix.blend_type = 'MULTIPLY'
    mix.use_alpha = False
    mix.use_clamp = False
    #Fac
    mix.inputs[0].default_value = 1.0

    #node Map Range
    map_range = scene_1.nodes.new("CompositorNodeMapRange")
    map_range.name = "Map Range"
    map_range.use_clamp = False
    #From Min
    map_range.inputs[1].default_value = min_map_range_distance
    #From Max
    map_range.inputs[2].default_value = max_map_range_distance
    #To Min
    map_range.inputs[3].default_value = 0.0
    #To Max
    map_range.inputs[4].default_value = 1.0

    #node Filter
    filter = scene_1.nodes.new("CompositorNodeFilter")
    filter.name = "Filter"
    filter.filter_type = 'SOBEL'
    #Fac
    filter.inputs[0].default_value = 1.0

    #node Color Ramp
    color_ramp = scene_1.nodes.new("CompositorNodeValToRGB")
    color_ramp.name = "Color Ramp"
    color_ramp.color_ramp.color_mode = 'RGB'
    color_ramp.color_ramp.hue_interpolation = 'NEAR'
    color_ramp.color_ramp.interpolation = 'CONSTANT'

    #initialize color ramp elements
    color_ramp.color_ramp.elements.remove(color_ramp.color_ramp.elements[0])
    color_ramp_cre_0 = color_ramp.color_ramp.elements[0]
    color_ramp_cre_0.position = 0.0
    color_ramp_cre_0.alpha = 0.0
    color_ramp_cre_0.color = (0.0, 0.0, 0.0, 0.0)

    color_ramp_cre_1 = color_ramp.color_ramp.elements.new(color_ramp_max_value)
    color_ramp_cre_1.alpha = 1.0
    color_ramp_cre_1.color = (1.0, 1.0, 1.0, 1.0)


    #node Dilate/Erode
    dilate_erode = scene_1.nodes.new("CompositorNodeDilateErode")
    dilate_erode.name = "Dilate/Erode"
    dilate_erode.distance = erode_filter_distance
    dilate_erode.edge = 0.0
    dilate_erode.falloff = 'SMOOTH'
    dilate_erode.mode = 'DISTANCE'

    #node Viewer.003
    viewer_003 = scene_1.nodes.new("CompositorNodeViewer")
    viewer_003.name = "Viewer.003"
    viewer_003.use_alpha = True
    #Alpha
    viewer_003.inputs[1].default_value = 1.0

    #node Render Layers.001
    render_layers_001 = scene_1.nodes.new("CompositorNodeRLayers")
    render_layers_001.name = "Render Layers.001"
    render_layers_001.layer = 'WithBox'

    #node Mix.001
    mix_001 = scene_1.nodes.new("CompositorNodeMixRGB")
    mix_001.name = "Mix.001"
    mix_001.blend_type = 'OVERLAY'
    mix_001.use_alpha = False
    mix_001.use_clamp = False
    #Image_001
    mix_001.inputs[2].default_value = (0.745, 0.122, 0.068, 1.0)

    #node Viewer
    viewer = scene_1.nodes.new("CompositorNodeViewer")
    viewer.name = "Viewer"
    viewer.use_alpha = True
    #Alpha
    viewer.inputs[1].default_value = 1.0


    #Set locations
    render_layers.location = (-386.3139953613281, 302.2930908203125)
    composite.location = (1750.862548828125, 275.450439453125)
    mix.location = (1040.0540771484375, 245.51763916015625)
    map_range.location = (15.223779678344727, 492.8875732421875)
    filter.location = (233.68077087402344, 342.34344482421875)
    color_ramp.location = (690.3121337890625, 402.66375732421875)
    dilate_erode.location = (464.23907470703125, 460.73236083984375)
    viewer_003.location = (1088.6759033203125, 554.6784057617188)
    render_layers_001.location = (-384.95684814453125, -205.5657958984375)
    mix_001.location = (1269.6842041015625, 57.699554443359375)
    viewer.location = (1578.9189453125, 345.1446838378906)

    #Set dimensions
    render_layers.width, render_layers.height = 240.0, 100.0
    composite.width, composite.height = 100.908447265625, 100.0
    mix.width, mix.height = 144.8297119140625, 100.0
    map_range.width, map_range.height = 140.0, 100.0
    filter.width, filter.height = 140.0, 100.0
    color_ramp.width, color_ramp.height = 240.0, 100.0
    dilate_erode.width, dilate_erode.height = 140.0, 100.0
    viewer_003.width, viewer_003.height = 140.0, 100.0
    render_layers_001.width, render_layers_001.height = 240.0, 100.0
    mix_001.width, mix_001.height = 140.0, 100.0
    viewer.width, viewer.height = 140.0, 100.0

    #initialize scene_1 links
    #render_layers.Depth -> map_range.Value
    scene_1.links.new(render_layers.outputs[2], map_range.inputs[0])
    #filter.Image -> dilate_erode.Mask
    scene_1.links.new(filter.outputs[0], dilate_erode.inputs[0])
    #map_range.Value -> filter.Image
    scene_1.links.new(map_range.outputs[0], filter.inputs[1])
    #color_ramp.Image -> viewer_003.Image
    scene_1.links.new(color_ramp.outputs[0], viewer_003.inputs[0])
    #dilate_erode.Mask -> color_ramp.Fac
    scene_1.links.new(dilate_erode.outputs[0], color_ramp.inputs[0])
    #render_layers.Alpha -> mix.Image
    scene_1.links.new(render_layers.outputs[1], mix.inputs[2])
    #color_ramp.Image -> mix.Image
    scene_1.links.new(color_ramp.outputs[0], mix.inputs[1])
    #mix_001.Image -> composite.Image
    scene_1.links.new(mix_001.outputs[0], composite.inputs[0])
    #mix_001.Image -> viewer.Image
    scene_1.links.new(mix_001.outputs[0], viewer.inputs[0])
    #mix.Image -> mix_001.Fac
    scene_1.links.new(mix.outputs[0], mix_001.inputs[0])
    #render_layers_001.Image -> mix_001.Image
    scene_1.links.new(render_layers_001.outputs[0], mix_001.inputs[1])
    return scene_1
