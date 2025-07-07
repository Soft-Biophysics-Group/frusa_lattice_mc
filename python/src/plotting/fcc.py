"""
Vincent Ouazan-Reboul, 2025
Functions to plot the results of simulations of FCC particles.
So far only supports one type of particles.
"""

from pathlib import Path

path_to_config = Path(__file__).parent.parent

ALL_CONTACTS = [(i, j) for i in range(24) for j in range(24)]
# DEFAULT_MATERIAL = (14, 0, 255, 0.1)
DEFAULT_MATERIAL = (14, 0, 255, 1.0)

src_plotting_dir = Path(__file__).parent

path_to_numbered_rhombic = (
    src_plotting_dir / "assets/oneRhombicDodecahedron/one_rhombic_numbered.obj"
)
path_to_one_axis_rhombic = (
    src_plotting_dir / "assets/oneRhombicDodecahedron/one_rhombic_one_axis.obj"
)
path_to_one_axis_colored_rhombic = (
    src_plotting_dir / "assets/oneRhombicDodecahedron/one_rhombic_one_axis_colored.obj"
)
path_to_one_axis_diagonal_rhombic = (
    src_plotting_dir / "assets/oneRhombicDodecahedron/one_rhombic_diagonal.obj"
)

path_to_one_axis_diagonal_colored_rhombic = (
    src_plotting_dir / "assets/oneRhombicDodecahedron/one_rhombic_diagonal_colored.obj"
)
path_to_one_axis_diagonal_colored_gray = (
    src_plotting_dir / "assets/oneRhombicDodecahedron/one_rhombic_diagonal_gray.obj"
)
paths = {
    "path_to_numbered_rhombic": path_to_numbered_rhombic,
    "path_to_one_axis_rhombic": path_to_one_axis_rhombic,
    "path_to_one_axis_diagonal_rhombic": path_to_one_axis_diagonal_rhombic,
    "path_to_one_axis_diagonal_colored_rhombic": path_to_one_axis_diagonal_colored_rhombic,
    "path_to_one_axis_diagonal_colored_gray": path_to_one_axis_diagonal_colored_gray,
}


def plot_boundary(
    site_1,
    site_2,
    orientation_1,
    orientation_2,
    lattice,
    boundary_contacts,
    material,
    collection_name="Rhombii",
):
    site_1_lattice_coords = lattice.lattice_site_to_lattice_coords(site_1)
    site_2_lattice_coords = lattice.lattice_site_to_lattice_coords(site_2)
    bond_lattice_coords = site_2_lattice_coords - site_1_lattice_coords
    # Case of periodic boundary conditions: put norm of bond back to 1
    # And invert its direction
    if np.abs(bond_lattice_coords.sum()) != 1.0:
        bond_lattice_coords //= -bond_lattice_coords.sum()
    bond_index = lattice.get_bond(bond_lattice_coords)
    face_1, face_2 = CubicParticle().get_faces_in_contact(
        orientation_1, orientation_2, bond_index
    )
    equiv_contacts = CubicParticle().get_equivalent_face_pairs(face_1, face_2)

    for equiv_face_1, equiv_face_2 in equiv_contacts:
        if (equiv_face_1, equiv_face_2) in boundary_contacts:
            site_1_coords = 2 * lattice.lattice_site_to_lattice_coords(site_1)
            # Factor of 2 is because of rescaling in blender
            face_cart_coords = site_1_coords + bond_lattice_coords
            plot_square(
                face_cart_coords,
                bond_lattice_coords,
                lattice,
                material,
                collection_name=collection_name,
            )

    return

def plot_lozenge(
    face_cart_coords,
    bond_lattice_coords,
    lattice,
    material = None,
    size=2.0,
    collection_name="Boundaries",
):
    collection = get_or_create_collection(collection_name)

    half_size = size / 2
    y_mid_vertices = 0
    # This is the center of the face.
    # ALong which direction of space is face oriented?
    face_dir = np.where(np.array(bond_lattice_coords) != 0)[0][0]
    other_dirs = [0, 1, 2]
    other_dirs.pop(face_dir)

    vertices = np.zeros((4, 3))
    vertices[0, other_dirs[0]] += half_size
    vertices[0, other_dirs[1]] += half_size
    vertices[1, other_dirs[0]] += half_size
    vertices[1, other_dirs[1]] -= half_size
    vertices[2, other_dirs[0]] -= half_size
    vertices[2, other_dirs[1]] -= half_size
    vertices[3, other_dirs[0]] -= half_size
    vertices[3, other_dirs[1]] += half_size

    # Define faces (Blender uses quads)
    faces = [(0, 1, 2, 3)]
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]

    # Create a new mesh and object
    mesh = bpy.data.meshes.new(name="SquareMesh")
    obj = bpy.data.objects.new(name="Square", object_data=mesh)
    collection.objects.link(obj)
    bpy.context.collection.objects.link(obj)
    mesh.from_pydata(vertices, edges, faces)
    mesh.update()

    # Create an edge-only object
    edge_mesh = bpy.data.meshes.new(name="EdgeMesh")
    edge_obj = bpy.data.objects.new(name="Edges", object_data=edge_mesh)
    bpy.context.collection.objects.link(edge_obj)
    # Create edges as a separate mesh
    edge_mesh.from_pydata(vertices, edges, [])
    edge_mesh.update()

    # Move the object to the given location
    obj.location = face_cart_coords
    obj.data.materials.append(material)

    # Create edge material and assign it
    edge_color = (0.0, 0.0, 0.0, 1.0)
    edge_mat = create_material(name="EdgeMaterial", color=edge_color)
    edge_obj.data.materials.append(edge_mat)
    edge_obj.location = face_cart_coords

    return obj, edge_obj

