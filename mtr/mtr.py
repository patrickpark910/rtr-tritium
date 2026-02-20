import openmc
import numpy as np

# ==============================================================================
# 1. MATERIALS DEFINITION
# ==============================================================================

# U3Si2 Fuel - Low Enriched Uranium (19.75% U-235)
# Density is approximately 12.2 g/cm3.
# Weight fractions: Uranium is ~92.7%, Silicon is ~7.3%.
fuel = openmc.Material(name='LEU U3Si2')
fuel.set_density('g/cm3', 12.2)
fuel.add_nuclide('U235', 0.927 * 0.1975, 'wo')
fuel.add_nuclide('U238', 0.927 * 0.8025, 'wo')
fuel.add_element('Si', 0.073, 'wo')

# Aluminum-6061 for Cladding and Structure
clad = openmc.Material(name='Aluminum Clad')
clad.set_density('g/cm3', 2.7)
clad.add_element('Al', 1.0)

# Light Water Coolant/Moderator/Reflector
water = openmc.Material(name='Light Water')
water.set_density('g/cm3', 0.998)
water.add_element('H', 2.0)
water.add_element('O', 1.0)
water.add_s_alpha_beta('c_H_in_H2O') # S(alpha, beta) for thermal scattering

materials = openmc.Materials([fuel, clad, water])
materials.export_to_xml()

# ==============================================================================
# 2. GEOMETRY DEFINITION
# ==============================================================================

# --- Dimensions (in cm) ---
# Plate dimensions
meat_thickness = 0.051
clad_thickness = 0.127
channel_thickness = 0.254
plate_pitch = clad_thickness + channel_thickness # 0.381 cm

# Fuel Assembly dimensions
num_plates = 21
assembly_width = num_plates * plate_pitch # X-dimension (8.001 cm)
assembly_depth = 7.6                      # Y-dimension (7.6 cm)
active_height = 60.0                      # Z-dimension (60.0 cm)

# --- Surfaces for a Single Fuel Plate Unit Cell ---
meat_x = openmc.XPlane(x0=-meat_thickness/2, name='Meat Left')
meat_x.x0 = -meat_thickness / 2
meat_x_right = openmc.XPlane(x0=meat_thickness/2, name='Meat Right')

clad_x_left = openmc.XPlane(x0=-clad_thickness/2, name='Clad Left')
clad_x_right = openmc.XPlane(x0=clad_thickness/2, name='Clad Right')

pitch_x_left = openmc.XPlane(x0=-plate_pitch/2, name='Pitch Left')
pitch_x_right = openmc.XPlane(x0=plate_pitch/2, name='Pitch Right')

# --- Cells for the Unit Plate ---
# Meat
meat_cell = openmc.Cell(fill=fuel, region=+meat_x & -meat_x_right)
# Cladding
clad_cell = openmc.Cell(fill=clad, region=(+clad_x_left & -meat_x) | (+meat_x_right & -clad_x_right))
# Coolant Channel
coolant_cell = openmc.Cell(fill=water, region=(+pitch_x_left & -clad_x_left) | (+clad_x_right & -pitch_x_right))

# Universe representing a single infinite plate along Y and Z
plate_univ = openmc.Universe(cells=[meat_cell, clad_cell, coolant_cell])

# --- Construct the MTR Assembly ---
# We use a 1D RectLattice in the X direction to stack the plates
assembly_lattice = openmc.RectLattice(name='MTR Assembly')
assembly_lattice.lower_left = [-assembly_width/2, -assembly_depth/2]
assembly_lattice.pitch = [plate_pitch, assembly_depth]
# Create a 1x21 grid
assembly_lattice.universes = [[plate_univ] * num_plates]

# Boundary for the assembly (Active Region)
assembly_box = openmc.model.RectangularPrism(width=assembly_width, height=assembly_depth)
# Create the top and bottom bounding planes
min_z = openmc.ZPlane(z0=-active_height/2)
max_z = openmc.ZPlane(z0=active_height/2)

# The active region is above min_z AND below max_z
active_z = +min_z & -max_z


assembly_cell = openmc.Cell(fill=assembly_lattice, region=-assembly_box & active_z)

# Fill the top/bottom empty regions of the assembly with water
assembly_water = openmc.Cell(fill=water, region=~(-assembly_box & active_z))
assembly_univ = openmc.Universe(cells=[assembly_cell, assembly_water])

# --- Construct the Core Lattice (3x3 Box) ---
core_lattice = openmc.RectLattice(name='3x3 Core Lattice')
core_lattice.lower_left = [-assembly_width * 1.5, -assembly_depth * 1.5]
core_lattice.pitch = [assembly_width, assembly_depth]

# A 3x3 grid of our fuel assemblies
core_lattice.universes = [
    [assembly_univ, assembly_univ, assembly_univ],
    [assembly_univ, assembly_univ, assembly_univ],
    [assembly_univ, assembly_univ, assembly_univ]
]

# --- Global Geometry / Reflector ---
# We'll put the core inside a 1 meter radius tank of water
# --- Global Geometry / Reflector ---
# We'll put the core inside a 1 meter radius tank of water, 2 meters tall
tank_radius = openmc.ZCylinder(r=100.0, boundary_type='vacuum')
tank_min_z = openmc.ZPlane(z0=-100.0, boundary_type='vacuum')
tank_max_z = openmc.ZPlane(z0=100.0, boundary_type='vacuum')

# Create the core cell, bounded tightly around the 3x3 lattice
core_box = openmc.model.RectangularPrism(width=assembly_width*3, height=assembly_depth*3)
core_cell = openmc.Cell(fill=core_lattice, region=-core_box & active_z)

# The reflector fills everything else inside the tank cylinder AND between the tank Z-planes
reflector_region = -tank_radius & +tank_min_z & -tank_max_z & ~(-core_box & active_z)
reflector_cell = openmc.Cell(fill=water, region=reflector_region)

# Create the core cell, bounded tightly around the 3x3 lattice
core_box = openmc.model.RectangularPrism(width=assembly_width*3, height=assembly_depth*3)
core_cell = openmc.Cell(fill=core_lattice, region=-core_box & active_z)

# The reflector fills everything else inside the tank
reflector_cell = openmc.Cell(fill=water, region=~(-core_box & active_z) & -tank_radius)

root_universe = openmc.Universe(cells=[core_cell, reflector_cell])
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# ==============================================================================
# 3. SETTINGS & EXECUTION
# ==============================================================================

settings = openmc.Settings()
settings.batches = 100
settings.inactive = 20
settings.particles = 10000

# Set a bounding box spatial distribution for the initial source
lower_left = [-assembly_width*1.5, -assembly_depth*1.5, -active_height/2]
upper_right = [assembly_width*1.5, assembly_depth*1.5, active_height/2]
uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=True)
settings.source = openmc.IndependentSource(space=uniform_dist)

settings.export_to_xml()

# ==============================================================================
# 4. PLOTTING
# ==============================================================================

# Define a common color scheme dictionary for the MTR materials
material_colors = {
    fuel: 'orange',
    clad: 'gray',
    water: 'cyan'
}

# --- Plot 1: XY Plane (Top-down view) ---
plot_xy = openmc.Plot()
plot_xy.filename = 'mtr_xy_plot'
plot_xy.basis = 'xy'
# Center exactly on the core midplane
plot_xy.origin = (0, 0, 0) 
# Zoomed to 40x40 cm to nicely frame the ~24x23 cm 3x3 core
plot_xy.width = (40, 40)      
plot_xy.pixels = (2000, 2000)   # High resolution needed to resolve thin plates
plot_xy.color_by = 'material'
plot_xy.colors = material_colors

# --- Plot 2: XZ Plane (Side view) ---
plot_xz = openmc.Plot()
plot_xz.filename = 'mtr_xz_plot'
plot_xz.basis = 'xz'
# Center exactly on the core midplane
plot_xz.origin = (0, 0, 0)
# 40 cm wide to show the core, 80 cm tall to capture the 60 cm active height + axial reflectors
plot_xz.width = (40, 80)      
plot_xz.pixels = (1000, 2000)   # Maintain aspect ratio
plot_xz.color_by = 'material'
plot_xz.colors = material_colors

# Add both plots to a Plots collection and export
plots = openmc.Plots([plot_xy, plot_xz])
plots.export_to_xml()

# Run the plot generation
openmc.plot_geometry()


print("Model successfully created. Run 'openmc' in the terminal to execute.")