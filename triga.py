import openmc
import math

# ==============================================================================
# 1. MATERIALS
# ==============================================================================

# --- Fuel: U-ZrH1.6 (20% enriched Uranium, 8.5 wt% U total) ---
# Note: Real TRIGA fuel is complex. We approximate with a homogeneous mixture.
fuel = openmc.Material(name='UZrH Fuel')
fuel.set_density('g/cm3', 6.0) # Approx density for standard TRIGA fuel
# Add Uranium (20% enriched)
fuel.add_element('U', 0.085, percent_type='wo', enrichment=20.0)
# Add Zirconium
fuel.add_element('Zr', 0.88, percent_type='wo') # Balance of weight
# Add Hydrogen (Stoichiometry approx ZrH1.6)
# We add H based on atomic ratio to Zr. 
# For simplicity in this basic script, we add H as a trace moderator 
# or calculate the exact mass fraction. Here we approximate:
fuel.add_element('H', 0.035, percent_type='wo') 
# Important: Add S(alpha, beta) thermal scattering for H in ZrH
fuel.add_s_alpha_beta('h_zrh')

# --- Cladding: Stainless Steel 304 ---
clad = openmc.Material(name='SS304 Clad')
clad.set_density('g/cm3', 7.9)
clad.add_element('Fe', 0.70, percent_type='wo')
clad.add_element('Cr', 0.19, percent_type='wo')
clad.add_element('Ni', 0.11, percent_type='wo')

# --- Central Rod: Zirconium ---
zr_rod = openmc.Material(name='Zr Rod')
zr_rod.set_density('g/cm3', 6.5)
zr_rod.add_element('Zr', 1.0)

# --- Moderator: Light Water ---
water = openmc.Material(name='Light Water')
water.set_density('g/cm3', 1.0)
water.add_element('H', 2.0)
water.add_element('O', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')

# --- Reflector: Graphite ---
graphite = openmc.Material(name='Graphite')
graphite.set_density('g/cm3', 1.7)
graphite.add_element('C', 1.0)
graphite.add_s_alpha_beta('c_Graphite')

# Export materials
materials = openmc.Materials([fuel, clad, zr_rod, water, graphite])
materials.export_to_xml()


# ==============================================================================
# 2. GEOMETRY
# ==============================================================================

# --- Dimensions (cm) ---
# Typical TRIGA dimensions
r_zr_rod = 0.3175   # Central Zr rod radius
r_fuel   = 1.82     # Outer radius of fuel meat
r_clad   = 1.87     # Outer radius of cladding (approx 0.05cm thickness)
pitch    = 4.0      # Pin pitch in lattice

# --- Surfaces ---
surf_zr_rod = openmc.ZCylinder(r=r_zr_rod)
surf_fuel   = openmc.ZCylinder(r=r_fuel)
surf_clad   = openmc.ZCylinder(r=r_clad)

# --- Cells for a Single Fuel Pin ---
# 1. Central Zirconium Rod
cell_center = openmc.Cell(name='Central Zr Rod')
cell_center.fill = zr_rod
cell_center.region = -surf_zr_rod

# 2. Fuel Meat (Annulus between Zr rod and Clad inner)
cell_fuel = openmc.Cell(name='Fuel Meat')
cell_fuel.fill = fuel
cell_fuel.region = +surf_zr_rod & -surf_fuel

# 3. Cladding (Annulus)
cell_clad = openmc.Cell(name='Cladding')
cell_clad.fill = clad
cell_clad.region = +surf_fuel & -surf_clad

# 4. Moderator (Outside clad)
cell_mod = openmc.Cell(name='Moderator')
cell_mod.fill = water
cell_mod.region = +surf_clad

# Create a Universe for the Fuel Pin
univ_fuel_pin = openmc.Universe(cells=[cell_center, cell_fuel, cell_clad, cell_mod])

# Create a Universe for Water (for empty lattice spots)
cell_all_water = openmc.Cell(fill=water)
univ_water = openmc.Universe(cells=[cell_all_water])

# --- Core Lattice ---
# We use a Hexagonal Lattice to simulate the core cluster.
# This creates a small core with 3 rings of fuel.
lattice = openmc.HexLattice()
lattice.center = (0.0, 0.0)
lattice.pitch = (pitch,)
lattice.outer = univ_water # Fill outside of lattice with water

# Define rings (From outside in)
# Ring 3: Water (Reflector/Edge)
# Ring 2: Fuel
# Ring 1: Fuel (Center)
ring3 = [univ_water]*12
ring2 = [univ_fuel_pin]*6
ring1 = [univ_fuel_pin]*1

lattice.universes = [ring3, ring2, ring1]
lattice.orientation = 'x'

# --- Root Geometry ---
# Fill a large cylinder with the lattice, surrounded by graphite reflector
core_region_surf = openmc.ZCylinder(r=15.0) # Boundary of the lattice region
reflector_surf   = openmc.ZCylinder(r=25.0, boundary_type='vacuum')

cell_core = openmc.Cell(name='Core Lattice')
cell_core.fill = lattice
cell_core.region = -core_region_surf

cell_reflector = openmc.Cell(name='Graphite Reflector')
cell_reflector.fill = graphite
cell_reflector.region = +core_region_surf & -reflector_surf

geometry = openmc.Geometry(root=[cell_core, cell_reflector])
geometry.export_to_xml()


# ==============================================================================
# 3. SETTINGS
# ==============================================================================

settings = openmc.Settings()
settings.batches = 50
settings.inactive = 10
settings.particles = 1000

# Define a starting source (Point source in the center)
source = openmc.IndependentSource()
source.space = openmc.stats.Point(xyz=(0, 0, 0))
settings.source = source

settings.export_to_xml()

print("Files generated: materials.xml, geometry.xml, settings.xml")
print("Run with: openmc.run()")