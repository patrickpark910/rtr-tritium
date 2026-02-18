import openmc
import math


def build_triga():

    # ==============================================================================
    # 1. MATERIALS
    # ==============================================================================
    # ... (Previous materials: Fuel, Clad, Zr Rod, Water, Graphite) ...
    fuel = openmc.Material(name='UZrH Fuel')
    fuel.set_density('g/cm3', 6.0)
    fuel.add_element('U', 0.085, percent_type='wo', enrichment=20.0)
    fuel.add_element('Zr', 0.88, percent_type='wo')
    fuel.add_element('H', 0.035, percent_type='wo')
    fuel.add_s_alpha_beta('c_H_in_ZrH')

    clad = openmc.Material(name='SS304 Clad')
    clad.set_density('g/cm3', 7.9)
    clad.add_element('Fe', 0.70, percent_type='wo')
    clad.add_element('Cr', 0.19, percent_type='wo')
    clad.add_element('Ni', 0.11, percent_type='wo')

    zr_rod = openmc.Material(name='Zr Rod')
    zr_rod.set_density('g/cm3', 6.5)
    zr_rod.add_element('Zr', 1.0)

    water = openmc.Material(name='Light Water')
    water.set_density('g/cm3', 1.0)
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')

    graphite = openmc.Material(name='Graphite')
    graphite.set_density('g/cm3', 1.7)
    graphite.add_element('C', 1.0)
    graphite.add_s_alpha_beta('c_Graphite')

    # --- NEW MATERIALS for Lazy Susan ---
    aluminum = openmc.Material(name='Aluminum Rack')
    aluminum.set_density('g/cm3', 2.7)
    aluminum.add_element('Al', 1.0)

    air = openmc.Material(name='Air')
    air.set_density('g/cm3', 0.0012)
    air.add_element('N', 0.78)
    air.add_element('O', 0.21)
    air.add_element('Ar', 0.01)

    materials = openmc.Materials([fuel, clad, zr_rod, water, graphite, aluminum, air])
    materials.export_to_xml()

    # ==============================================================================
    # 2. GEOMETRY
    # ==============================================================================

    # --- Core Dimensions ---
    r_zr_rod = 0.3175
    r_fuel   = 1.82
    r_clad   = 1.87
    pitch    = 4.5

    # --- Standard Fuel Pin (Same as before) ---
    surf_zr_rod = openmc.ZCylinder(r=r_zr_rod)
    surf_fuel   = openmc.ZCylinder(r=r_fuel)
    surf_clad   = openmc.ZCylinder(r=r_clad)

    cell_center = openmc.Cell(fill=zr_rod, region=-surf_zr_rod)
    cell_fuel   = openmc.Cell(fill=fuel, region=+surf_zr_rod & -surf_fuel)
    cell_clad   = openmc.Cell(fill=clad, region=+surf_fuel & -surf_clad)
    cell_mod    = openmc.Cell(fill=water, region=+surf_clad)
    univ_fuel_pin = openmc.Universe(cells=[cell_center, cell_fuel, cell_clad, cell_mod])

    cell_all_water = openmc.Cell(fill=water)
    univ_water = openmc.Universe(cells=[cell_all_water])

    # --- Core Lattice ---
    lattice = openmc.HexLattice()
    lattice.center = (0.0, 0.0)
    lattice.pitch = (pitch,)
    lattice.outer = univ_water
    # ringW = [univ_water]*42
    ringG = [univ_fuel_pin]*36
    ringF = [univ_fuel_pin]*30
    ringE = [univ_fuel_pin]*24
    ringD = [univ_fuel_pin]*18
    ringC = [univ_fuel_pin]*12
    ringB = [univ_fuel_pin]*6
    ringA = [univ_fuel_pin]*1
    lattice.universes = [ringG, ringF, ringE, ringD, ringC, ringB, ringA] # ringW, 
    lattice.orientation = 'x'

    # ==============================================================================
    # NEW SECTION: Rotary Specimen Rack (Lazy Susan)
    # ==============================================================================
    # The rack is a circular array of 40 tubes in the reflector.
    # We model this by creating a Universe filled with air (the "watertight assembly")
    # and placing the 40 tubes into it.

    # 1. Define the Specimen Tube Geometry
    # Inner Diameter ~ 25mm -> Radius = 1.25 cm
    # Assume Wall Thickness ~ 1mm -> Outer Radius = 1.35 cm
    r_tube_inner = 1.25
    r_tube_outer = 1.35

    surf_tube_in  = openmc.ZCylinder(r=r_tube_inner)
    surf_tube_out = openmc.ZCylinder(r=r_tube_outer)

    # Tube Contents: Air/Void (Empty for now)
    cell_tube_void = openmc.Cell(fill=air, region=-surf_tube_in)
    # Tube Wall: Aluminum
    cell_tube_wall = openmc.Cell(fill=aluminum, region=+surf_tube_in & -surf_tube_out)
    # Outside Tube: Air (Inside the rack housing)
    cell_tube_ext  = openmc.Cell(fill=air, region=+surf_tube_out)

    univ_sample_tube = openmc.Universe(cells=[cell_tube_void, cell_tube_wall, cell_tube_ext])


    # 2. Place 40 Tubes in a Ring
    ls_radius = 34.0  # Radius where the rack sits (cm)
    n_positions = 40
    rack_cells = []

    # We create the background "Air" cell for the rack universe
    # We will cut holes in this air cell to place the tubes
    rack_region = openmc.ZCylinder(r=60.0) # Arbitrary large bound for this universe
    background_air_cell = openmc.Cell(fill=air, region=-rack_region)
    rack_universe = openmc.Universe(cells=[background_air_cell])

    # Loop to generate 40 tubes
    for i in range(n_positions):
        angle = 2 * math.pi * i / n_positions
        x_pos = ls_radius * math.cos(angle)
        y_pos = ls_radius * math.sin(angle)
        
        # Create a copy of the tube universe shifted to this position
        # OpenMC requires a 'fill' cell to place a universe
        
        # Define the hole for this specific tube position
        # We must use a cylinder centered at x_pos, y_pos
        cyl_surf = openmc.ZCylinder(x0=x_pos, y0=y_pos, r=r_tube_outer)
        
        # Cell containing the tube universe
        c_tube = openmc.Cell(fill=univ_sample_tube)
        c_tube.region = -cyl_surf
        c_tube.translation = (x_pos, y_pos, 0)
        
        # Add this cell to the rack universe
        rack_universe.add_cell(c_tube)
        
        # Update background air cell to exclude this region
        background_air_cell.region &= +cyl_surf


    # ==============================================================================
    # ROOT GEOMETRY
    # ==============================================================================

    # Surfaces defining radial regions
    surf_core_bound = openmc.ZCylinder(r=30.0)  # Core/Lattice boundary
    surf_rack_inner = openmc.ZCylinder(r=32.0)  # Start of Lazy Susan Well
    surf_rack_outer = openmc.ZCylinder(r=36.0)  # End of Lazy Susan Well
    surf_refl_outer = openmc.ZCylinder(r=40.0)  # Edge of graphite reflector
    surf_pool_bound = openmc.ZCylinder(r=60.0, boundary_type='vacuum') # edge of pool water in model

    # 1. Central Core Region
    cell_core = openmc.Cell(name='Core Lattice')
    cell_core.fill = lattice
    cell_core.region = -surf_core_bound

    # 2. Inner Graphite Reflector (Between core and rack)
    cell_refl_inner = openmc.Cell(name='Graphite Inner')
    cell_refl_inner.fill = graphite
    cell_refl_inner.region = +surf_core_bound & -surf_rack_inner

    # 3. The Lazy Susan Well (Filled with our Rack Universe)
    cell_rack = openmc.Cell(name='Lazy Susan Rack')
    cell_rack.fill = rack_universe
    cell_rack.region = +surf_rack_inner & -surf_rack_outer

    # 4. Outer Graphite Reflector
    cell_refl_outer = openmc.Cell(name='Graphite Outer')
    cell_refl_outer.fill = graphite
    cell_refl_outer.region = +surf_rack_outer & -surf_refl_outer

    cell_pool = openmc.Cell(name='Reactor Pool', fill=water)
    cell_pool.region = +surf_refl_outer & -surf_pool_bound

    geometry = openmc.Geometry(root=[cell_core, cell_refl_inner, cell_rack, cell_refl_outer, cell_pool])
    geometry.export_to_xml()

    # Settings (Same as before)
    settings = openmc.Settings()
    settings.batches = 50
    settings.inactive = 10
    settings.particles = 1000
    settings.source = openmc.IndependentSource(space=openmc.stats.Point())
    settings.export_to_xml()

    print("TRIGA model with Lazy Susan generated.")


    # ==============================================================================
    # 4. PLOTTING
    # ==============================================================================

    # Create a plot of the XY plane (Top-down view)
    plot = openmc.Plot()
    plot.filename = 'triga_geometry_plot'
    plot.basis = 'xy'
    plot.origin = (0, 0, 0)
    plot.width = (82, 82)       # 65x65 cm view to see the 30cm radius reflector
    plot.pixels = (1640, 1640)      # High resolution

    # Color the materials for easier identification
    plot.color_by = 'material'
    plot.colors = {
        fuel: 'orange',
        clad: 'black',
        water: 'cyan',
        graphite: 'gray',
        aluminum: 'green',        # The Rack Housing
        air: 'white',              # Inside the sample tubes
        zr_rod: 'red'
    }

    # Add the plot to a Plots collection and export
    plots = openmc.Plots([plot])
    plots.export_to_xml()

    # Run the plot generation
    openmc.plot_geometry()

if __name__ == "__main__":
    build_triga()
