import openmc
import math


def build_triga():

    # ==================================================================
    # MATERIALS
    # ==================================================================

    # U-ZrH fuel meat (10 wt% U in U-ZrH / 20 wt% enriched U = "10 by 20" fuel)
    # Suppose this fuel is 10 wt% U, 90 wt% ZrH...
    # ...of 10 wt% U, we have 20 wt% enrichment = 2 wt% U-235 and 8 wt% U-238
    # ...of 90 wt% ZrH, we have an atom ratio H/Zr of 1.6 = 1.74 wt%/98.26 wt% = 1.57 wt% H, 88.43 wt% Zr in the fuel
    fuel = openmc.Material(name='Fuel')
    fuel.set_density('g/cm3', 6.0)
    fuel.add_nuclide('U235', 0.02, percent_type='wo')
    fuel.add_nuclide('U238', 0.08, percent_type='wo')
    fuel.add_element('Zr', 0.8843, percent_type='wo')
    fuel.add_element('H', 0.0157, percent_type='wo')
    fuel.temperature = 600.0  # [K] -- ranges 600-800 K, peak 1000 K centerline
    fuel.add_s_alpha_beta('c_H_in_ZrH')
    fuel.add_s_alpha_beta('c_Zr_in_ZrH')  # Using H in H2O S(a,b) data for the remaining H in ZrH, which is a common approximation
    # Normally, you would define the fuel using openmc.Material.mix_materials
    # but OpenMC currently does not support mixing materials with S(a,b) data.
    # So, we will define the fuel composition directly.

    # Zr rod - runs through center of fuel element
    zr_rod = openmc.Material(name='Zr rod')
    zr_rod.set_density('g/cm3', 6.5)
    zr_rod.add_element('Zr', 1.0)

    # Light water
    water = openmc.Material(name='Water')
    water.set_density('g/cm3', 1.0)
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')

    """ Graphite """
    # Pure carbon has density 2.2 g/cm3, but reactor graphite is about 1.6 g/cm3
    # So 1.6/2.2 = 0.72 = ~30% porosity
    graphite = openmc.Material(name='Graphite')
    graphite.set_density('g/cm3', 1.6)
    graphite.add_element('C', 1.0)
    graphite.add_s_alpha_beta('c_Graphite_20p')  # density 1.6/2.2 = 0.72 = 30% porosity

    # ss304
    ss304 = openmc.Material(name='SS304')
    ss304.set_density('g/cm3', 7.8)
    ss304.add_element('Fe', 0.7, percent_type='wo')  # from PNNL-15870
    ss304.add_element('Cr', 0.2, percent_type='wo')
    ss304.add_element('Ni', 0.1, percent_type='wo')
    # ss304.add_s_alpha_beta('c_Fe56')  

    # Air
    air = openmc.Material(name='Air')
    air.set_density('g/cm3', 0.0012)
    air.add_element('N',  0.755, percent_type='wo')  # from PNNL-15870
    air.add_element('O',  0.232, percent_type='wo')
    air.add_element('Ar', 0.013, percent_type='wo')

    # Helium-3
    helium = openmc.Material(name='Helium-3')
    helium.set_density('g/cm3', 0.008375) # at 1000 psi  # 0.00016 g/cm3 at 21 C
    helium.add_nuclide('He3', 1.0)

    # Add materials
    materials = openmc.Materials([fuel, zr_rod, water, graphite, ss304, air, helium])
    materials.export_to_xml()


    # ==================================================================
    # GEOMETRY
    # ==================================================================

    # ------------------------------------------------------------------
    # Surfaces (1 inch = 1" = 2.54 cm)
    # ------------------------------------------------------------------

    """ Radial surfaces """
    r_zr        = 0.3              # [cm]  Zr pin      -- OD 0.25" = outer radius 0.3175 cm
    r_meat      = 1.9              # [cm]  Meat        -- OD 1.40" = outer radius 1.778 cm
    r_clad      = r_meat + 0.05    # [cm]  Clad        -- OD 1.50" = outer radius 1.905 cm
    r_tube_in   = r_meat           # [cm]  Guide tube  
    r_tube_out  = r_clad           # [cm]  Guide tube  
    r_rack_in   = 29.0             # [cm]  Rotary rack -- ID 23.8" = inner radius 30.226 cm 
    r_rack_out  = 34.0             # [cm]  Rotary rack -- OD 28.9" = outer radius 36.703 cm 
    r_refl_in   = 24.0             # [cm]  Reflector   -- ID 17.9" = inner radius 22.733 cm
    r_refl_out  = 48.0             # [cm]  Reflector   -- OD 41.8" = outer radius 53.086 cm
    r_tube_in   = 1.25             # [cm]  Sample tube 
    r_tube_out  = r_tube_in + 0.25 # [cm]  Sample tube 
    r_vial_in   = 1.15             # [cm]  Sample vial 
    r_vial_out  = r_vial_in + 0.10 # [cm]  Sample vial
    r_void      = 100.0            # [cm]  Model boundary
        
    """ Axial surfaces """
    z_meat_bot = 0.0                # [cm]  Meat bot
    z_meat_top = 38.0               # [cm]  Meat top -- H 15.0" = 38.1 cm
    z_grph_bot = z_meat_bot - 9.0   # [cm]  Graphite cap bot -- H 3.7" = 9.398 cm  
    z_grph_top = z_meat_top + 6.0   # [cm]  Graphite cap top -- H 2.6" = 6.604 cm
    z_pin_bot  = z_grph_bot - 0.05  # [cm]  Fuel pin bot -- Thk 0.02" = 0.0508 cm
    z_pin_top  = z_grph_top + 0.05  # [cm]  Fuel pin top -- Thk 0.02" = 0.0508 cm
    z_refl_bot = z_pin_bot          # [cm]  Reflector bot 
    z_refl_top = z_pin_top          # [cm]  Reflector top 
    z_rack_top = z_refl_top         # [cm]  Rotary rack top -- 0.02" = 0.0508 cm thick on all sides
    z_rack_bot = z_refl_top - 25.0  # [cm]  Rotary rack bot -- H 10.2" = 25.908 cm
    z_vial_bot = z_rack_bot         # [cm]  Sample tube bot
    z_vial_top = z_vial_bot + 15.1  # [cm]  Sample tube top
    z_ct_bot   = 19.0 - (15.1/2.0)  # [cm]  Central thimble bot
    z_ct_top   = 19.0 + (15.1/2.0)  # [cm]  Central thimble top
    z_grid_bot = z_pin_bot - 0.5    # [cm]  Grid plate bot -- H 0.75" = 1.905 cm
    z_grid_top = z_pin_top + 0.5    # [cm]  Grid plate top -- H 0.75" = 1.905 cm


    
    """ Radial surface definitions """

    cz_zr       = openmc.ZCylinder(r=r_zr)
    cz_meat     = openmc.ZCylinder(r=r_meat)
    cz_clad     = openmc.ZCylinder(r=r_clad)
    cz_tube_in  = openmc.ZCylinder(r=r_tube_in )
    cz_tube_out = openmc.ZCylinder(r=r_tube_out)
    cz_rack_in  = openmc.ZCylinder(r=r_rack_in )
    cz_rack_out = openmc.ZCylinder(r=r_rack_out)
    cz_refl_in  = openmc.ZCylinder(r=r_refl_in )
    cz_refl_out = openmc.ZCylinder(r=r_refl_out)
    cz_vial_in  = openmc.ZCylinder(r=r_vial_in)
    cz_vial_out = openmc.ZCylinder(r=r_vial_out)

    """ Axial surface definitions """
    pz_pin_bot  = openmc.ZPlane(z0=z_pin_bot)
    pz_grph_bot = openmc.ZPlane(z0=z_grph_bot)
    pz_meat_bot = openmc.ZPlane(z0=z_meat_bot)
    pz_meat_top = openmc.ZPlane(z0=z_meat_top)
    pz_grph_top = openmc.ZPlane(z0=z_grph_top)
    pz_pin_top  = openmc.ZPlane(z0=z_pin_top)
    pz_refl_bot = openmc.ZPlane(z0=z_refl_bot)
    pz_refl_top = openmc.ZPlane(z0=z_refl_top)
    pz_rack_bot = openmc.ZPlane(z0=z_rack_bot)
    pz_rack_top = openmc.ZPlane(z0=z_rack_top)
    pz_grid_bot = openmc.ZPlane(z0=z_grid_bot)
    pz_grid_top = openmc.ZPlane(z0=z_grid_top)
    pz_vial_bot = openmc.ZPlane(z0=z_vial_bot)
    pz_vial_top = openmc.ZPlane(z0=z_vial_top)
    pz_ct_bot   = openmc.ZPlane(z0=z_ct_bot)
    pz_ct_top   = openmc.ZPlane(z0=z_ct_top)


    # ------------------------------------------------------------------
    # Cells
    # ------------------------------------------------------------------

    cell_cap_bot   = openmc.Cell(fill=ss304, region=-cz_clad & +pz_pin_bot & -pz_grph_bot)
    cell_bot_gr    = openmc.Cell(fill=graphite, region=-cz_meat & +pz_grph_bot & -pz_meat_bot)
    cell_bot_clad  = openmc.Cell(fill=ss304, region=+cz_meat & -cz_clad & +pz_grph_bot & -pz_meat_bot)
    cell_fuel_zr   = openmc.Cell(fill=zr_rod,   region=-cz_zr & +pz_meat_bot & -pz_meat_top)
    cell_fuel_meat = openmc.Cell(fill=fuel,     region=+cz_zr & -cz_meat & +pz_meat_bot & -pz_meat_top)
    cell_fuel_clad = openmc.Cell(fill=ss304, region=+cz_meat & -cz_clad & +pz_meat_bot & -pz_meat_top)
    cell_top_gr    = openmc.Cell(fill=graphite, region=-cz_meat & +pz_meat_top & -pz_grph_top)
    cell_top_clad  = openmc.Cell(fill=ss304, region=+cz_meat & -cz_clad & +pz_meat_top & -pz_grph_top)
    cell_cap_top   = openmc.Cell(fill=ss304, region=-cz_clad & +pz_grph_top & -pz_pin_top)

    # Water Surround (Fills everything in the universe outside the pin)
    pin_region = -cz_clad & +pz_pin_bot & -pz_pin_top
    cell_water = openmc.Cell(fill=water, region=~pin_region)

    univ_fuel_pin = openmc.Universe(cells=[
        cell_cap_bot, cell_bot_gr, cell_bot_clad, 
        cell_fuel_zr, cell_fuel_meat, cell_fuel_clad, 
        cell_top_gr, cell_top_clad, cell_cap_top, 
        cell_water,
        ])

    univ_water = openmc.Universe(cells=[openmc.Cell(fill=water)])


    # ------------------------------------------------------------------
    # Guide tube universe
    # ------------------------------------------------------------------

    # We will use the same axial Z-bounds as the fuel pin to keep it neat
    guide_region = -cz_tube_out & +pz_pin_bot & -pz_pin_top

    # Cells defining the ss304 guide tube
    cell_guide_in  = openmc.Cell(fill=water,    region=-cz_tube_in & +pz_pin_bot & -pz_pin_top)
    cell_guide_wall  = openmc.Cell(fill=ss304, region=+cz_tube_in & -cz_tube_out & +pz_pin_bot & -pz_pin_top)
    
    # Pool water filling everything outside the tube in this universe
    cell_guide_out = openmc.Cell(fill=water, region=~guide_region)

    univ_guide_tube = openmc.Universe(cells=[cell_guide_in , cell_guide_wall, cell_guide_out])


    # ------------------------------------------------------------------
    # Central He-3 Holder Universe (Ring A)
    # ------------------------------------------------------------------
    # The 8cm He-3 region
    cell_ct = openmc.Cell(fill=helium, region=-cz_vial_in & +pz_ct_bot & -pz_ct_top)
    
    # Fill the rest of the lattice position with water
    cell_ct_water = openmc.Cell(fill=water, region=~(cell_ct.region))
    
    univ_ct = openmc.Universe(cells=[cell_ct, cell_ct_water])

    # ------------------------------------------------------------------
    # Core concentric ring universe
    # ------------------------------------------------------------------
    core_universe = openmc.Universe(name='Concentric Core')
    
    # Fill the background of the core with pool water. 
    # We will "punch holes" in this water for each pin location.
    core_bg_cell = openmc.Cell(name='Core Background', fill=water)
    core_universe.add_cell(core_bg_cell)
    
    pitch = 4.2       # [cm] Standard TRIGA pin-to-pin pitch
    pin_radius = 1.95  # [cm] Bounding radius for the translation cell (must be > clad r=1.905 and < pitch/2)

    def place_element(universe, x, y):
        # Create the bounding cylinder directly at the (x, y) coordinates
        pin_boundary = openmc.ZCylinder(x0=x, y0=y, r=pin_radius)
        
        # Fill the region inside this boundary with the target universe
        cell = openmc.Cell(fill=universe, region=-pin_boundary)
        
        # Translate the coordinates of the internal universe so it centers on (x,y)
        cell.translation = (x, y, 0)
        core_universe.add_cell(cell)
        
        # Punch a hole in the core background water at this exact location
        if core_bg_cell.region is None:
            core_bg_cell.region = +pin_boundary
        else:
            core_bg_cell.region &= +pin_boundary

    # Ring A (Center)
    place_element(univ_ct, 0.0, 0.0)

    # Define the number of elements per ring (TRIGA standard: 6, 12, 18, 24, 30)
    ring_counts = [6, 12, 18, 24, 30]
    
    # Maintain your exact layout of guide tubes vs fuel pins
    ring_elements = {
        0: [univ_fuel_pin] * 6,                                                                                # Ring B
        1: [univ_fuel_pin]*3 + [univ_guide_tube] + [univ_fuel_pin]*5 + [univ_guide_tube] + [univ_fuel_pin]*2,  # Ring C
        2: ([univ_guide_tube] + [univ_fuel_pin]*8) * 2,                                                        # Ring D
        3: [univ_fuel_pin] * 24,                                                                               # Ring E
        4: [univ_fuel_pin] * 30                                                                                # Ring F
    }

    # Loop through and calculate polar coordinates for each pin
    for ring_idx, count in enumerate(ring_counts):
        radius = (ring_idx + 1) * pitch
        elements = ring_elements[ring_idx]
        
        for i in range(count):
            angle = 2 * math.pi * i / count
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            place_element(elements[i], x, y)

    # ------------------------------------------------------------------
    # Rotary rack universe
    # ------------------------------------------------------------------
    univ_sample_tube = openmc.Universe()
    # Holes are 38 mm (1.9 cm radius). Using 1.8 cm for inner radius to simulate 1mm wall.
    cell_tube_wall = openmc.Cell(fill=ss304, region=+openmc.ZCylinder(r=1.8) & -openmc.ZCylinder(r=1.9))
    
    # He-3 sample region
    cell_sample    = openmc.Cell(fill=helium, region=-cz_vial_in & +pz_vial_bot & -pz_vial_top) # fill=helium
    
    # SS304 radial cladding around the He-3
    cell_vial_clad = openmc.Cell(fill=ss304,  region=+cz_vial_in & -cz_vial_out & +pz_vial_bot & -pz_vial_top)
    
    # Define the combined volume of the vial (He-3 + cladding) to cut it out of the air cell
    vial_volume = -cz_vial_out & +pz_vial_bot & -pz_vial_top
    
    # Fill the remaining empty space inside the 1.8 cm tube with air
    cell_tube_inner_air = openmc.Cell(fill=air, region=-openmc.ZCylinder(r=1.8) & ~vial_volume)
    
    cell_tube_ext  = openmc.Cell(fill=air, region=+openmc.ZCylinder(r=1.9)) 
    
    # Make sure to include the new cell_vial_clad in the universe
    univ_sample_tube.add_cells([cell_tube_wall, cell_sample, cell_vial_clad, cell_tube_inner_air, cell_tube_ext])

    ls_radius = (r_rack_in  + r_rack_out) / 2.0
    rack_universe = openmc.Universe()
    bg_cell = openmc.Cell(fill=air, region=-openmc.ZCylinder(r=60)) 
    rack_universe.add_cell(bg_cell)

    for i in range(40):
        angle = 2 * math.pi * i / 40
        x = ls_radius * math.cos(angle)
        y = ls_radius * math.sin(angle)
        cyl = openmc.ZCylinder(x0=x, y0=y, r=1.9)
        c = openmc.Cell(fill=univ_sample_tube, region=-cyl)
        c.translation = (x, y, 0)
        rack_universe.add_cell(c)
        bg_cell.region &= +cyl

    # ==============================================================================
    # 6. ROOT GEOMETRY ASSEMBLY
    # ==============================================================================

    # Derived Regions
    region_refl_bulk = (+cz_refl_in  & -cz_refl_out & +pz_refl_bot & -pz_refl_top)
    region_ls_cutout = (+cz_rack_in  & -cz_rack_out & +pz_rack_bot & -pz_refl_top) 

    # 1-4. Solid Components
    cell_core        = openmc.Cell(name='Core',               fill=core_universe, region=-cz_refl_in & +pz_pin_bot  & -pz_pin_top)
    cell_reflector   = openmc.Cell(name='Graphite Reflector', fill=graphite,      region=region_refl_bulk & ~region_ls_cutout)
    cell_rack        = openmc.Cell(name='Lazy Susan',         fill=rack_universe, region=region_ls_cutout)
    cell_grid_bot    = openmc.Cell(name='Bottom Grid',        fill=ss304,      region=-cz_refl_in & +pz_grid_bot & -pz_pin_bot)
    cell_grid_top    = openmc.Cell(name='Top Grid',           fill=ss304,      region=-cz_refl_in & +pz_pin_top  & -pz_grid_top)

    # 5. Pool Water boundaries
    pool             = openmc.ZCylinder(r=100.0, boundary_type='vacuum')
    pool_top         = openmc.ZPlane(z0=200.0,   boundary_type='vacuum')
    pool_bot         = openmc.ZPlane(z0=-100.0,  boundary_type='vacuum')

    cell_pool        = openmc.Cell(name='Pool Water', fill=water)
    cell_pool.region = (+pool_bot & -pool_top & -pool) \
                     & ~(cell_core.region)             \
                     & ~(cell_reflector.region)        \
                     & ~(cell_rack.region)             \
                     & ~(cell_grid_bot.region)         \
                     & ~(cell_grid_top.region)

    geometry         = openmc.Geometry(root=[cell_core, cell_reflector, cell_rack, cell_grid_bot, cell_grid_top, cell_pool])
    # geometry.export_to_xml()


    # ==================================================================
    # SETTINGS
    # ==================================================================
    settings = openmc.Settings()
    settings.batches = 50
    settings.inactive = 10
    settings.particles = int(1e5)
    # Define a spatial bounding box for the starting source particles inside the core
    bounds = [-20, -20, z_meat_bot, 20, 20, z_meat_top]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
    settings.source = openmc.IndependentSource(space=uniform_dist)
    # settings.export_to_xml()

    print("Full 3D TRIGA Geometry with Simplified Pins Generated.")


    # ==================================================================
    # TALLIES
    # ==================================================================
    tallies = openmc.Tallies()
    cell_filter = openmc.CellFilter([cell_ct, cell_sample])

    # He-3(n,p)T
    he3_tally = openmc.Tally(name='He3_np_RotaryRack')  
    he3_tally.filters = [cell_filter]
    he3_tally.nuclides = ['He3']
    he3_tally.scores = ['(n,p)']
    tallies.append(he3_tally)

    # Flux
    flux_tally = openmc.Tally(name='Flux_RotaryRack')  
    flux_tally.filters = [cell_filter]
    flux_tally.scores = ['flux']
    tallies.append(flux_tally)

    # tallies.export_to_xml()
    
    print("Tallies exported successfully.")


    # ==============================================================================
    # 4. PLOTTING
    # ==============================================================================
    
    # Define a common color scheme dictionary
    material_colors = {
        fuel: 'orange',
        water: 'cyan',
        graphite: 'gray',
        ss304: 'green',        # Rack housing, cladding, grid plates
        air: 'white',             # Inside sample tubes, lazy susan groove
        zr_rod: 'red',
        helium: 'magenta',
    }

    # --- Plot 1: XY Plane (Top-down view) ---
    plot_xy = openmc.Plot()
    plot_xy.filename = 'triga_xy_plot'
    plot_xy.basis = 'xy'
    # Z=19.05 cm is the vertical midplane of the 38.1 cm active fuel
    plot_xy.origin = (0, 0, 24) 
    plot_xy.width = (110, 110)      # 110cm wide to see the full 53cm radius reflector
    plot_xy.pixels = (2200, 2200)   # High resolution
    plot_xy.color_by = 'material'
    plot_xy.colors = material_colors

    # --- Plot 2: XZ Plane (Side view) ---
    plot_xz = openmc.Plot()
    plot_xz.filename = 'triga_xz_plot'
    plot_xz.basis = 'xz'
    # Center on the core midplane
    plot_xz.origin = (0, 0, 19.05)
    plot_xz.width = (110, 110)      # 110cm wide and tall to see the axial reflectors and grids
    plot_xz.pixels = (2200, 2200)
    plot_xz.color_by = 'material'
    plot_xz.colors = material_colors

    # Add both plots to a Plots collection and export
    plots = openmc.Plots([plot_xy, plot_xz])
    plots.export_to_xml()


    # ==============================================================================
    # MODEL EXPORT
    # ==============================================================================
    
    # Combine all objects into a single OpenMC Model
    model = openmc.Model(
        materials=materials,
        geometry=geometry,
        settings=settings,
        tallies=tallies,
        plots=plots
    )
    
    # Export everything to a single model.xml file
    model.export_to_xml()


    # Run the plot generation
    openmc.plot_geometry()

if __name__ == "__main__":
    build_triga()
