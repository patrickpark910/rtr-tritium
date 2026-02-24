import os, glob
import openmc
import openmc.deplete
import numpy as np
import matplotlib.pyplot as plt

from python.utilities import *

def build_mtr(pattern='R', he3_pressure=1000.0, he3_units='psi'):

    # ================================================================================
    # MATERIALS
    # ================================================================================

    # U3Si2 (20% enriched)
    fuel = openmc.Material(name='U3Si2-Al')
    fuel.set_density('g/cm3', 6.732)
    fuel.add_nuclide('U235', 0.713 * 0.1975, 'wo')
    fuel.add_nuclide('U238', 0.713 * 0.8025, 'wo')
    fuel.add_element('Si', 0.056, 'wo')
    fuel.add_element('Al', 0.231, 'wo')
    fuel.temperature = 600 # [K]
    fuel.depletable = True
    fuel.volume = 7821 # 0.0550*7.9*60*20*15=

    # Aluminum
    clad = openmc.Material(name='Aluminum')
    clad.set_density('g/cm3', 2.7)
    clad.add_element('Al', 1.0)
    clad.depletable = False

    # Light water
    water = openmc.Material(name='Water')
    water.set_density('g/cm3', 1.0)
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')
    water.depletable = False

    if pattern not in ['R']:
        steel = openmc.Material(name='316L stainless steel')
        steel.set_density('g/cm3', 8.00)
        steel.add_element('Fe', 65.5, percent_type='wo')
        steel.add_element('Cr', 17.0, percent_type='wo')
        steel.add_element('Ni', 12.0, percent_type='wo')
        steel.add_element('Mo', 2.5,  percent_type='wo')
        steel.add_element('Mn', 2.0,  percent_type='wo')
        steel.add_element('Si', 1.0,  percent_type='wo')

        he3_gcm3 = rho_he3(T=300.0, P=he3_pressure, units=he3_units)
        helium = openmc.Material(name='Helium-3')
        helium.set_density('g/cm3', he3_gcm3)
        helium.add_element('He', 1.0)
        helium.depletable = True
        helium.volume = 5890.5 # cm^3 -- volume of ALL helium in the core for depletion

        materials = openmc.Materials([fuel, clad, water, steel, helium])

    else:
        materials = openmc.Materials([fuel, clad, water])

    # ================================================================================
    # GEOMETRY
    # ================================================================================

    # --------------------------------------------------------------------------------
    # Dimensions
    # --------------------------------------------------------------------------------

    # Core dimensions
    n_fa_per_side = 5

    # Fuel assembly dimensions
    n_plates = 20
    fa_width_active = 7.9
    fa_width  =  8.0
    fa_height = 60.0 

    # Fuel plate dimensions
    tk_meat = 0.0550
    tk_clad = 0.0375
    tk_cool = (fa_width_active - n_plates * (tk_meat + tk_clad*2)) / n_plates / 2  # = 0.1325 cm for 20 plates -- thickness of water on EACH SIDE of plate
    fp_pitch = tk_meat + 2*tk_clad + 2*tk_cool  # fp_pitch * n_plates = fa_width_active = 7.9 cm -- total width of all plates + channels = fuel assembly width minus structure on either side

    # Helium tank dimensions
    tk_steel = 1.0
    r_helium = 2.5 
    r_steel  = r_helium + tk_steel  # thickness 1 cm of austenitic 316L stainless steel

    # Pool
    r_pool = 100.0
    hh_pool = r_pool # half-height
    

    # --------------------------------------------------------------------------------
    # Surfaces
    # --------------------------------------------------------------------------------

    px_meat_l = openmc.XPlane(x0= -tk_meat/2)
    px_meat_r = openmc.XPlane(x0= +tk_meat/2)

    px_clad_l = openmc.XPlane(x0= -tk_clad)
    px_clad_r = openmc.XPlane(x0= +tk_clad)

    px_pitch_l = openmc.XPlane(x0= -fp_pitch/2)
    px_pitch_r = openmc.XPlane(x0= +fp_pitch/2)

    # RectangularPrism is an _infinite_ rectangular prism bounded by four planar surfaces.
    rpp_fa_in  = openmc.model.RectangularPrism(width=fa_width_active, height=fa_width_active)
    rpp_fa_out = openmc.model.RectangularPrism(width=fa_width, height=fa_width)

    pz_core_bot = openmc.ZPlane(z0= -fa_height/2)
    pz_core_top = openmc.ZPlane(z0= +fa_height/2)

    pz_steel_bot = openmc.ZPlane(z0= -fa_height/2+tk_steel)
    pz_steel_top = openmc.ZPlane(z0= +fa_height/2-tk_steel)

    pz_pool_bot = openmc.ZPlane(z0= -hh_pool)
    pz_pool_top = openmc.ZPlane(z0= +hh_pool)

    rpp_core = openmc.model.RectangularPrism(width= fa_width*n_fa_per_side, height=fa_width*n_fa_per_side)
    rcc_pool = openmc.model.RightCircularCylinder( (0, 0, -hh_pool), 2*hh_pool, r_pool, boundary_type='vacuum')

    rcc_helium = openmc.model.RightCircularCylinder( (0, 0, -hh_pool), 2*hh_pool, r_helium)
    rcc_steel   = openmc.model.RightCircularCylinder( (0, 0, -hh_pool), 2*hh_pool, r_steel)


    # --------------------------------------------------------------------------------
    # Cells
    # --------------------------------------------------------------------------------

    # Fuel plate
    cell_meat  = openmc.Cell(fill=fuel, region= +px_meat_l & -px_meat_r)
    cell_clad  = openmc.Cell(fill=clad, region=(+px_clad_l & -px_meat_r) | (+px_meat_l & -px_clad_r))
    cell_cool  = openmc.Cell(fill=water, region=(+px_pitch_l & -px_clad_l) | (+px_clad_r & -px_pitch_r))
    univ_plate = openmc.Universe(cells=[cell_meat, cell_clad, cell_cool])

    # Fuel assembly
    lat_fa = openmc.RectLattice(name='Fuel assembly lattice')
    lat_fa.lower_left = [-fa_width_active/2, -fa_width_active/2]
    lat_fa.pitch      = [fp_pitch, fa_width_active] 
    lat_fa.universes  = [[univ_plate] * n_plates]

    cell_active = openmc.Cell(fill=lat_fa, region= -rpp_fa_in)  
    cell_fa     = openmc.Cell(fill=clad,   region= +rpp_fa_in & -rpp_fa_out & +pz_core_bot & -pz_core_top)  # Structural region around the active fuel, filled with clad material for simplicity
    cell_cool   = openmc.Cell(fill=water,  region= +rpp_fa_out & +pz_core_bot & -pz_core_top)  # Coolant region outside the fuel assembly
    univ_fa     = openmc.Universe(cells=[cell_active, cell_fa, cell_cool])  

    if pattern not in ['R']:
        # Helium tank
        cell_helium = openmc.Cell(fill=helium, region= -rcc_helium & +pz_steel_bot & -pz_steel_top)
        cell_steel  = openmc.Cell(fill=steel,  region= +rcc_helium & -rcc_steel    | (-rcc_steel & (-pz_steel_bot | +pz_steel_top)))
        cell_water  = openmc.Cell(fill=water,  region= +rcc_steel | (+pz_steel_top | -pz_steel_bot))
        univ_he = openmc.Universe(cells=[cell_helium, cell_steel, cell_water])

    # Core 
    lat_core = openmc.RectLattice(name='Core lattice')
    lat_core.lower_left = [-fa_width * (n_fa_per_side/2), -fa_width * (n_fa_per_side/2)]
    lat_core.pitch = [fa_width, fa_width]

    lat_core.universes = [ [univ_fa, univ_fa, univ_fa, univ_fa, univ_fa],
                           [univ_fa, univ_fa, univ_fa, univ_fa, univ_fa],
                           [univ_fa, univ_fa, univ_fa, univ_fa, univ_fa],
                           [univ_fa, univ_fa, univ_fa, univ_fa, univ_fa],
                           [univ_fa, univ_fa, univ_fa, univ_fa, univ_fa],
                            ]

    if pattern == 'A':
        lat_core.universes = [ [univ_fa, univ_fa, univ_fa, univ_fa, univ_fa],
                               [univ_fa, univ_he, univ_fa, univ_he, univ_fa],
                               [univ_fa, univ_fa, univ_he, univ_fa, univ_fa],
                               [univ_fa, univ_he, univ_fa, univ_he, univ_fa],
                               [univ_fa, univ_fa, univ_fa, univ_fa, univ_fa], ]


    cell_core = openmc.Cell(fill=lat_core, region= -rpp_core & +pz_core_bot & -pz_core_top)
    cell_refl = openmc.Cell(fill=water, region= -rcc_pool & (+rpp_core | -pz_core_bot | +pz_core_top) & +pz_pool_bot & -pz_pool_top )  # Reflector region outside the core, filled with water for simplicity

    univ_root = openmc.Universe(cells=[cell_core, cell_refl])
    geometry = openmc.Geometry(univ_root)


    # ================================================================================
    # SETTINGS
    # ================================================================================

    settings = openmc.Settings()
    settings.batches = 120
    settings.inactive = 20
    settings.particles = int(1e5)

    # Set a bounding box spatial distribution for the initial source
    lower_left   = [-fa_width*(n_fa_per_side/2), -fa_width*(n_fa_per_side/2), -fa_height/2]
    upper_right  = [+fa_width*(n_fa_per_side/2), +fa_width*(n_fa_per_side/2), +fa_height/2]
    uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=True)
    settings.source = openmc.IndependentSource(space=uniform_dist)


    # ================================================================================
    # TALLIES
    # ================================================================================

    tallies = openmc.Tallies()

    # Tally 1: Total (n,p) reactions and flux (neutron population) at each depth
    tally_spatial = openmc.Tally(name='Spatial tally')
    tally_spatial.scores = ['flux', '(n,p)']
    tallies.append(tally_spatial)

    # Tally 2: Fine energy grid (n,p) reactions at each depth for histogram
    # 500 logarithmically spaced bins from 1e-5 eV to 10 eV
    energy_bins = np.logspace(-5, 1, 501)
    filter_energy = openmc.EnergyFilter(energy_bins)

    tally_energy = openmc.Tally(name='Energy-spatial tally')
    tally_energy.filters = [filter_energy]
    tally_energy.scores = ['flux', '(n,p)']
    tallies.append(tally_energy)


    # ================================================================================
    # PLOTTING
    # ================================================================================

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
    # Zoomed to 50x50 cm to nicely frame the 30x30 cm 3x3 core
    plot_xy.width = (50, 50)      
    plot_xy.pixels = (6000, 6000)   # High resolution needed to resolve thin plates
    plot_xy.color_by = 'material'
    plot_xy.colors = material_colors

    # --- Plot 2: XZ Plane (Side view) ---
    plot_xz = openmc.Plot()
    plot_xz.filename = 'mtr_xz_plot'
    plot_xz.basis = 'xz'
    # Center exactly on the core midplane
    plot_xz.origin = (0, 0, 0)
    # 50 cm wide to show the core, 80 cm tall to capture the 60 cm active height + axial reflectors
    plot_xz.width = (50, 80)      
    plot_xz.pixels = (1000, 2000)   # Maintain aspect ratio
    plot_xz.color_by = 'material'
    plot_xz.colors = material_colors

    # Add both plots to a Plots collection and export
    plots = openmc.Plots([plot_xy, plot_xz])

    return openmc.model.Model(geometry=geometry, materials=materials, settings=settings, tallies=tallies, plots=plots) 


def plot_np_tallies(path):
    # 1. Find all statepoint files in the directory and sort them sequentially
    search_pattern = os.path.join(path, "openmc_simulation_n*.h5")
    statepoints = sorted(glob.glob(search_pattern))
    
    if not statepoints:
        print(f"No statepoint files found in {path}")
        return

    np_rates = []
    np_errors = []  # Added list to store standard deviations
    
    # 2. Iterate through each depletion step
    for sp_file in statepoints:
        with openmc.StatePoint(sp_file) as sp:
            try:
                # Retrieve the specific tally by the name defined in mtr.py
                tally = sp.get_tally(name='Spatial tally')
                
                # Convert the tally data to a Pandas DataFrame for easy filtering
                df = tally.get_pandas_dataframe()
                
                # Filter out the 'flux' score, isolate '(n,p)'
                df_np = df[df['score'] == '(n,p)']
                
                # Get the total mean
                step_np_rate = df_np['mean'].sum()
                np_rates.append(step_np_rate)
                
                # Calculate propagated standard deviation (square root of the sum of variances)
                step_np_err = np.sqrt((df_np['std. dev.'] ** 2).sum())
                np_errors.append(step_np_err)

            except KeyError:
                print(f"Error: 'Spatial tally' not found in {sp_file}.")
                print("Make sure 'tallies=tallies' is included in your mtr.py model export.")
                return

    steps = list(range(len(np_rates)))

    # 3. Export the extracted data to CSV using np.savetxt
    # Stack steps, rates, and errors side-by-side into a 2D array
    export_data = np.column_stack((steps, np_rates, np_errors))
    csv_filepath = os.path.join(path, "np_reaction_rates.csv")
    
    np.savetxt(
        csv_filepath, 
        export_data, 
        delimiter=",", 
        header="Depletion_Step,NP_Reaction_Rate,NP_Reaction_Error", 
        comments="",             # Prevents adding a '#' at the start of the header
        fmt=["%d", "%.8e", "%.8e"] # Format steps as integers, rates/errors as scientific notation
    )
    print(f"Data successfully exported to: {csv_filepath}")

    # 4. Plot the extracted data with error bars
    plt.figure(figsize=(8, 5))
    
    # Use errorbar instead of plot to visualize the uncertainty
    plt.errorbar(
        steps, np_rates, yerr=np_errors, 
        marker='o', linestyle='-', color='darkred', 
        linewidth=2, capsize=4, capthick=1.5
    )
    
    plt.title("Total (n,p) Reaction Rate Across Depletion Steps")
    plt.xlabel("Depletion Step")
    plt.ylabel("(n,p) Reaction Rate (reactions / source neutron)")
    plt.xticks(steps) # Ensure the x-axis only shows whole integer steps
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    plt.show()



if __name__ == "__main__":
    build_mtr()