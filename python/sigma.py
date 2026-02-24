import openmc
import numpy as np
import matplotlib.pyplot as plt

from python.utilities import *


def build_sigma(pattern='R', he3_pressure=1000.0, he3_units='psi'):

    # ================================================================================
    # MATERIALS
    # ================================================================================

    he3_gcm3 = rho_he3(T=300.0, P=he3_pressure, units=he3_units)

    he3 = openmc.Material(name='He-3')
    he3.add_nuclide('He3', 1.0)
    he3.set_density('g/cm3', he3_gcm3)

    materials = openmc.Materials([he3])


    # ================================================================================
    # GEOMETRY
    # ================================================================================
    
    pxl = openmc.XPlane(x0= -0.5, boundary_type='reflective')
    pxr = openmc.XPlane(x0= +0.5, boundary_type='reflective')
    pyl = openmc.YPlane(y0= -0.5, boundary_type='reflective')
    pyr = openmc.YPlane(y0= +0.5, boundary_type='reflective')
    pz_min = openmc.ZPlane(z0=0.0, boundary_type='vacuum')
    pz_max = openmc.ZPlane(z0=4.0, boundary_type='vacuum')

    cell_prism = openmc.Cell(name='Sigma Pile', region=(+pxl & -pxr & +pyl & -pyr & +pz_min & -pz_max), fill=he3)

    geometry = openmc.Geometry([cell_prism])


    # ================================================================================
    # SETTINGS
    # ================================================================================
    
    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.batches = 100
    settings.particles = 10000

    # Define a Maxwellian thermal source at 0.0253 eV (room temperature)
    # Emitting from the z=0 plane, traveling monodirectionally down the +z axis
    source_dist = openmc.stats.Box([-0.5, -0.5, 0.0], [0.5, 0.5, 0.0])
    energy_dist = openmc.stats.Maxwell(theta=0.0253)
    angle_dist = openmc.stats.Monodirectional(reference_uvw=(0.0, 0.0, 1.0))

    source = openmc.IndependentSource(space=source_dist, energy=energy_dist, angle=angle_dist)
    settings.source = source


    # ================================================================================
    # TALLIES
    # ================================================================================

    tallies = openmc.Tallies()

    # Create a mesh for spatial resolution: 1x1 in XY, 80 bins in Z ( 4cm / 0.025cm)
    mesh = openmc.RegularMesh()
    mesh.dimension = [1, 1, 160]
    mesh.lower_left = [-0.5, -0.5, 0.0]
    mesh.upper_right = [0.5, 0.5, 4.0]
    filter_mesh = openmc.MeshFilter(mesh)

    # Tally 1: Total (n,p) reactions and flux (neutron population) at each depth
    tally_spatial = openmc.Tally(name='Spatial tally')
    tally_spatial.filters = [filter_mesh]
    tally_spatial.scores = ['flux', '(n,p)']
    tallies.append(tally_spatial)

    # Tally 2: Fine energy grid (n,p) reactions at each depth for histogram
    # 500 logarithmically spaced bins from 1e-5 eV to 10 eV
    energy_bins = np.logspace(-5, 1, 501)
    filter_energy = openmc.EnergyFilter(energy_bins)

    tally_energy = openmc.Tally(name='Energy-spatial tally')
    tally_energy.filters = [filter_mesh, filter_energy]
    tally_energy.scores = ['(n,p)']
    tallies.append(tally_energy)


    # ================================================================================
    # EXPORT MODEL
    # ================================================================================

    return openmc.model.Model(geometry, materials, settings, tallies)


def plot_sigma(path):
    # 1. Load the statepoint file generated after running the model
    # (Assuming settings.batches = 100)
    sp = openmc.StatePoint(f'{path}/statepoint.100.h5')

    # 2. Extract the spatial tally by its name
    tally_spatial = sp.get_tally(name='Spatial tally')

    # 3. Get the data as a Pandas DataFrame for easy manipulation
    df = tally_spatial.get_pandas_dataframe()

    # 4. Filter for the '(n,p)' score and sort by the Z-axis mesh index
    # The DataFrame typically has 'x', 'y', and 'z' columns representing the mesh bin indices
    # 4. Filter for the '(n,p)' score
    df_np = df[df['score'] == '(n,p)'].copy()

    # Dynamically find the z-coordinate column (e.g., 'mesh 1-z')
    # Dynamically find the z-coordinate column (handles both strings and tuples)
    z_col = [
        col for col in df_np.columns 
        if (isinstance(col, tuple) and col[-1] == 'z') or (isinstance(col, str) and col.endswith('-z'))
    ][0]

    # Sort using the correct column name
    df_np = df_np.sort_values(by=z_col)

    # 5. Extract the mean reaction rates
    reactions = df_np['mean'].values

    # 6. Calculate 1 - (cumulative fraction of reactions)
    # Using np.cumsum to get the running total of reactions up to each depth
    cumulative_reactions = np.cumsum(reactions)
    total_reactions = cumulative_reactions[-1]
    y_values = 1.0 - (cumulative_reactions / total_reactions)

    # 7. Calculate depth values (Z-axis)
    # The mesh has 160 bins from 0.0 to 4.0 cm. 
    # We map the cumulative sum to the upper boundary of each mesh bin.
    num_bins = 160
    z_max = 4.0
    depths = np.linspace(z_max / num_bins, z_max, num_bins)

    # Prepend the (0, 1) data point to anchor the start of the beam
    depths = np.insert(depths, 0, 0.0)
    y_values = np.insert(y_values, 0, 1.0)

    # 9. Export the data to a CSV file
    # Stack the 1D arrays as columns
    export_data = np.column_stack((depths, y_values))

    # Save to CSV with a clean header (comments='' removes the default '#' before the header)
    np.savetxt('sigma_pile_spatial_results.csv', export_data, delimiter=',', 
            header='Depth_cm,Unreacted_Fraction', comments='')

    # 8. Plot the results
    plt.figure(figsize=(8, 6))

    # A semi-log y-axis is usually best for viewing beam attenuation
    plt.plot(depths, y_values, linestyle='-', marker='.', color='b')
    # plt.yscale('log')

    plt.xlim(0, 1.0)

    plt.xlabel('Depth Z (cm)', fontsize=12)
    plt.ylabel('1 - (Cumulative Fraction of (n,p) Reactions)', fontsize=12)
    plt.title('Unreacted Neutron Beam Fraction vs. Depth in He-3', fontsize=14)
    plt.grid(True, which="both", ls="--", alpha=0.7)

    plt.tight_layout()
    plt.show()