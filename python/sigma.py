import openmc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from python.utilities import *

# ================================================================================
# PARAMETERS FOR SIGMA PILE-RELATED FUNCTIONS
# ================================================================================

SIGMA_BOX_DEPTH = 2.0  # [cm]
SIGMA_BOX_BINS  = 200  # [bins] over SIGMA_BOX_DEPTH
SIGMA_E_BINS    = 20000


def build_sigma(he3_pressure=1000.0, he3_units='psi'):

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
    pz_max = openmc.ZPlane(z0=SIGMA_BOX_DEPTH, boundary_type='vacuum')

    cell_prism = openmc.Cell(name='Sigma Pile', region=(+pxl & -pxr & +pyl & -pyr & +pz_min & -pz_max), fill=he3)

    geometry = openmc.Geometry([cell_prism])


    # ================================================================================
    # SETTINGS
    # ================================================================================
    
    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.batches = 25
    settings.particles = int(4e10)

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
    mesh.dimension = [1, 1, SIGMA_BOX_BINS]
    mesh.lower_left = [-0.5, -0.5, 0.0]
    mesh.upper_right = [0.5, 0.5, SIGMA_BOX_DEPTH]
    filter_mesh = openmc.MeshFilter(mesh)

    # Tally 1: Total (n,p) reactions and flux (neutron population) at each depth
    tally_spatial = openmc.Tally(name='Spatial tally')
    tally_spatial.filters = [filter_mesh]
    tally_spatial.scores = ['flux', '(n,p)']
    tallies.append(tally_spatial)

    # Tally 2: Fine energy grid (n,p) reactions at each depth for histogram
    # 500 logarithmically spaced bins from 1e-5 eV to 10 eV
    energy_bins = np.linspace(0.0025, 0.1, SIGMA_E_BINS+1)
    filter_energy = openmc.EnergyFilter(energy_bins)

    tally_energy = openmc.Tally(name='Energy-spatial tally')
    tally_energy.filters = [filter_mesh, filter_energy]
    tally_energy.scores = ['flux', '(n,p)']
    tallies.append(tally_energy)


    # ================================================================================
    # EXPORT MODEL
    # ================================================================================

    return openmc.model.Model(geometry, materials, settings, tallies)


def extract_sigma(path):
    # 1. Load the statepoint file
    sp = openmc.StatePoint(f"{path}/statepoint.25.h5")
    tally = sp.get_tally(name='Energy-spatial tally')

    # 2. Get the DataFrame
    df = tally.get_pandas_dataframe()

    """
    tally.get_pandas_dataframe() looks like this:
            mesh 1         energy low [eV] energy high [eV] nuclide  score     mean std. dev.
                x  y    z                                                                   
    0            1  1    1        1.00e-03         1.00e-03   total  (n,p) 2.66e-06  1.14e-07
    1            1  1    1        1.00e-03         1.00e-03   total  (n,p) 2.61e-06  1.06e-07
    ...
    1599998      1  1  160        9.99e-01         9.99e-01   total  (n,p) 0.00e+00  0.00e+00
    1599999      1  1  160        9.99e-01         1.00e+00   total  (n,p) 0.00e+00  0.00e+00
    """

    # --- FLATTEN THE MULTIINDEX COLUMNS ---
    # This converts ('mesh 1', 'z') to 'mesh 1_z' and ('mean', '') to 'mean'
    df.columns = df.columns.map(lambda x: x[0] if x[1] == '' else f"{x[0]}_{x[1]}")

    # Define our new flattened string column names
    z_col     = 'mesh 1_z'
    score_col = 'score'
    elow_col  = 'energy low [eV]'
    ehigh_col = 'energy high [eV]'
    mean_col  = 'mean'
    std_col   = 'std. dev.'

    # 3. Filter for the (n,p) reaction score
    df_np = df[df[score_col] == '(n,p)'].copy()
    df_np.to_csv(f"{path}/raw_data.csv", index=False)


    # 4. Group by the Z mesh index and find the index of the max reaction rate
    idx_max = df_np.groupby(z_col)[mean_col].idxmax()
    modal_df = df_np.loc[idx_max].copy()

    # 5. Calculate the modal energy
    modal_df['mode_eV'] = (modal_df[elow_col] + modal_df[ehigh_col]) / 2.0
    modal_df['z_cm'] = (modal_df[z_col] - 0.5) * (SIGMA_BOX_DEPTH / SIGMA_BOX_BINS)

    # 6. Clean up the DataFrame for the final CSV
    final_csv_df = modal_df[[
        z_col, 'z_cm', elow_col, ehigh_col, 'mode_eV', mean_col, std_col
    ]].copy()

    # Rename the columns to final clean names for the CSV
    final_csv_df.columns = [
        'z_idx', 'z_cm', 'min_eV', 'max_eV', 
        'mode_eV', '(n,p)_rate', 'stdev'
    ]

    # Sort by depth sequentially
    final_csv_df = final_csv_df.sort_values('z_idx')

    # 7. Export to CSV
    output_filepath = f"{path}/tally_data.csv"
    final_csv_df.to_csv(output_filepath, index=False)

    print(f"Successfully extracted modal energies to {output_filepath}")


def plot_sigma(path):

    df = pd.read_csv(f"{path}/tally_data.csv")
    sub_df = df[['z_cm', '(n,p)_rate']]
    sub_df = sub_df[sub_df['z_cm'] <= 1]

    start = pd.DataFrame({'z_cm': [0.0], '(n,p)_rate': [0.0]})
    df = pd.concat([start, sub_df], ignore_index=True)

    print(df)
    
    rr_tot = df['(n,p)_rate'].sum()
    df['1-(n,p)_cum'] = 1.0 - (df['(n,p)_rate'].cumsum() / rr_tot)

    # 8. Plot the results
    plt.figure(figsize=(8, 6))

    # A semi-log y-axis is usually best for viewing beam attenuation
    plt.plot(df['z_cm'], df['1-(n,p)_cum'], linestyle='-', marker='.', color='b')
    # plt.yscale('log')

    # Add the theoretical exp(-8.9x) curve
    z_vals = np.linspace(0, 1.05, 100) # Generate smooth x-values (z_cm)
    exp_decay = np.exp(-8.9 * z_vals)
    plt.plot(z_vals, exp_decay, linestyle='--', color='r', label=r'$\exp(-8.9z)$')

    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)

    ax = plt.gca()
    ax.xaxis.set_ticks_position('both')
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    # ax.xaxis.set_minor_locator(MultipleLocator(v['x_minor']))
    
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    # ax.yaxis.set_minor_locator(MultipleLocator(0.01))

    plt.xlabel("Depth [cm]")
    plt.ylabel("Unreacted neutron fraction")
    # plt.title('Unreacted Neutron Beam Fraction vs. Depth in He-3', fontsize=14)
    plt.grid(True, which="both", ls="-", alpha=0.5)

    plt.tight_layout()
    plt.savefig(f'./_figures/fig_sigma.pdf', bbox_inches='tight', format='pdf')
    plt.savefig(f'./_figures/fig_sigma.png', bbox_inches='tight', format='png')
    plt.show()