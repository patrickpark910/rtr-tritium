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
SIGMA_PARTICLES = int(4e10)
SIGMA_BATCHES   = 25
SIGMA_TEMP_K    = 300.0


def build_sigma(he3_pressure=1000.0, he3_units='psi'):

    # ================================================================================
    # MATERIALS
    # ================================================================================

    he3_gcm3 = rho_he3(T=SIGMA_TEMP_K, P=he3_pressure, units=he3_units)

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
    settings.batches = SIGMA_BATCHES
    settings.particles = SIGMA_PARTICLES

    """
    Emitting from the z=0 plane, traveling monodirectionally down the +z axis.
    Define a Maxwellian thermal source with most probable velocity corresponding to 0.0253 eV.
    """
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


def extract_sigma(path, he3_pressure=1000.0, he3_units='psi'):
    
    sp = openmc.StatePoint(f"{path}/statepoint.{SIGMA_BATCHES}.h5")
    tally1 = sp.get_tally(name='Energy-spatial tally')
    df1 = tally1.get_pandas_dataframe()
    
    """
    >>> print(df1.head(4))
      mesh 1       energy low [eV] energy high [eV] nuclide  score          mean     std. dev.
           x  y  z                                                                            
    0      1  1  1        0.002500         0.002505   total   flux  5.443157e-07  6.302758e-09
    1      1  1  1        0.002500         0.002505   total  (n,p)  1.532470e-05  1.774481e-07
    2      1  1  1        0.002505         0.002510   total   flux  5.427214e-07  8.261216e-09
    3      1  1  1        0.002505         0.002510   total  (n,p)  1.526500e-05  2.323632e-07
    """

    # Flatten the multi-index columns into single string column names for easier handling
    # This converts ('mesh 1', 'z') to 'mesh 1_z' and ('mean', '') to 'mean'
    df1.columns = df1.columns.map(lambda x: x[0] if x[1] == '' else f"{x[0]}_{x[1]}")
    df1['z mid [cm]'] = (df1['mesh 1_z'] - 0.5) * SIGMA_BOX_DEPTH / SIGMA_BOX_BINS

    # Filter for the (n,p) reaction score
    df_np    = df1[df1['score'] == '(n,p)'].copy()
    df_fluxV = df1[df1['score'] == 'flux'].copy()

    """ 
    >>> print(df_np.head(4))
       mesh 1_x  mesh 1_y  mesh 1_z  energy low [eV]  energy high [eV] nuclide  score      mean     std. dev.  z mid [cm]
    1         1         1         1         0.002500          0.002505   total  (n,p)  0.000015  1.774481e-07       0.005
    3         1         1         1         0.002505          0.002510   total  (n,p)  0.000015  2.323632e-07       0.005
    5         1         1         1         0.002510          0.002515   total  (n,p)  0.000015  2.254258e-07       0.005
    7         1         1         1         0.002515          0.002520   total  (n,p)  0.000015  2.874775e-07       0.005

    >>> print(df_fluxV.head(4))
       mesh 1_x  mesh 1_y  mesh 1_z  energy low [eV]  energy high [eV] nuclide score          mean     std. dev.  z mid [cm]
    0         1         1         1         0.002500          0.002505   total  flux  5.443157e-07  6.302758e-09       0.005
    2         1         1         1         0.002505          0.002510   total  flux  5.427214e-07  8.261216e-09       0.005
    4         1         1         1         0.002510          0.002515   total  flux  5.421151e-07  8.022359e-09       0.005
    6         1         1         1         0.002515          0.002520   total  flux  5.382456e-07  1.024056e-08       0.005
    """

    # Group by the Z mesh index and find the index of the max reaction rate
    df_np_modes_idx = df_np.groupby("mesh 1_z")['mean'].idxmax()
    df_fluxV_modes_idx = df_fluxV.groupby("mesh 1_z")['mean'].idxmax()

    df_np_modes     = df_np.loc[df_np_modes_idx].copy()
    df_fluxV_modes  = df_fluxV.loc[df_fluxV_modes_idx].copy()

    df_np_modes['energy mid [eV]']    = df_np_modes['energy low [eV]'] + (df_np_modes['energy high [eV]'] - df_np_modes['energy low [eV]']) / 2.0
    df_fluxV_modes['energy mid [eV]'] = df_fluxV_modes['energy low [eV]'] + (df_fluxV_modes['energy high [eV]'] - df_fluxV_modes['energy low [eV]']) / 2.0

    """ 
    >>> print(df_np_modes.head(4))
            mesh 1_x  mesh 1_y  mesh 1_z  energy low [eV]  energy high [eV] nuclide  score      mean     std. dev.  z mid [cm]  energy mid [eV]
    119            1         1         1         0.002788          0.002792   total  (n,p)  0.000016  2.145889e-07       0.005         0.002790
    40515          1         1         2         0.003753          0.003758   total  (n,p)  0.000012  1.514260e-07       0.015         0.003755
    81919          1         1         3         0.007175          0.007180   total  (n,p)  0.000010  1.440632e-07       0.025         0.007178
    122425         1         1         4         0.008408          0.008413   total  (n,p)  0.000008  1.389684e-07       0.035         0.008411
    """

    # =========================================================================
    # Now we extract the total (n,p) and flux in each bin
    # =========================================================================

    tally2 = sp.get_tally(name='Spatial tally')
    df2 = tally2.get_pandas_dataframe()
    df2.columns = df2.columns.map(lambda x: x[0] if x[1] == '' else f"{x[0]}_{x[1]}")
    df2['z mid [cm]'] = (df2['mesh 1_z'] - 0.5) * SIGMA_BOX_DEPTH / SIGMA_BOX_BINS

    """
    >>> print(df2.head(4))
       mesh 1_x  mesh 1_y  mesh 1_z nuclide  score      mean     std. dev.  z mid [cm]
    0         1         1         1   total   flux  0.009527  1.446081e-07       0.005
    1         1         1         1   total  (n,p)  0.092982  4.800248e-06       0.005
    2         1         1         2   total   flux  0.008658  2.754392e-07       0.015
    3         1         1         2   total  (n,p)  0.081297  4.892962e-06       0.015
    """

    df_totals = df2.pivot(index=['mesh 1_z', 'z mid [cm]'], 
                          columns='score', 
                          values=['mean', 'std. dev.'])

    # Flatten the multi-level columns by changing (mean, flux) to 'flux' and (std. dev., flux) to 'flux sd'
    df_totals.columns = [f'{score} sd' if val == 'std. dev.' else score 
                    for val, score in df_totals.columns]

    # Reset the index to make mesh 1_z and z mid [cm] regular columns again
    df_totals = df_totals.reset_index()
    df_totals = df_totals.rename(columns={"(n,p)":"(n,p) [rxns]", 'flux': 'flux*V [trks-cm]', 'flux sd': "flux*V sd"})

    # Calculate macroscopic cross section of (n,p) reactions in each bin: Sigma = R / (flux * V)
    df_totals["(n,p) macro-xs [1/cm]"] = df_totals["(n,p) [rxns]"] / df_totals["flux*V [trks-cm]"]
    
    # Calculate microscopic cross section of (n,p) reactions in each bin: sigma = Sigma / N, where N is the number density of He-3 atoms
    N_he3 = rho_he3(T=SIGMA_TEMP_K, P=he3_pressure, units=he3_units, output_units='at/cm^3')
    df_totals["(n,p) micro-xs [b]"] = df_totals["(n,p) macro-xs [1/cm]"] / N_he3 * 1e24

    """ 
    >>> print(df_totals.head(4))
       mesh 1_z  z mid [cm]  (n,p) [rxns]  flux*V [trks-cm]      (n,p) sd     flux*V sd  (n,p) macro-xs [1/cm]  (n,p) micro-xs [b]
    0         1       0.005      0.092973          0.009527  5.612909e-08  2.093049e-09               9.759244         5862.891812
    1         2       0.015      0.081291          0.008658  6.548187e-08  3.717248e-09               9.389394         5640.703272
    2         3       0.025      0.071843          0.007894  7.729814e-08  4.600615e-09               9.101407         5467.694801
    3         4       0.035      0.063941          0.007216  7.308039e-08  5.001062e-09               8.861406         5323.513257
    """

    # Merge the mode energy data with the total (n,p) and flux data based on the Z mesh index
    df_totals = df_totals.merge( df_np_modes[['mesh 1_z', 'energy mid [eV]']], on='mesh 1_z', how='left')
    df_totals = df_totals.rename(columns={'energy mid [eV]': '(n,p) mode [eV]'})

    df_totals = df_totals.merge( df_fluxV_modes[['mesh 1_z', 'energy mid [eV]']], on='mesh 1_z', how='left')
    df_totals = df_totals.rename(columns={'energy mid [eV]': 'flux mode [eV]'})

    print(df_totals.head(4))

    """ 
    >>> print(df_totals.head(4))
       mesh 1_z  z mid [cm]  (n,p) [rxns]  flux*V [trks-cm]      (n,p) sd     flux*V sd  (n,p) macro-xs [1/cm]  (n,p) micro-xs [b]  (n,p) mode [eV]  flux mode [eV]
    0         1       0.005      0.092973          0.009527  5.612909e-08  2.093049e-09               9.759244         5862.891812         0.002502        0.013525
    1         2       0.015      0.081291          0.008658  6.548187e-08  3.717248e-09               9.389394         5640.703272         0.004194        0.014763
    2         3       0.025      0.071843          0.007894  7.729814e-08  4.600615e-09               9.101407         5467.694801         0.005720        0.016299
    3         4       0.035      0.063941          0.007216  7.308039e-08  5.001062e-09               8.861406         5323.513257         0.007207        0.017283
    """

    # =========================================================================
    # Save all our hard work :3 >.<
    # =========================================================================

    try:
        df_np_modes.to_csv(f"{path}/tally_np_modes.csv", index=False)
        df_fluxV_modes.to_csv(f"{path}/tally_fluxV_modes.csv", index=False)
        df_totals.to_csv(f"{path}/tally_totals.csv", index=False)
        print(f"Comment. <sigma.py/extract_sigma()> Successfully extracted tally data to CSVs in: {path}")
        return
    except:
        print(f"Warning. <sigma.py/extract_sigma()> Failed to export tally data to: {path}")
        return


def plot_sigma(path):

    try:
        df = pd.read_csv(f"{path}/tally_np_data.csv")
    except FileNotFoundError:
        print(f"Fatal. <sigma.py/plot_sigma()> Could not find: {path}/tally_np_data.csv")
        return
    
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