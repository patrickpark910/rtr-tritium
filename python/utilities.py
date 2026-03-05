import os 
"""
MATPLOTLIB SETTINGS
"""
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# Fonts
# try:
#     font_path = './Python/fonts/DIN-Regular.ttf' # DIN-Regular.ttf' # './Python/arial.ttf'
#     fm.fontManager.addfont(font_path)
#     prop = fm.FontProperties(fname=font_path)
#     plt.rcParams['font.family'] = prop.get_name()
# except:
font_path = './python/fonts/arial.ttf' # './arial.ttf'
fm.fontManager.addfont(font_path)
prop = fm.FontProperties(fname=font_path)
plt.rcParams['font.family'] = prop.get_name()

# Ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 0.5   # major x‐tick line width
plt.rcParams['ytick.major.width'] = 0.5   # major y‐tick line width
plt.rcParams['xtick.major.size'] = 6      # major x‐tick line length
plt.rcParams['ytick.major.size'] = 6      # major y‐tick line length
plt.rcParams['xtick.minor.size'] = 3      # minor x‐tick line length
plt.rcParams['ytick.minor.size'] = 3      # minor y‐tick line length
plt.rcParams['axes.linewidth']    = 0.5   # axes spines linewidth (use this instead of messing with spine.set_linewidth(1) )

# Font sizes
plt.rcParams['font.size']          = 14   # default text size for labels, legends, etc.
plt.rcParams['axes.titlesize']     = 14   # axes title
plt.rcParams['axes.labelsize']     = 16   # x- and y-axis labels
plt.rcParams['xtick.labelsize']    = 14   # x-tick labels
plt.rcParams['ytick.labelsize']    = 14   # y-tick labels
plt.rcParams['legend.fontsize']    = 12   # legend text
plt.rcParams['figure.titlesize']   = 14   # figure title


def has_statepoint(directory_path):
    """
    Check if any file starting with 'statepoint' exists in the given directory.
    
    Args:
        directory_path (str): Path to the directory to search
    
    Returns:
        bool: True if a file starting with 'statepoint' is found, False otherwise
    """
    found = False
    for filename in os.listdir(directory_path):
        if filename.startswith("statepoint"):
            found = True
    return found


def rho_he3(T, P, units='psi', output_units='g/cc'):
    """
    Calculate density of helium-3 gas given temperature and pressure.

    Parameters:
        T (float): Temperature in Kelvin
        P (float): Pressure in Pascals or psi

    Returns:
        float: density of helium-3 gas in g/cm^3 or at/cm^3 depending on output_units
    """

    # Convert pressure to Pascals
    if units not in ['psi', 'atm', 'Pa']:
        raise ValueError("Warning. <utilities.py/rho_he3()> Unsupported units for pressure. Use 'psi', 'atm', or 'Pa'.")

    elif units in ['psi']:  P = P * 6894.75729  # 1 psi = 6894.75729 Pa
    elif units in ['atm']:  P = P * 101325.0    # 1 atm = 101325 Pa
    
    if output_units not in ['g/cm^3', 'g/cc', 'at/cm^3', 'at/cc']:
        raise ValueError("Warning. <utilities.py/rho_he3()> Unsupported units for He-3 density output units. Use 'g/cm^3', 'g/cc', 'at/cm^3', 'at/cc'.")

    amu_he3 = 3.0160293e-3  # kg/mol
    R       = 8.314462618   # J/(mol*K)

    # Calculate density using the ideal gas law

    if output_units in ['g/cm^3', 'g/cc']:
        rho_kgm3 = (P * amu_he3) / (R * T)
        rho_gcm3 = rho_kgm3 / 1000.0  # Convert kg/m^3 to g/cm^3
        return rho_gcm3
    
    elif output_units in ['at/cm^3', 'at/cc']:
        n_m3 = (P * 6.022e23) / (R * T)
        n_cm3 = n_m3 / 1e6  # Convert from atoms/m^3 to atoms/cm^3
        return n_cm3



if __name__ == "__main__":

    for t, p, u in [(300.0, 1.0, 'atm'), (300.0, 1000.0, 'psi')]:
        r = rho_he3(t, p, units=u)
        print(f"rho_he3 at T={t} K and P={p} {u}: {r:.6f} g/cm^3")

    """
    rho_he3 at T=300.0 K and P=  14.7 psi  (1 atm): 0.000123 g/cm^3
    rho_he3 at T=300.0 K and P=1000.0 psi (68 atm): 0.008337 g/cm^3
    """