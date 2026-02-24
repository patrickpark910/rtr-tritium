

def rho_he3(T, P, units='psi'):
    """
    Calculate the rho_he3 of helium-3 gas given temperature and pressure.

    Parameters:
        T (float): Temperature in Kelvin
        P (float): Pressure in Pascals

    Returns:
        float: rho_he3 of helium-3 gas in kg/m^3
    """

    if units not in ['psi', 'atm', 'Pa']:
        raise ValueError("Warning. <utilities.py/rho_he3()> Unsupported units for pressure. Use 'psi', 'atm', or 'Pa'.")
    
    elif units in ['psi']:
        P = P * 6894.75729  # 1 psi = 6894.75729 Pa

    elif units in ['atm']:
        P = P * 101325.0  # 1 atm = 101325 Pa

    # Molar mass of helium-3 in kg/mol
    amu_he3 = 3.0160293e-3  # kg/mol

    # Universal gas constant in J/(mol*K)
    R = 8.314462618  # J/(mol*K)

    # Calculate the rho_he3 using the ideal gas law: ρ = (P * M) / (R * T)
    rho_kgm3 = (P * amu_he3) / (R * T)
    rho_gcm3 = rho_kgm3 / 1000.0  # Convert kg/m^3 to g/cm^3

    return rho_gcm3


if __name__ == "__main__":

    for t, p, u in [(300.0, 1.0, 'atm'), (300.0, 1000.0, 'psi')]:
        r = rho_he3(t, p, units=u)
        print(f"rho_he3 at T={t} K and P={p} {u}: {r:.6f} g/cm^3")

    """
    rho_he3 at T=300.0 K and P=  14.7 psi  (1 atm): 0.000123 g/cm^3
    rho_he3 at T=300.0 K and P=1000.0 psi (68 atm): 0.008337 g/cm^3
    """