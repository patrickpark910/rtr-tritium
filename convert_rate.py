def calculate_tritium_production(
    power_mw: float, 
    rxn_per_neutron: float = 0.012495,
    energy_per_event_mev: float = 200.0,
    neutrons_per_event: float = 2.43
) -> float:
    """
    Converts a He3(n,p)T reaction rate into grams of Tritium produced per second
    based on the reactor's thermal power.
    
    Args:
        power_mw: Reactor thermal power in Megawatts (MW).
        rxn_per_neutron: Tally result (reactions per source neutron).
        energy_per_event_mev: Energy released per fission/source event in MeV 
                              (Default: 200.0 for U-235 fission).
        neutrons_per_event: Source neutrons produced per event 
                            (Default: 2.43 for U-235 thermal fission).
                            
    Returns:
        float: Grams of Tritium produced per second.
    """
    # Fundamental Physical Constants
    JOULES_PER_MEV = 1.602176634e-13
    AVOGADRO_NUMBER = 6.02214076e23  # atoms/mol
    MOLAR_MASS_TRITIUM = 3.0160492   # g/mol for Tritium-3
    
    # 1. Convert thermal power from MW to Watts (Joules/sec)
    power_watts = power_mw * 1e6
    
    # 2. Calculate source events (fissions) per second
    energy_per_event_joules = energy_per_event_mev * JOULES_PER_MEV
    events_per_sec = power_watts / energy_per_event_joules
    
    # 3. Calculate total source neutrons produced per second
    neutrons_per_sec = events_per_sec * neutrons_per_event
    print(f"{power_mw} MW-th = {neutrons_per_sec:.2e} neutrons/second: ")
    
    # 4. Calculate Tritium atoms produced per second 
    # (1 He3(n,p)T reaction yields exactly 1 atom of Tritium)
    atoms_t_per_sec = neutrons_per_sec * rxn_per_neutron
    
    # 5. Convert atoms of Tritium to grams
    moles_t_per_sec = atoms_t_per_sec / AVOGADRO_NUMBER
    grams_t_per_sec = moles_t_per_sec * MOLAR_MASS_TRITIUM
    
    return grams_t_per_sec

# --- Example Usage ---
reactor_power = 1 # Let's assume a 1 MWth reactor
grams_per_sec = calculate_tritium_production(reactor_power)

print(f"Reactor Power: {reactor_power} MWth")
print(f"Tritium Produced: {grams_per_sec:.6e} grams/second")
# To get annual production, multiply by seconds in a year (31,536,000)
print(f"Annual Tritium: {grams_per_sec * 31536000:.4f} grams/year")