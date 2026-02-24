import openmc.data

# Load the HDF5 file for a major nuclide
filepath = 'U235.h5'
u235 = openmc.data.IncidentNeutron.from_hdf5(filepath)

# In the ENDF-6 format, MT=301 is total heating and MT=318 is fission heating
has_total_heating = 301 in u235.reactions
has_fission_heating = 318 in u235.reactions

print(f"Total heating (MT=301) present: {has_total_heating}")
print(f"Fission heating (MT=318) present: {has_fission_heating}")