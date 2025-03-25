import sys
sys.path.append('LAMMPS location')
from lammps import lammps  # Import the LAMMPS class
import numpy as np
from scipy.spatial import cKDTree
# Initialize the LAMMPS object
k_B = 1.38e-23  # Boltzmann constant [J/K]
Na = 6.022e23  # Avogadro [mol-1]
e_elem = 1.602e-19  # Elementary charge [C]
epsilon = (106.1 * k_B)/1000*Na  # [kJ/mol]
sigma = 3.1668  # [AA]
r_cut = 3*sigma

lmp = lammps()

# Establish base parameters and read the file (CIF to LMP via atomsk)
lmp.command('units real')
lmp.command('boundary p p p')
lmp.command('dimension 3')
lmp.command('atom_style charge')
lmp.command(r"read_data location_supercell")


central_O_id = 42729

#print(pos_slab[central_O_id-1])

UC = r"new_name_location_supercell"
data_UC = np.loadtxt(UC, skiprows=18)
id_UC = data_UC[:, 0].astype(int)         # Atom IDs (1st column)
type_UC = data_UC[:, 1].astype(int)       # Atom Types (2nd column)


# # Construct the LAMMPS command to include all IDs in other_molecules
lmp.command('group O type 2')  # Oxygen atoms
lmp.command('group H type 1')  # Hydrogen atoms
lmp.command('group M type 3')  # Hydrogen atomslaa
lmp.command(f'group central_molecule id {central_O_id} {central_O_id+1} {central_O_id+2} {central_O_id+3}')
lmp.command('group UC subtract all central_molecule')
lmp.command('minimize 1.0e-8 1.0e-8 10000 1000')
lmp.command('pair_style lj/cut/coul/long 8.5')
lmp.command('pair_modify table 0')
lmp.command('pair_coeff 1 1 0.0 0.0')
lmp.command('pair_coeff 1 2 0.0 0.0')
lmp.command('pair_coeff 1 3 0.0 0.0')
lmp.command('pair_coeff 2 2 0.21084 3.1668')
lmp.command('pair_coeff 2 3 0.0 0.0')
lmp.command('pair_coeff 3 1 0.0 0.0')
lmp.command('pair_coeff 3 2 0.0 0.0')
lmp.command('pair_coeff 3 3 0.0 0.0')

lmp.command('kspace_style ewald 1.0e-8')
lmp.command('compute energy central_molecule group/group UC pair yes kspace yes')

# Output the interaction energy
lmp.command('thermo_style custom c_energy')
lmp.command('thermo 1')
lmp.command("run 0")  # This triggers the dump without performing a simulation


lmp.command("run 0")  # This trigge