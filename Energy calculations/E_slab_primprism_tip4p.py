import sys
sys.path.append('location LAMMPS')
from lammps import lammps  # Import the LAMMPS class
from ase.io import read
from ase.visualize import view
# Initialize the LAMMPS object
k_B = 1.38e-23  # Boltzmann constant [J/K]
Na = 6.022e23  # Avogadro [mol-1]
e_elem = 1.602e-19  # Elementary charge [C]
sigma = 3.1668  # [AA]
r_cut = 3*sigma

lmp = lammps()
# atoms = read(r"C:\Users\miepv\Documents\Master\Master\Thesis\Crystalgrower\CIF files\orthogonal_ice.cif")
# view(atoms)  # Interactive visualization

# Establish base parameters and read the file (CIF to LMP via atomsk)
lmp.command('units real')
lmp.command('boundary p p p')
lmp.command('dimension 3')
lmp.command('atom_style charge')


lmp.command(r"read_data location_file")


lmp.command('group O type 2')  # Oxygen atoms
lmp.command('group H type 1')  # Hydrogen atoms
lmp.command('group M type 3')  # Hydrogen atoms
lmp.command('group central_molecule id 7740 7739 7738 7737')  # Central molecule group
# lmp.command('group central_molecule id 2041 2042 2043')
lmp.command('group other_molecules subtract all central_molecule')  # All other molecules
lmp.command('minimize 1.0e-8 1.0e-8 1000 10')
lmp.command('pair_style lj/cut/coul/long 8.5')
lmp.command('pair_coeff 1 1 0.0 0.0')
lmp.command('pair_coeff 1 2 0.0 0.0')
lmp.command('pair_coeff 1 3 0.0 0.0')
lmp.command('pair_coeff 2 1 0.0 0.0')
lmp.command('pair_coeff 2 2 0.21084 3.1668')
lmp.command('pair_coeff 2 3 0.0 0.0')
lmp.command('pair_coeff 3 1 0.0 0.0')
lmp.command('pair_coeff 3 2 0.0 0.0')
lmp.command('pair_coeff 3 3 0.0 0.0')

lmp.command('kspace_style ewald 1.0e-8')
lmp.command('compute central_energy central_molecule group/group other_molecules pair yes kspace yes')


lmp.command('timestep 1')
# Output the interaction energy
lmp.command('thermo_style custom c_central_energy')
lmp.command('thermo 1')
lmp.command('run 0')
