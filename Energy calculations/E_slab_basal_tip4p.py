import sys
sys.path.append('Location LAMMPS')
from lammps import lammps  # Import the LAMMPS class

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
lmp.command('atom_style full')


lmp.command(r"read_data location_file")

lmp.command('group central_molecule molecule 1494')  # Central molecule group
lmp.command('group other_molecules subtract all central_molecule')  # All other molecules
lmp.command('minimize 1.0e-4 1.0e-4 10000 100')
lmp.command('fix rigid_molecules all rigid/nvt molecule temp 270 270 1000')
#

lmp.command('pair_style lj/cut/coul/long 12')
lmp.command('pair_coeff 1 1 0.0 1.0')
lmp.command('pair_coeff 1 2 0.0 1.0')
lmp.command('pair_coeff 1 3 0.0 1.0')
lmp.command('pair_coeff 2 3 0.0 1.0')
lmp.command('pair_coeff 3 3 0.0 1.0')
lmp.command('pair_coeff 2 2 0.21084 3.1668')


lmp.command('kspace_style ewald 1.0e-4')
lmp.command('compute central_energy central_molecule group/group other_molecules')
lmp.command('dump myDump all atom 100 dump.lammpstrj')
lmp.command('timestep 1')
# Output the interaction energy
lmp.command('thermo_style custom step c_central_energy press temp')
lmp.command('thermo 1')
lmp.command('run 0')
#
