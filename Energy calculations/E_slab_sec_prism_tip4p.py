import sys
sys.path.append('LAMMPS location')
from lammps import lammps  # Import the LAMMPS class

# Initialize the LAMMPS object
k_B = 1.38e-23  # Boltzmann constant [J/K]
Na = 6.022e23  # Avogadro [mol-1]
e_elem = 1.602e-19  # Elementary charge [C]
epsilon = (106.1 * k_B)/1000*Na  # [kJ/mol]
sigma = 3.1668  # [AA]
r_cut = 3*sigma

lmp = lammps()

lmp.command('units real')
lmp.command('boundary p p p')
lmp.command('dimension 3')
lmp.command('atom_style charge')


lmp.command(r"location_file")

lmp.command('group O type 2')  # Oxygen atoms
lmp.command('group H type 1')  # Hydrogen atoms
lmp.command('group M type 3')  # Hydrogen atoms
lmp.command('group central_molecule id 3381 3382 3383 3384 ')  # Central molecule group  #3381 3382 3383 3384 #885 886 887 888
lmp.command('group other_molecules subtract all central_molecule')  # All other molecules
lmp.command('minimize 0 0 10000 100000')
lmp.command('fix mynpt all npt temp 270 270 100 iso 1.0 1.0 1000')
lmp.command('pair_style lj/cut/coul/long 8.5')
lmp.command('pair_coeff 1 1 0.0 1.0')
lmp.command('pair_coeff 2 2 0.21084 3.1668')
lmp.command('pair_coeff 3 3 0.0 1.0')
lmp.command('pair_coeff 1 3 0.0 1.0')
lmp.command('pair_coeff 1 2 0.0 1.0')
lmp.command('pair_coeff 2 3  0.0 1.0')

lmp.command('kspace_style ewald 1.0e-4')
lmp.command('compute central_energy central_molecule group/group other_molecules')

# Output the interaction energy
lmp.command('thermo_style custom step c_central_energy')
lmp.command('thermo 1')
lmp.command('timestep 1')
lmp.command('run 0')
