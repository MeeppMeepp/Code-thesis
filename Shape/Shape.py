import itertools
import numpy


from wulffpack import (SingleCrystal,
                       Decahedron,
                       Icosahedron)
from ase.build import bulk
from ase.io import write
prim = bulk('X', crystalstructure='hcp',a=7.82,c=7.36)
surface_energies = {(0,0,0,1): 0.5*1,
                    (1,0,-1,0): 2.137,
                    (1, 1, -2,0): 3.438,
                    (2, 0, -2, 1): 1.577}
particle = SingleCrystal(surface_energies,primitive_structure=prim)
# colors = {(0,0,0,1): '#2980B9',
#          (1,0,-1,0): '#E92866',
#          (1, 1, -2,0): '#FFE82C',
#          (2, 0, -2, 1): '#7E86D2'}

colors = {(0,0,0,1): '#F8C2E4',
         (1,0,-1,0): '#FFFFB2',
         (1, 1, -2,0): '#DEF1CA',
         (2, 0, -2, 1): '#C7D4ED'}
particle.view(colors=colors)
from ase.visualize import view
view(particle.atoms)

write('icebas5.xyz', particle.atoms)

# prim = bulk('H2O', crystalstructure='hcp')
# surface_energies = {(1, 0, -1, 0): 1.1,
#                     (0, 0, 0, 1): 1.0,
#                     (1, 1, -2, 0): 1.0}
# particle = SingleCrystal(surface_energies,
#                          primitive_structure=prim,
#                          natoms=5000)
# particle.view()
# write('hcp_single_crystal.xyz', particle.atoms)
# edges = particle.edge_length
# print(edges)

# from ase.io import write
# write('atoms.xyz', particle.atoms) ## VOOR AFSTANDEN / RATIO'S