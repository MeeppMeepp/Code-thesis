# Code-thesis
This repository contains the files used within my thesis for the devlopment of the KMC algorithm for crystal growth. 
There are 4 different folders containing:

- Crystal structures
  * .cif files of the hexagonal and orthogonal unit cells of Ice Ih
  * Python code to apply the TIP4P model to the structures, and a code that rotates the pyramidal slab to be upright
  * Slab structures with TIP4P applied
  * Proton ordered slab structures (.lmp) for the basal, primary prism, secondary prism, and pyramidal plane


- Energy calculations
  * The LAMMPS code to perform energy calculations on the slabs

- The KMC algorithm
  * The grid set-ups
  * The KMC code per plane and per grid

- The Wulff algorithm
