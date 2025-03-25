import numpy as np
filename = r"location_of.lmp_file"
output_filename = r"new_name_location_of.lmp_file"
# Tolerance for identifying O-H bonds
  # Scaled tolerance for bond length

# Parse the LAMMPS data
atoms = []
header = []
new_atoms = []
with open(filename, "r") as file:
    lines = file.readlines()
    for line in lines:
        if line.strip().startswith("Atoms"):
            break
        header.append(line)
    atom_data_start = lines.index(line) + 2
    atoms = [line.split() for line in lines[atom_data_start:]]

# Convert atom data to structured format
atom_data = []
for atom in atoms:
    atom_id, atom_type, x, y, z = atom
    atom_data.append({
        "id": int(atom_id),
        "type": int(atom_type),
        # "charge":float(charge),
        "position": np.array([float(x), float(y), float(z)]),
    })

new_atoms = []
for i in range(0, len(atom_data),3):
    type_O = atom_data[i]# First type 1 atom
    type_H1 = atom_data[i+1] # Second type 1 atom
    type_H2 = atom_data[i+2]    # Type 2 atom

    # Update type 2 charge
    type_O['charge'] = 0
    type_H1['charge'] = 0.5897
    type_H2['charge'] = 0.5897

    new_atoms.append(type_O)
    new_atoms.append(type_H1)
    new_atoms.append(type_H2)
    # Calculate the midpoint
    a = type_O['position']
    midpoint = (type_H2['position'] + type_H1['position']) / 2
    direction = midpoint - type_O['position']
    unit_direction = direction/np.linalg.norm(direction)
    # Calculate the new position for type 3
    type3_pos = type_O['position'] + (0.1577*unit_direction)
   
    new_atoms.append({
        "id": 0,
        "type": 3,  # Type 3 atom
        "charge": -1.1794,  # Charge for type 3 atom (example value)
        "position": type3_pos
    })

with open(output_filename, "w") as f:
    for line in header:
        f.write(line)
    f.write("\nAtoms # charge \n\n")
    for atom in new_atoms:
        f.write(
            f"{atom['id']} {atom['type']} {atom['charge']}  "
            f"{atom['position'][0]:.6f} {atom['position'][1]:.6f} {atom['position'][2]:.6f}\n"
        )
    f.write("\n")