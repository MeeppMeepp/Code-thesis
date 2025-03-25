import numpy as np

# Define rotation matrix for -45 degrees around the z-axis
theta = 0.46364# Rotation angle in radians
# rotation_matrix = np.array([
#     [np.cos(theta), -np.sin(theta), 0],
#     [np.sin(theta),  np.cos(theta), 0],
#     [0,              0,             1]
# ])

rotation_matrix = np.array([
    [np.cos(theta), 0, np.sin(theta)],
    [0,  1, 0],
    [-np.sin(theta),              0,       np.cos(theta)]
])


# Input and output file paths
input_file = r"location_pyramidal_slab"
output_file =  r"new_name_location_pyramidal_slab"

# Read the data from the file
with open(input_file, 'r') as file:
    lines = file.readlines()

# Extract original simulation box dimensions from the LAMMPS header
xlo, xhi = float(lines[5].split()[0]), float(lines[5].split()[1])
ylo, yhi = float(lines[6].split()[0]), float(lines[6].split()[1])
zlo, zhi = float(lines[7].split()[0]), float(lines[7].split()[1])

# Find the start of the atom section
atom_start = next(i for i, line in enumerate(lines) if line.startswith("Atoms"))
atom_start += 1  # Move to the first line of atom data

# Process atom data
rotated_data = []
rotated_positions = []
for line in lines[atom_start:]:
    tokens = line.split()
    if len(tokens) < 6:  # Skip lines that don't contain atom data
        rotated_data.append(line.strip())  # Add them back as-is
        continue

    # Extract x, y, z positions (columns 4, 5, 6)
    x, y, z = map(float, tokens[3:6])
    position = np.array([x, y, z])

    # Apply the rotation
    rotated_position = np.dot(rotation_matrix, position)
    rotated_positions.append(rotated_position)  # Store for bounds calculation

    # Replace the original coordinates with the rotated ones
    tokens[3:6] = [f"{val:.6f}" for val in rotated_position]
    rotated_data.append(" ".join(tokens))

# Convert list of rotated positions to a numpy array
rotated_positions = np.array(rotated_positions)

# Calculate new box dimensions
min_bounds = rotated_positions.min(axis=0)  # [min_x, min_y, min_z]
max_bounds = rotated_positions.max(axis=0)  # [max_x, max_y, max_z]

# Update the box bounds in the header
lines[5] = f"{min_bounds[0]:.6f} {max_bounds[0]:.6f} xlo xhi\n"
lines[6] = f"{min_bounds[1]:.6f} {max_bounds[1]:.6f} ylo yhi\n"
lines[7] = f"{min_bounds[2]:.6f} {max_bounds[2]:.6f} zlo zhi\n"

# Write the rotated data to a new file
with open(output_file, 'w') as file:
    # Write everything before the atom data section
    file.writelines(lines[:atom_start])
    # Write rotated atom data
    file.writelines("\n".join(rotated_data) + "\n")

print(f"Rotated data saved to {output_file}")