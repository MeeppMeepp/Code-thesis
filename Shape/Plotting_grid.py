import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Define file path
excelfile = r"location_grid"

# Read the data from the Excel file
data = pd.read_excel(excelfile)

# Extract X, Y, and Height values from the Excel file
time = data.iloc[:, 0]
xpos1, ypos1, height1 = data.iloc[:, 1], data.iloc[:, 2], data.iloc[:, 3]
xpos2, ypos2, height2 = data.iloc[:, 4], data.iloc[:, 5], data.iloc[:, 6]

# Initialize parameters
maxX, maxY = 80, 80  # Grid dimensions
numFrames = len(time)  # Number of frames

# Create a folder to save images
folderPath = os.path.join(os.getcwd(), 'frames_height')
os.makedirs(folderPath, exist_ok=True)

# Initialize grid with ones (1 represents empty cells)
grid = np.full((maxY, maxX), 3)

# Define a colormap
cmap = plt.cm.turbo.reversed()
boundaries = [1, 2, 3, 4, 5, 6]
norm = mcolors.BoundaryNorm(boundaries, cmap.N)
# Compute tick positions as the midpoints of each bin
tick_positions = [(boundaries[i] + boundaries[i + 1]) / 2 for i in range(len(boundaries) - 1)]
tick_labels = [str(i) for i in range(1, len(boundaries))]

# --- SAVE INITIAL EMPTY GRID ---
plt.figure(figsize=(8, 8))
img = plt.imshow(grid, cmap=cmap, norm=norm)

# Create colorbar
cbar = plt.colorbar(img, ticks=tick_positions, label='Height')
cbar.set_ticklabels(tick_labels)

for i in range(maxX):
    plt.axvline(i, color='k', linewidth=0.5)
for j in range(maxY):
    plt.axhline(j, color='k', linewidth=0.5)


# Set axis limits and ticks
plt.xticks(np.arange(0, maxX, 10))
plt.yticks(np.arange(0, maxY, 10))
plt.xlim(0, maxX - 1)
plt.ylim(maxY - 1, 0)  # Invert Y-axis for correct visualization



# Labels and title
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('t = 0.0')
plt.grid(False)

# Save the figure
initialFileName = os.path.join(folderPath, 'initial_grid.pdf')
plt.savefig(initialFileName)
plt.close()


# Process each timestep
# for t in range(numFrames):
#     x1, y1, h1 = (xpos1[t]), (ypos1[t]), height1[t]
#     x2, y2, h2 = (xpos2[t]), (ypos2[t]), height2[t]
#
#     # Update grid values
#     grid[y1, x1] = h1  # Note: Matplotlib follows (row, col) indexing
#     grid[y2, x2] = h2
#
#     # Save image only at specific timesteps
#     if t % 100 == 0 or t == 0 or t == numFrames - 1:
#         plt.figure(figsize=(8, 8))
#
#         # Use extent to align the grid with the image cells
#         plt.imshow(grid, cmap=cmap, norm=norm, origin = 'upper', extent=[0, maxX, maxY, 0], interpolation= 'None')
#
#         # Add colorbar
#         plt.colorbar(ticks=range(1, 6), label='Height')
#
#         # Set axis limits and ticks
#         plt.xticks(np.arange(0, maxX, 10))
#         plt.yticks(np.arange(0, maxY, 10))
#         plt.xlim(0, maxX)
#         plt.ylim(maxY, 0)  # Invert Y-axis for correct visualization
#
#         # Draw grid lines
#         for i in range(maxX):
#             plt.axvline(i, color='k', linewidth=0.5)
#         for j in range(maxY):
#             plt.axhline(j, color='k', linewidth=0.5)
#
#         # Labels and title
#         plt.xlabel('X Position')
#         plt.ylabel('Y Position')
#         plt.title(f't = {time[t]}')
#         plt.grid(False)
#
#         # Save the figure
#         frameFileName = os.path.join(folderPath, f'frame_{t:03d}.pdf')
#         plt.savefig(frameFileName)
#         plt.close()