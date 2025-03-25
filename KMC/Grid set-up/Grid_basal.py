import numpy as np
import matplotlib.pyplot as plt

# Define unit vectors
a = 1#1.303205
b = 1.732#2.25721

# Define grid size
row = 80#10
column= 96#12
# Define the unit cell
  # these points have to be copied

# Create figure
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')  # For rectangular grid

grid_points = set() # Each index, i..e, grid_points[i], is an array containing the xpos and ypos of a point [xpos,ypos]

# Draw the grid (each rectangle has width a and length b)
for i in range(row + 5):
    ax.plot([i * a, i * a], [0, row * b], color='gray', linestyle='--', linewidth=0.5)
for j in range(row + 5):
    ax.plot([0, column * a], [j * b, j * b], color='gray', linestyle='--', linewidth=0.5)
unit_cell = np.array([[1,2],[2,1],[4,1],[5,2],[7,2]]) * np.array([a, b])
#ax.scatter(unit_cell[:, 0], unit_cell[:, 1], color='red', s=50, label="Unit Cell Points")
for i in range(row):
    for j in range(column):
            translated_cell = unit_cell + np.array([i * 6 * a, j * 2 * b])  # Shift in steps of (5a, 2b)
            for point in translated_cell:
                if (point[0] <= column * a) and (point[1] <= row * b):
                    grid_points.add(tuple(point))  # Add as a tuple to the set


grid_points = np.array(sorted(grid_points, key=lambda point: (point[1], point[0])))
print(len(grid_points))
# colors = []
# for index, (x,y) in enumerate(grid_points):
#     row_index = y/b  # Determine the row number
#     col_index = x/a  # Determine the column number
#
#     if (row_index) % 2 == 0 and index %2 ==0:  # Even row: Start with red
#         colors.append('pink')
#     elif (row_index) % 2 == 0 and index %2 !=0:
#         colors.append('red')
#     else:
#         if (row_index) % 2 != 0 and index % 2 == 0:  # Even row: Start with red
#             colors.append('red')
#         elif (row_index / b) % 2 != 0 and index % 2 != 0:
#             colors.append('pink')
#
#
# # Scatter plot with alternating row colors
# ax.scatter(grid_points[:, 0], grid_points[:, 1], color=colors, s=20)

# ax.scatter(grid_points[:, 0], grid_points[:, 1], color='red', s=20)

# ax.set_aspect('equal')
# grid_points = sorted(grid_points, key=lambda point: (point[1], point[0]))
# print(len(grid_points))
# print(grid_points)
# # print(grid_points[-1])
# # af = 1
# # print(af)
# plt.xticks(np.arange(0, (column+1)*a,a))
# plt.yticks(np.arange(0, (row+1)*b,b))
#
# ax.set_xlim(0, column * a)
# ax.set_ylim(0, row * b)
# plt.show()
# ax.set_xlabel("X-axis")
# ax.set_ylabel("Y-axis")
# ax.set_title("10x10 Rectangular Grid with Properly Placed Unit Cell")

#grid_points = np.array(grid_points, dtype=object)
#
# # Print the stored unit cell points
# lateral_neighbors = []  # Dictionary to store neighbors
# for i, (x1, y1) in enumerate(grid_points):
#     for j, (x2, y2) in enumerate(grid_points):
#         if i == j:
#             continue  # Skip self
#
#         # Check if the points follow the neighbor distance rules
#         dx = abs(x1 - x2)
#         dy = abs(y1 - y2)
#
#         if (dx == 2 * a and dy == 0) or (dx == a and dy == b):
#             lateral_neighbors.append((x2, y2))
#             print(lateral_neighbors[0][0])