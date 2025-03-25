import numpy as np
import matplotlib.pyplot as plt

# Define unit vectors
a = 1#0.1667
b = 0.3749#0.0625

# Define grid size
row = 72#8
column= 72#8
# Define the unit cell
  # these points have to be copied

# Create figure
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_aspect('equal')  # For rectangular grid

grid_points = set() # Each index, i..e, grid_points[i], is an array containing the xpos and ypos of a point [xpos,ypos]
unit_cell = np.array([[1,2],[2,1],[2,2],[1,1],[1,5],[2,5],[2,6],[1,6]]) * np.array([a, b])
# Draw the grid (each rectangle has width a and length b)
for i in range(row + 1):
    ax.plot([i * a, i * a], [0, row * b], color='gray', linestyle='--', linewidth=0.5)
for j in range(row + 1):
    ax.plot([0, row * a], [j * b, j * b], color='gray', linestyle='--', linewidth=0.5)
#ax.scatter(unit_cell[:, 0], unit_cell[:, 1], color='red', s=50, label="Unit Cell Points")
for i in range(row):
    for j in range(column):
            translated_cell = unit_cell + np.array([i * 2 * a, j *8 * b])  # Shift in steps of (5a, 2b)
            for point in translated_cell:
                if (point[0] <= column * a) and (point[1] <= row * b):
                    grid_points.add(tuple(point))  # Add as a tuple to the set

grid_points = np.array(sorted(grid_points, key=lambda point: (point[1], point[0])))
print(len(grid_points))

# colors = []
#
# for index, (x,y) in enumerate(grid_points):
#     row_index = y/b  # Determine the row number
#     col_index = x/a  # Determine the column number
#     if index < 16:
#         if (row_index) % 2 == 0:
#             if index % 2 == 0:
#                 colors.append('pink')
#             else:
#                 colors.append('red')
#
#         elif (row_index) % 2 != 0:
#             if index % 2 == 0:
#                 colors.append('maroon')
#             else:
#                 colors.append('pink')
#
#     if (16 <= index) and (index < 32):
#         if (row_index) % 2 == 0:
#             if index % 2 == 0:
#                 colors.append('maroon')
#             else:
#                 colors.append('pink')
#
#         elif (row_index) % 2 != 0:
#             if index % 2 == 0:
#                 colors.append('pink')
#             else:
#                 colors.append('red')
#
#
#
# ax.scatter(grid_points[:, 0], grid_points[:, 1], color=colors, s=20)
#
# plt.xticks(np.arange(0, (column+1)*a,a))
# plt.yticks(np.arange(0, (row+1)*b,b))
#
# ax.set_xlim(0, column * a)
# ax.set_ylim(0, row * b)
# plt.show()
# ax.set_xlabel("X-axis")
# ax.set_ylabel("Y-axis")
# ax.set_title("10x10 Rectangular Grid with Properly Placed Unit Cell")
