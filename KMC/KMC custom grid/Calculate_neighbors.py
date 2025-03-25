import numpy
import json


#-----------------------------------BASAL--------------------------------------------------------

# row = 80 #10
# column = 96 #12
# k_B = 1.38e-23
#
# a_grid = 1.303205
# b_grid = 2.25721
# grid_points = set()
# unit_cell = numpy.array([[1,2],[2,1],[4,1],[5,2],[7,2]]) * numpy.array([a_grid, b_grid])
# for i in range(row):
#     for j in range(column):
#         translated_cell = unit_cell + numpy.array([i * 6 * a_grid, j * 2 * b_grid])
#         for point in translated_cell:
#             if (point[0] <= column * a_grid) and (point[1] <= row * b_grid):
#                 grid_points.add(tuple(point))  # Add as a tuple to the set
#
# grid_points = numpy.array(list(grid_points))
# grid_points = sorted(grid_points, key=lambda point: (point[1], point[0]))
#
#
# a_grid = 1.303205
# b_grid = 2.25721
# lateral_neighbors = {}
# for i,pos_ref in enumerate(grid_points):
#     xpos_point = pos_ref[0]
#     ypos_point = pos_ref[1]
#     neighbors = {'horizontal': None, 'above': None, 'below': None}
#     for j,pos in enumerate(grid_points):
#         if i == j:  # Skip self
#             continue
#         dx_new = 0
#         dy_new = 0
#         dx = (xpos_point - pos[0])
#         dy = (ypos_point - pos[1])
#
#         if dx < -(column*a_grid)/2:
#             dx_new =  dx + (column*a_grid)
#         elif dx >= (column*a_grid)/2:
#             dx_new = dx - (column*a_grid)
#         else:
#             dx = dx
#
#         if dy < -(row *b_grid) / 2:
#             dy_new = dy + (row *b_grid)
#
#         elif dy >= (row *b_grid) / 2:
#             dy_new = dy - (row *b_grid)
#         else:
#             dy = dy
#
# Check if the distance equals that of a neighboring molecule
#         if (round(abs(dx),3) == round(2*a_grid,3) and dy == 0) or (round(abs(dx_new),3) == round(2*a_grid,3) and dy == 0):
#             neighbors['horizontal'] = j
#
#         elif ((round(abs(dx),3) == round(a_grid,3)) and (round(abs(dy_new),3) == round(b_grid,3))) or ((round(abs(dx),3) == round(a_grid,3)) and (round(abs(dy),3) == round(b_grid,3))):
#             if (round(dy,3) >= row*b_grid/2):
#                 neighbors['above'] = j
#             elif (round(dy,3) < -row*b_grid/2):
#                 neighbors['below'] = j
#             elif round(b_grid + dy,2) == 0:
#                 neighbors['above'] = j
#             else:
#                 neighbors['below'] = j
#
#
#
#         if all(val is not None for val in neighbors.values()):
#             break
#
#     lateral_neighbors[i] = neighbors
#
# # print(lateral_neighbors)
#
# # Save to a JSON file
# with open("neighbors.json", "w") as f:
#     json.dump(lateral_neighbors, f)
#
# print("Neighbors computed and saved to neighbors.json")
#




#-----------------------------------PRIMARY PRISM--------------------------------------------------------

# row = 144 #16
# column = 108 #12
# k_B = 1.38e-23
#
# a_grid = 1
# b_grid = 0.7508
# grid_points = set()
# unit_cell = numpy.array([[1, 2], [2, 1], [1, 5], [2, 6], [2, 9]]) * numpy.array([a_grid, b_grid])
# for i in range(row):
#     for j in range(column):
#         translated_cell = unit_cell + numpy.array([i * 3 * a_grid, j * 8 * b_grid])
#         for point in translated_cell:
#             if (point[0] <= column * a_grid) and (point[1] <= row * b_grid):
#                 grid_points.add(tuple(point))  # Add as a tuple to the set
#
# grid_points = numpy.array(list(grid_points))
# grid_points = sorted(grid_points, key=lambda point: (point[1], point[0]))
#
#
# a_grid = 1
# b_grid = 0.7508
# lateral_neighbors = {}
# for i,pos_ref in enumerate(grid_points):
#     xpos_point = pos_ref[0]
#     ypos_point = pos_ref[1]
#     neighbors = {'vertical': None, 'short': None, 'long': None}
#     for j,pos in enumerate(grid_points):
#         if i == j:  # Skip self
#             continue
#
#         dx_new = 0
#         dy_new = 0
#         dx = (xpos_point - pos[0])
#         dy = (ypos_point - pos[1])
#
#         if dx < -(column*a_grid)/2:
#             dx_new =  dx + (column*a_grid)
#         elif dx >= (column*a_grid)/2:
#             dx_new = dx - (column*a_grid)
#         else:
#             dx = dx
#
#         if dy < -(row *b_grid) / 2:
#             dy_new = dy + (row *b_grid)
#
#         elif dy >= (row *b_grid) / 2:
#             dy_new = dy - (row *b_grid)
#         else:
#             dy = dy
#
#
#         if (round(abs(dy),2) == round(3*b_grid,2) and dx == 0) or (round(abs(dy_new),2) == round(3*b_grid,2) and dx == 0):
#             neighbors['vertical'] = j
#
#         elif ((round(abs(dx),3) == round(a_grid,3)) and (round(abs(dy_new),3) == round(b_grid,3))) or ((round(abs(dx),3) == round(a_grid,3)) and (round(abs(dy),3) == round(b_grid,3)))\
#                 or ((round(abs(dx_new),3) == round(a_grid,3)) and (round(abs(dy_new),3) == round(b_grid,3))) or ((round(abs(dx_new),3) == round(a_grid,3)) and (round(abs(dy),3) == round(b_grid,3))):
#             neighbors['short'] = j
#
#         elif ((round(abs(dx), 3) == round(2*a_grid, 3)) and (round(abs(dy_new), 3) == round(b_grid, 3))) or (
#                     (round(abs(dx), 3) == round(2*a_grid, 3)) and (round(abs(dy), 3) == round(b_grid, 3))) or((round(abs(dx_new), 3) == round(2*a_grid, 3)) and (round(abs(dy_new), 3) == round(b_grid, 3))) or (
#                          (round(abs(dx_new), 3) == round(2*a_grid, 3)) and (round(abs(dy), 3) == round(b_grid, 3))):
#             neighbors['long'] = j
#
#         if all(val is not None for val in neighbors.values()):
#             break
#
#     lateral_neighbors[i] = neighbors
#
# # print(lateral_neighbors)
#
# # Save to a JSON file
# with open("neighborspp.json", "w") as f:
#     json.dump(lateral_neighbors, f)
#
# print("Neighbors computed and saved to neighbors.json")


#-----------------------------------SECONDARY PRISM--------------------------------------------------------
# row = 72
# column = 72
# k_B = 1.38e-23
#
# a_grid = 1
# b_grid = 0.3749
# grid_points = set()
# unit_cell = numpy.array([[1,2],[2,1],[2,2],[1,1],[1,5],[2,5],[2,6],[1,6]]) * numpy.array([a_grid, b_grid])
# for i in range(row):
#     for j in range(column):
#         translated_cell = unit_cell + numpy.array([i * 2 * a_grid, j * 8 * b_grid])
#         for point in translated_cell:
#             if (point[0] <= column * a_grid) and (point[1] <= row * b_grid):
#                 grid_points.add(tuple(point))  # Add as a tuple to the set
#
# grid_points = numpy.array(list(grid_points))
# grid_points = sorted(grid_points, key=lambda point: (point[1], point[0]))
#
# lateral_neighbors = {}
# for i,pos_ref in enumerate(grid_points):
#     xpos_point = pos_ref[0]
#     ypos_point = pos_ref[1]
#     neighbors = {'vertical': None, 'left': None, 'right': None}
#     for j,pos in enumerate(grid_points):
#         if i == j:  # Skip self
#             continue
#
#         dx_new = 0
#         dy_new = 0
#         dx = (xpos_point - pos[0])
#         dy = (ypos_point - pos[1])
#
#         if dx < -(column*a_grid)/2:
#             dx_new = dx + (column*a_grid)
#         elif dx >= (column*a_grid)/2:
#             dx_new = dx - (column*a_grid)
#         else:
#             dx = dx
#
#         if dy < -(row *b_grid) / 2:
#             dy_new = dy + (row *b_grid)
#
#         elif dy >= (row *b_grid) / 2:
#             dy_new = dy - (row *b_grid)
#         else:
#             dy = dy
#
#
#         if (round(abs(dy),2) == round(b_grid,2) and dx == 0) or (round(abs(dy_new),2) == round(3*b_grid,2) and dx == 0):
#             neighbors['vertical'] = j
#
#
#         elif ((round(abs(dx),2) == round(a_grid,2)) and dy==0) or ((round(abs(dx_new),2) == round(a_grid,2)) and dy==0):
#             if (round(dx,3) >= column*a_grid/2) or (round(dx,3) == -a_grid):
#                 neighbors['right'] = j
#             elif (round(dx,3) < -column*a_grid/2) or (round(dx,3) == a_grid):
#                 neighbors['left'] = j
#
#
#         if all(val is not None for val in neighbors.values()):
#             break
#
#     lateral_neighbors[i] = neighbors
#
# # print(lateral_neighbors)
#
# # Save to a JSON file
# with open("neighborssp.json", "w") as f:
#     json.dump(lateral_neighbors, f)
#
# print("Neighbors computed and saved to neighbors.json")
#
