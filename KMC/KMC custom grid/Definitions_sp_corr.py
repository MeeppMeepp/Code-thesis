import numpy
import json

row = 72
column = 72
k_B = 1.38e-23
a_grid = 1
b_grid = 0.3749

# Function to create the grid points based on the lattice constants
def create_grid():
    grid_points = set()
    unit_cell = numpy.array([[1,2],[2,1],[2,2],[1,1],[1,5],[2,5],[2,6],[1,6]]) * numpy.array([a_grid, b_grid])
    # Generate grid points by translating the unit cell over the grid
    for i in range(row):
        for j in range(column):
            translated_cell = unit_cell + numpy.array([i * 2 * a_grid, j * 8 * b_grid])
            for point in translated_cell:
                if (point[0] <= column * a_grid) and (point[1] <= row * b_grid):
                    grid_points.add(tuple(point))  # Add as a tuple to the set

    grid_points = numpy.array(list(grid_points))
    grid_points = sorted(grid_points, key=lambda point: (point[1], point[0]))
    return grid_points

# Function to initialize arrays related to block heights, neighbors, and rate constants
def create_arrays(grid_points):
    block_height = 3 * numpy.ones((len(grid_points)), dtype=int)
    block_height_ref = 3 * numpy.ones((len(grid_points)), dtype=int)
    neighbor = numpy.zeros((len(grid_points), 1), dtype=int)  # Initialize with zeros (no solid neighbors)

    k_des = numpy.zeros((len(grid_points), 1), dtype=float)  # desorption and adsorption rate constant
    k_ads = numpy.zeros((len(grid_points), 1), dtype=float)

    k_diff = numpy.zeros((len(grid_points),), dtype=float)
    k_diff_hor = numpy.zeros((len(grid_points),), dtype=float)
    k_diff_ab = numpy.zeros((len(grid_points),), dtype=float)
    k_diff_below = numpy.zeros((len(grid_points),), dtype=float)

    k_x_y = numpy.zeros((len(grid_points), 1), dtype=float)  # frequency of an event in a location (x,y)

    # Initialize block heights and reference heights based on the grid pattern
    for index, (x, y) in enumerate(grid_points):
        row_index = y / b_grid  # Determine the row number
        if (index < 160) or ((160*2 <= index) and (index < 160*3)) or ((160*4 <= index) and (index < 160*5)) or ((160*6 <= index) and (index < 160*7)) or ((160*8 <= index) and (index < 160*9)) or ((160*11 <= index) and (index < 160*12)) or ((160*13 <= index) and (index < 160*14)) or ((160*15 <= index) and (index < 160*16)) or ((160*17 <= index) and (index < 160*18)) or ((160*19 <= index) and (index < 160*20)):
            if (row_index) % 2 == 0:
                if index % 2 == 0:
                    block_height[index] = 1
                    block_height_ref[index] = 1
                else:
                    block_height[index] = 2
                    block_height_ref[index] = 2

            elif (row_index) % 2 != 0:
                if index % 2 == 0:
                    block_height[index] = 3
                    block_height_ref[index] = 3
                else:
                    block_height[index] = 1
                    block_height_ref[index] = 1

        if ((160 <= index) and (index < 160*2)) or ((160*3 <= index) and (index < 160*4)) or ((160*5 <= index) and (index < 160*6)) or ((160*7 <= index) and (index < 160*8)) or ((160*9 <= index) and (index < 160*10)) or ((160*12 <= index) and (index < 160*13)) or ((160*14 <= index) and (index < 160*15)) or ((160*16 <= index) and (index < 160*17)) or ((160*18 <= index) and (index < 160*19)):
            if (row_index) % 2 == 0:
                if index % 2 == 0:
                    block_height[index] = 3
                    block_height_ref[index] = 3
                else:
                    block_height[index] = 1
                    block_height_ref[index] = 1

            elif (row_index) % 2 != 0:
                if index % 2 == 0:
                    block_height[index] = 1
                    block_height_ref[index] = 1
                else:
                    block_height[index] = 2
                    block_height_ref[index] = 2
    # k_x = numpy.zeros(len(g), dtype=float)  # frequency per row x

    return block_height, block_height_ref, neighbor, k_des, k_ads, k_diff, k_diff_hor, k_diff_ab, k_diff_below, k_x_y


# Initialise the neighbors array (index equals the # of neighbors of that block)
def initialise_neighbors(neighbor):
    neighbor[:] = 3  # PBC + every H2O molecule has 3 NB
    return neighbor

# Create rate constants (for adsorption, desorption, and diffusion) based on grid configuration
def create_rates(grid_points, neighbor, block_height, alpha, beta, k_des, k_ads, k_diff, k_diff_hor, k_diff_ab,
                 k_diff_below, k_x_y, lateral_neighbors):
    xs = 7.4e-10
    a = 7.82e-10
    ftrial_ads_des = 1
    ftrial_diff = 1
    for x in range(len(grid_points)):
        i = neighbor[x]
        k_des[x] = ftrial_ads_des * numpy.exp(alpha / 4 * (2 - i) - beta / 2)  # [1/s]
        k_ads[x] = ftrial_ads_des * numpy.exp(-alpha / 4 * (2 - i) + beta / 2)  # [1/s]
        vertical, left, right = lateral_neighbors[x].values()
        if (block_height[vertical] == block_height[x]) or (block_height[vertical] + 1 == block_height[x]):
            k_diff_hor[x] = ftrial_diff * (xs / a) ** 2 * (
                numpy.exp(alpha / 4 * (2 + neighbor[vertical] - i) - beta / 2))

        if (block_height[left] == block_height[x]) or (block_height[left] + 1 == block_height[x]):
            k_diff_ab[x] = ftrial_diff * (xs / a) ** 2 * (numpy.exp(alpha / 4 * (2 + neighbor[left] - i) - beta / 2))

        if (block_height[right] == block_height[x]) or (block_height[right] + 1 == block_height[x]):
            k_diff_below[x] = ftrial_diff * (xs / a) ** 2 * (
                numpy.exp(alpha / 4 * (2 + neighbor[right] - i) - beta / 2))

        k_diff[x] = k_diff_hor[x] + k_diff_ab[x] + k_diff_below[x]
        k_x_y[x] = k_des[x] + k_ads[x] + k_diff[x]
    k_tot = numpy.sum(k_x_y)

    return k_des, k_ads, k_diff, k_diff_hor, k_diff_ab, k_diff_below, k_x_y, k_tot

# Function to check the position for an event, based on random number and event frequencies
def check_position( rnd_pos, k_x_y):
    sum_k_x_y = 0
    for x,freq in enumerate(k_x_y):
        sum_k_x_y += k_x_y[x]
        if sum_k_x_y > rnd_pos:
            return x


# Function to check which event (adsorption, desorption, or diffusion) occurs at a grid point
def check_event(rnd, k_ads, k_des, k_x_y, xpos):
    if rnd <= k_ads[xpos] / (k_x_y[xpos]):
        event = "adsorption"
    elif rnd <= (k_ads[xpos] + k_des[xpos]) / (k_x_y[xpos]):
        event = "desorption"
    else:
        event = "diffusion"
        rnd -= (k_ads[xpos] + k_des[xpos]) / (k_x_y[xpos])

    return event,rnd

# Function to update block heights and neighbors based on the chosen event (adsorption, desorption, diffusion)
def update_height_nb(rnd, xpos, k_diff_hor, k_diff_ab, k_diff_below, block_height, neighbor, event, lateral_neighbors):
    if event == "adsorption":
        block_height[xpos] += 1
        neighbor[xpos] = 0
        block_height, neighbor = adsorption_nb(xpos, block_height, neighbor, lateral_neighbors)
        diff_pos = [1]
    elif event == "desorption":
        block_height[xpos] -= 1
        neighbor[xpos] = 0
        block_height, neighbor = desorption_nb(xpos, block_height, neighbor, lateral_neighbors)
        diff_pos = [1]
    else:  # random direction is chosen
        block_height, neighbor, migrated_position = diffusion_nb(rnd, xpos, k_diff_hor, k_diff_ab, k_diff_below, block_height, neighbor, lateral_neighbors)
        diff_pos = migrated_position

    return block_height, neighbor, diff_pos


def adsorption_nb(xpos, block_height, neighbor, lateral_neighbors):
    new_neighbors = lateral_neighbors[xpos].values()

    for posx in new_neighbors:
        if block_height[posx] == block_height[xpos]:
            neighbor[posx] += 1
            neighbor[xpos] += 1

        elif block_height[posx] == block_height[xpos] + 1:
            neighbor[xpos] += 1

    return block_height, neighbor


def desorption_nb(xpos, block_height, neighbor, lateral_neighbors):
    new_neighbors = lateral_neighbors[xpos].values()

    for posx in new_neighbors:
        if block_height[posx] == block_height[xpos]:
            neighbor[xpos] += 1

        elif block_height[posx] == block_height[xpos] + 1:
            neighbor[xpos] += 1
            neighbor[posx] -= 1

    return block_height, neighbor


def diffusion_nb(rnd, xpos, k_diff_hor, k_diff_ab, k_diff_below, block_height, neighbor, lateral_neighbors):
    vertical, left, right = lateral_neighbors[xpos].values()

    if rnd <= k_diff_hor[xpos]:  # horizontal
        if (block_height[vertical] == block_height[xpos]) or (block_height[vertical] + 1 == block_height[xpos]):
            neighbor[xpos] = 0
            block_height[vertical] += 1
            block_height[xpos] -= 1
            neighbor[vertical] = 0
            migrated_position = [xpos, vertical]
            block_height, neighbor = neighbor_change_diffusion(lateral_neighbors, migrated_position, block_height, xpos,
                                                              neighbor)

        else:
            migrated_position = [xpos]

    elif rnd <= (k_diff_hor[xpos] + k_diff_ab[xpos]):  # up
        if (block_height[left] == block_height[xpos]) or (block_height[left] + 1 == block_height[xpos]):
            neighbor[xpos] = 0
            block_height[left] += 1
            block_height[xpos] -= 1
            neighbor[left] = 0
            migrated_position = [xpos, left]
            block_height, neighbor = neighbor_change_diffusion(lateral_neighbors, migrated_position, block_height, xpos,
                                                              neighbor)

        else:
            migrated_position = [xpos]

    elif rnd <= (k_diff_hor[xpos] + k_diff_ab[xpos] + k_diff_below[xpos]):  # down
        if (block_height[right] == block_height[xpos]) or (block_height[right] + 1 == block_height[xpos]):
            neighbor[xpos] = 0
            block_height[right] += 1
            block_height[xpos] -= 1
            neighbor[right] = 0
            migrated_position = [xpos, right]
            block_height, neighbor = neighbor_change_diffusion(lateral_neighbors, migrated_position, block_height, xpos,
                                                              neighbor)
        else:
            migrated_position = [xpos]

    return block_height, neighbor, migrated_position


def neighbor_change_diffusion(nb_diff, migrated_position, block_height, xpos, neighbor):
    old_pos = migrated_position[0]
    old_neighbors = nb_diff[old_pos].values()

    new_pos = migrated_position[1]
    new_neighbors = nb_diff[new_pos].values()

    # For the neighbors of the original block on the original layer
    # Ensure neigbors lose 1 nb as the block moved down a layer
    for pos in old_neighbors:
        if block_height[pos] == block_height[xpos] + 1 or (block_height[pos] == block_height[xpos]):
            neighbor[pos] += 1
    # For the block on the lower layer, add up its neighbors again (+1 per nb) if there are more on the same level

    # For the block that diffused to a higher layer
    # Update the nb amount for the new neighbors and itself by +1 for every nb
    for pos in new_neighbors:
        if block_height[pos] == block_height[xpos]:
            neighbor[pos] += 1
            neighbor[xpos] += 1

        elif block_height[pos] == block_height[xpos] + 1:
            neighbor[xpos] += 1

    return block_height, neighbor

# Update the arrays after the event has taken place
def update_arrays(grid_points, neighbor, block_height, alpha, beta, k_des, k_ads, k_diff, k_diff_hor, k_diff_ab,
                  k_diff_below, k_x_y, lateral_neighbors):
    xs = 7.4e-10
    a = 7.82e-10
    ftrial_ads_des = 1
    ftrial_diff = 1
    for x in range(len(grid_points)):
        i = neighbor[x]
        k_des[x] = ftrial_ads_des * numpy.exp(alpha / 4 * (2 - i) - beta / 2)  # [1/s]
        k_ads[x] = ftrial_ads_des * numpy.exp(-alpha / 4 * (2 - i) + beta / 2)  # [1/s]
        vertical, left, right = lateral_neighbors[x].values()
        if (block_height[vertical] == block_height[x]) or (block_height[vertical] + 1 == block_height[x]):
            k_diff_hor[x] = ftrial_diff * (xs / a) ** 2 * (
                numpy.exp(alpha / 4 * (2 + neighbor[vertical] - i) - beta / 2))

        if (block_height[left] == block_height[x]) or (block_height[left] + 1 == block_height[x]):
            k_diff_ab[x] = ftrial_diff * (xs / a) ** 2 * (
                numpy.exp(alpha / 4 * (2 + neighbor[left] - i) - beta / 2))

        if (block_height[right] == block_height[x]) or (block_height[right] + 1 == block_height[x]):
            k_diff_below[x] = ftrial_diff * (xs / a) ** 2 * (
                numpy.exp(alpha / 4 * (2 + neighbor[right] - i) - beta / 2))

        k_diff[x] = k_diff_hor[x] + k_diff_ab[x] + k_diff_below[x]
        k_x_y[x] = k_des[x] + k_ads[x] + k_diff[x]
    k_tot = numpy.sum(k_x_y)

    return k_des, k_ads, k_diff, k_diff_hor, k_diff_ab, k_x_y, k_tot


def average_growth(block_height, block_height_ref, time,d_int):
    numpy.sum(block_height)
    numpy.sum(block_height_ref)
    growth = (d_int * 10 ** -10) * (numpy.sum(block_height) - numpy.sum(block_height_ref)) / (
                time * numpy.sum(block_height_ref))
    print(growth)
    return growth