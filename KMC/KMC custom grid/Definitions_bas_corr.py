import numpy
import json
row = 80 #10
column = 96  #12
k_B = 1.38e-23
a_grid = 1.303205
b_grid = 2.25721

# Function to create the grid points based on the lattice constants
def create_grid():
    a_grid = 1.303205
    b_grid = 2.25721
    grid_points = set()
    unit_cell = numpy.array([[1,2],[2,1],[4,1],[5,2],[7,2]]) * numpy.array([a_grid, b_grid])
    # Generate grid points by translating the unit cell over the grid
    for i in range(row):
        for j in range(column):
            translated_cell = unit_cell + numpy.array([i * 6 * a_grid, j * 2 * b_grid])
            for point in translated_cell:
                if (point[0] <= column * a_grid) and (point[1] <= row * b_grid):
                    grid_points.add(tuple(point))  # Add as a tuple to the set

    grid_points = numpy.array(list(grid_points))
    grid_points = sorted(grid_points, key=lambda point: (point[1], point[0]))
    return grid_points

# Function to initialize arrays related to block heights, neighbors, and rate constants
def create_arrays(grid_points):
    block_height = 3*numpy.ones((len(grid_points)), dtype=int)
    block_height_ref = 3 * numpy.ones((len(grid_points)), dtype=int)
    neighbor = numpy.zeros((len(grid_points),1),  dtype=int)  # Initialize with zeros (no solid neighbors)

    k_des = numpy.zeros((len(grid_points),1), dtype=float)  # desorption and adsorption rate constant
    k_ads = numpy.zeros((len(grid_points),1), dtype=float)

    k_diff = numpy.zeros((len(grid_points),),  dtype=float)
    k_diff_hor = numpy.zeros((len(grid_points),), dtype=float)
    k_diff_ab = numpy.zeros((len(grid_points),), dtype=float)
    k_diff_below = numpy.zeros((len(grid_points),), dtype=float)

    k_x_y = numpy.zeros((len(grid_points), 1), dtype=float)  # Event frequency at each grid point

    # Initialize block heights and reference heights based on the grid pattern
    for index, (x, y) in enumerate(grid_points):
        row_index = y / b_grid  # Determine the row number

        if (row_index) % 2 == 0 and index % 2 == 0:  # Even row: Start with red
            block_height[index] = 3
            block_height_ref[index] = 3
        elif (row_index) % 2 == 0 and index % 2 != 0:
            block_height[index] = 2
            block_height_ref[index] = 2
        else:
            if (row_index) % 2 != 0 and index % 2 == 0:  # Even row: Start with red
                block_height[index] = 2
                block_height_ref[index] = 2
            elif (row_index) % 2 != 0 and index % 2 != 0:
                block_height[index] = 3
                block_height_ref[index] =3
    #k_x = numpy.zeros(len(g), dtype=float)  # frequency per row x

    return block_height, block_height_ref,neighbor, k_des, k_ads, k_diff, k_diff_hor, k_diff_ab, k_diff_below, k_x_y


# Initialise the neighbors array (index equals the # of neighbors of that block)
def initialise_neighbors(neighbor):
    neighbor[:] = 3 # PBC + every H2O molecule has 3 NB
    return neighbor

# Create rate constants (for adsorption, desorption, and diffusion) based on grid configuration
def create_rates(grid_points,neighbor,block_height, alpha, beta, k_des, k_ads, k_diff, k_diff_hor, k_diff_ab, k_diff_below, k_x_y,lateral_neighbors):
    xs = 7.4e-10
    a = 7.82e-10
    ftrial_ads_des = 1
    ftrial_diff = 1
    for x in range(len(grid_points)):
            i = neighbor[x]
            k_des[x] = ftrial_ads_des*numpy.exp(alpha/4 * (2-i) - beta/2)  # [1/s]
            k_ads[x] = ftrial_ads_des * numpy.exp(-alpha / 4 * (2-i) + beta/2 )  # [1/s]
            horizontal, above, below = lateral_neighbors[x].values()
            if(block_height[horizontal] == block_height[x]) or (block_height[horizontal] + 1 == block_height[x]):
                k_diff_hor[x] = ftrial_diff * (xs / a) ** 2 * (numpy.exp(alpha / 4 * (2 + neighbor[horizontal] - i) - beta / 2))

            if (block_height[above] == block_height[x]) or (block_height[above] + 1 == block_height[x]):
                k_diff_ab[x] = ftrial_diff * (xs / a) ** 2 * (numpy.exp(alpha / 4 * (2 + neighbor[above] - i) - beta / 2))

            if (block_height[below]== block_height[x]) or (block_height[below] + 1 == block_height[x]):
                k_diff_below[x] = ftrial_diff * (xs / a) ** 2 * (numpy.exp(alpha / 4 * (2 + neighbor[below] - i) - beta / 2))

            k_diff[x] = k_diff_hor[x]  + k_diff_ab[x] + k_diff_below[x]
            k_x_y[x] = k_des[x] + k_ads[x] + k_diff[x]
    k_tot = numpy.sum(k_x_y)



    return k_des, k_ads, k_diff, k_diff_hor, k_diff_ab, k_diff_below, k_x_y, k_tot

# Function to check the position for an event, based on random number and event frequencies
def check_position(rnd_pos, k_x_y):
    sum_k_x_y = 0
    for x,freq in enumerate(k_x_y):
        sum_k_x_y += k_x_y[x]
        if sum_k_x_y > rnd_pos:
            return x


# Function to check which event (adsorption, desorption, or diffusion) occurs at a grid point
def check_event(rnd, k_ads, k_des, k_x_y, xpos):
    if rnd < k_ads[xpos]/(k_x_y[xpos]):
        event = "adsorption"
    elif rnd < (k_ads[xpos]+k_des[xpos])/(k_x_y[xpos]):
        event = "desorption"
    else:
        event = "diffusion"
        rnd -= (k_ads[xpos] + k_des[xpos]) / (k_x_y[xpos])

    return event, rnd

# Function to update block heights and neighbors based on the chosen event (adsorption, desorption, diffusion)
def update_height_nb(rnd,xpos, k_diff_hor, k_diff_ab, k_diff_below,block_height, neighbor, event,lateral_neighbors):
    if event == "adsorption":
        block_height[xpos] += 1
        neighbor[xpos] = 0
        block_height, neighbor = adsorption_nb(xpos,block_height,neighbor,lateral_neighbors)
        diff_pos = [1]
    elif event == "desorption":
        block_height[xpos] -= 1
        neighbor[xpos] = 0
        block_height, neighbor = desorption_nb(xpos, block_height, neighbor,lateral_neighbors)
        diff_pos = [1]
    else: # random direction is chosen
        block_height, neighbor, migrated_position = diffusion_nb(rnd, xpos, k_diff_hor, k_diff_ab, k_diff_below,
                                                                 block_height, neighbor, lateral_neighbors)
        diff_pos = migrated_position

    return block_height, neighbor, diff_pos


def adsorption_nb(xpos,block_height,neighbor,lateral_neighbors):
    new_neighbors = lateral_neighbors[xpos].values()

    for posx in new_neighbors:
        if block_height[posx] == block_height[xpos]:
            neighbor[posx] += 1
            neighbor[xpos] += 1

        elif block_height[posx] == block_height[xpos]+1:
            neighbor[xpos] +=1



    return block_height, neighbor


def desorption_nb(xpos, block_height,neighbor, lateral_neighbors):
    new_neighbors = lateral_neighbors[xpos].values()

    for posx in new_neighbors:
        if block_height[posx] == block_height[xpos]:
            neighbor[xpos] += 1
            
        elif block_height[posx] == block_height[xpos]+1:
            neighbor[xpos] +=1
            neighbor[posx] -=1

    return block_height, neighbor



def diffusion_nb(rnd,xpos, k_diff_hor, k_diff_ab, k_diff_below,block_height,neighbor,lateral_neighbors):
    horizontal,above,below= lateral_neighbors[xpos].values()

    if rnd <= k_diff_hor[xpos]:   # horizontal
        if (block_height[horizontal] == block_height[xpos]) or (block_height[horizontal] + 1 == block_height[xpos]):
            neighbor[xpos] = 0
            block_height[horizontal] += 1
            block_height[xpos] -= 1
            neighbor[horizontal] = 0
            migrated_position = [xpos, horizontal]
            block_height, neighbor = neighbor_change_diffusion(lateral_neighbors,migrated_position, block_height, xpos, neighbor)

        else:
            migrated_position = [xpos]

    elif rnd <= ( k_diff_hor[xpos] + k_diff_ab[xpos] ):  # up
        if (block_height[above] == block_height[xpos]) or (block_height[above] + 1 == block_height[xpos]):
            neighbor[xpos] = 0
            block_height[above] += 1
            block_height[xpos] -= 1
            neighbor[above] = 0
            migrated_position = [xpos, above]
            block_height, neighbor = neighbor_change_diffusion(lateral_neighbors,migrated_position, block_height, xpos, neighbor)

        else:
            migrated_position = [xpos]

    elif rnd >  (k_diff_hor[xpos] + k_diff_ab[xpos]):  # down
        if (block_height[below] == block_height[xpos]) or (block_height[below] + 1 == block_height[xpos]):
            neighbor[xpos] = 0
            block_height[below ] += 1
            block_height[xpos] -= 1
            neighbor[below ] = 0
            migrated_position = [xpos,below]
            block_height, neighbor = neighbor_change_diffusion(lateral_neighbors,migrated_position, block_height, xpos, neighbor)

        else:
            migrated_position = [xpos]

    return block_height, neighbor, migrated_position

def neighbor_change_diffusion(nb_diff,migrated_position,block_height,xpos,neighbor):
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
def update_arrays(grid_points,neighbor,block_height, alpha, beta, k_des, k_ads, k_diff, k_diff_hor, k_diff_ab, k_diff_below, k_x_y, lateral_neighbors):
    xs = 7.4e-10
    a = 7.82e-10
    ftrial_ads_des = 1
    ftrial_diff = 1
    for x in range(len(grid_points)):
            i = neighbor[x]
            k_des[x] = ftrial_ads_des*numpy.exp(alpha/4 * (2-i) - beta/2)  # [1/s]
            k_ads[x] = ftrial_ads_des * numpy.exp(-alpha / 4 * (2-i) + beta/2 )  # [1/s]
            horizontal, above, below = lateral_neighbors[x].values()
            if (block_height[horizontal] == block_height[x]) or (block_height[horizontal] + 1 == block_height[x]):
                k_diff_hor[x] = ftrial_diff * (xs / a) ** 2 * (
                    numpy.exp(alpha / 4 * (2 + neighbor[horizontal] - i) - beta / 2))

            if (block_height[above] == block_height[x]) or (block_height[above] + 1 == block_height[x]):
                k_diff_ab[x] = ftrial_diff * (xs / a) ** 2 * (
                    numpy.exp(alpha / 4 * (2 + neighbor[above] - i) - beta / 2))

            if (block_height[below] == block_height[x]) or (block_height[below] + 1 == block_height[x]):
                k_diff_below[x] = ftrial_diff * (xs / a) ** 2 * (
                    numpy.exp(alpha / 4 * (2 + neighbor[below] - i) - beta / 2))

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