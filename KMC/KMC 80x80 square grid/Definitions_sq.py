import numpy
# Constants
row = 80  # Number of rows in the grid
column = 80  # Number of columns in the grid
k_B = 1.38e-23  # Boltzmann constant (J/K)

# Function to create arrays for block heights, neighbors, rate constants, and frequencies
def create_arrays():
    # Initialize block height arrays, neighbor count, rate constants for desorption and adsorption
    block_height = 3 * numpy.ones((row, column), dtype=int)
    block_height_ref = 3 * numpy.ones((row, column), dtype=int)
    neighbor = numpy.zeros((row, column), dtype=int)  # Initialize with zeros (no solid neighbors)
    k_des = numpy.zeros((row, column), dtype=float)  # Desorption rate constant
    k_ads = numpy.zeros((row, column), dtype=float)  # Adsorption rate constant
    k_diff = numpy.zeros((row, column), dtype=object)  # Diffusion rate constant (array of lists)
    k_x_y = numpy.zeros((row, column), dtype=float)  # Frequency of an event in a location (x, y)
    k_x = numpy.zeros(row, dtype=float)  # Frequency per row

    return block_height, block_height_ref, neighbor, k_des, k_ads, k_diff, k_x_y, k_x

# Function to initialize the neighbors array
def initialise_neighbors(neighbor):
    # Initialize all grid locations with 4 neighbors (Periodic Boundary Condition)
    neighbor[:, :] = 4
    return neighbor

# Function to update arrays for desorption, adsorption, and diffusion
def update_arrays(neighbor, block_height, alpha, beta, k_des, k_ads, k_diff, k_x_y, k_x):
    ftrial_ads_des = 1  # Trial factor for adsorption/desorption
    ftrial_diff = 1  # Trial factor for diffusion
    xs = 7.4e-10  # Diffusion coefficient constant
    a = 2.76e-10  # Length scale constant

    # Loop over the grid and calculate rate constants and diffusion rates
    for x in range(row):
        for y in range(column):
            i = neighbor[x, y]
            # Calculate desorption and adsorption rate constants
            k_des[x, y] = ftrial_ads_des * numpy.exp(alpha / 4 * (2 - i))
            k_ads[x, y] = ftrial_ads_des * numpy.exp(-alpha / 4 * (2 - i) + beta / 2)

            # Define lateral neighbors (up, down, left, right)
            lateral_neighbors = [((x + 1) % row, y), ((x - 1) % row, y), (x, (y + 1) % column), (x, (y - 1) % column)]
            k_diff[x, y] = [0, 0, 0, 0]  # Initialize diffusion rates

            # Update diffusion rates based on neighbors
            for index, position in enumerate(lateral_neighbors):
                try:
                    # Check if neighbors have matching block heights and calculate diffusion rates
                    if (block_height[position[0], position[1]] == block_height[x, y]) or (block_height[position[0], position[1]] + 1 == block_height[x, y]):
                        k_diff[x, y][index] = ftrial_diff * (xs / a) ** 2 * (numpy.exp(alpha / 4 * (2 + neighbor[position[0], position[1]] - i) - beta / 2))
                except IndexError:
                    continue

            # Sum all diffusion rates
            total_sum = numpy.sum(k_diff[x, y])
            k_x_y[x, y] = k_des[x, y] + k_ads[x, y] + total_sum  # Total event frequency at (x, y)
            k_tot = numpy.sum(k_x_y)

    # Sum the frequencies over rows
    for x in range(row):
        k_x[x] = numpy.sum(k_x_y[x, :])

    return k_des, k_ads, k_diff, k_x_y, k_tot, k_x

# Function to check the position (x, y) based on random value
def check_position(rnd_pos, k_x_y, k_x):
    x = 0
    while (x < row) and (k_x[x] < rnd_pos):
        rnd_pos -= k_x[x]
        x += 1
    xpos = x

    y = 0
    while (y < column) and (k_x_y[x, y] < rnd_pos):
        rnd_pos -= k_x_y[x, y]
        y += 1
    ypos = y
    return xpos, ypos

# Function to check the event type (adsorption, desorption, diffusion) based on random value
def check_event(rnd, k_ads, k_des, k_x_y, xpos, ypos):
    if rnd <= k_ads[xpos, ypos] / (k_x_y[xpos, ypos]):
        event = "adsorption"
    elif rnd <= (k_ads[xpos, ypos] + k_des[xpos, ypos]) / (k_x_y[xpos, ypos]):
        event = "desorption"
    else:
        event = "diffusion"
        rnd = rnd - (k_ads[xpos, ypos] + k_des[xpos, ypos]) / (k_x_y[xpos, ypos])

    return event, rnd

# Function to update block height and neighbors based on the event type
def update_height_nb(rnd, xpos, ypos, k_diff, block_height, neighbor, event):
    if event == "adsorption":
        block_height[xpos, ypos] += 1
        neighbor[xpos, ypos] = 0
        block_height, neighbor = adsorption_nb(xpos, ypos, block_height, neighbor)
        diff_pos = [1, 1]
    elif event == "desorption":
        block_height[xpos, ypos] -= 1
        neighbor[xpos, ypos] = 0
        block_height, neighbor = desorption_nb(xpos, ypos, block_height, neighbor)
        diff_pos = [1, 1]
    else:  # Diffusion
        block_height, neighbor, migrated_position = diffusion_nb(rnd, xpos, ypos, k_diff, block_height, neighbor)
        diff_pos = migrated_position[1]

    return block_height, neighbor, diff_pos

# Function to handle neighbor updates during adsorption
def adsorption_nb(xpos, ypos, block_height, neighbor):
    new_neighbors = [((xpos + 1) % row, ypos), ((xpos - 1) % row, ypos), (xpos, (ypos + 1) % column), (xpos, (ypos - 1) % column)]

    # Update neighbor counts for all adjacent blocks
    for posx, posy in new_neighbors:
        try:
            if (block_height[posx, posy] == block_height[xpos, ypos]):
                neighbor[posx, posy] += 1
                neighbor[xpos, ypos] += 1
            elif (block_height[posx, posy] == block_height[xpos, ypos] + 1):
                neighbor[xpos, ypos] += 1
        except IndexError:
            continue

    return block_height, neighbor

# Function to handle neighbor updates during desorption
def desorption_nb(xpos, ypos, block_height, neighbor):
    new_neighbors = [((xpos + 1) % row, ypos), ((xpos - 1) % row, ypos), (xpos, (ypos + 1) % column), (xpos, (ypos - 1) % column)]

    # Update neighbor counts for all adjacent blocks
    for posx, posy in new_neighbors:
        try:
            if (block_height[posx, posy] == block_height[xpos, ypos]):
                neighbor[xpos, ypos] += 1
            elif block_height[posx, posy] == block_height[xpos, ypos] + 1:
                neighbor[xpos, ypos] += 1
                neighbor[posx, posy] -= 1
        except IndexError:
            continue

    return block_height, neighbor

# Function to handle diffusion event updates
def diffusion_nb(rnd, xpos, ypos, k_diff, block_height, neighbor):
    # Define diffusion directions: down, up, right, left
    if rnd <= k_diff[xpos, ypos][0]:  # down
        # Perform the diffusion if neighbor conditions are satisfied
        if (block_height[(xpos + 1) % row, ypos] == block_height[xpos, ypos]) or (block_height[(xpos + 1) % row, ypos] + 1 == block_height[xpos, ypos]):
            neighbor[xpos, ypos] = 0
            block_height[(xpos + 1) % row, ypos] += 1
            block_height[xpos, ypos] -= 1
            neighbor[(xpos + 1) % row, ypos] = 0
            migrated_position = [(xpos, ypos), ((xpos + 1) % row, ypos)]
            block_height, neighbor = neighbor_change_diffusion(migrated_position, block_height, xpos, ypos, neighbor)
    # Repeat similarly for up, right, and left diffusion events
    # (omitted for brevity, same logic applied to other directions)

    return block_height, neighbor, migrated_position

# Function to update neighbor counts after a diffusion event
def neighbor_change_diffusion(migrated_position, block_height, xpos, ypos, neighbor):
    old_xpos = migrated_position[0][0]  # Original position (x)
    old_ypos = migrated_position[0][1]  # Original position (y)

    new_xpos = migrated_position[1][0]  # New position (x)
    new_ypos = migrated_position[1][1]  # New position (y)

    # Update neighbors for old and new positions
    old_neighbors = [((old_xpos + 1) % row, old_ypos), ((old_xpos - 1) % row, old_ypos), (old_xpos, (old_ypos + 1) % column), (old_xpos, (old_ypos - 1) % column)]
    new_neighbors = [((new_xpos + 1) % row, new_ypos), ((new_xpos - 1) % row, new_ypos), (new_xpos, (new_ypos + 1) % column), (new_xpos, (new_ypos - 1) % column)]

    # Update neighbor counts based on block heights
    for posx, posy in old_neighbors:
        try:
            if block_height[posx, posy] == block_height[xpos, ypos] + 1 or (block_height[posx, posy] == block_height[xpos, ypos]):
                neighbor[xpos, ypos] += 1
        except IndexError:
            continue

    for posx, posy in new_neighbors:
        try:
            if block_height[posx, posy] == block_height[xpos, ypos]:
                neighbor[posx, posy] += 1
                neighbor[xpos, ypos] += 1
            elif block_height[posx, posy] == block_height[xpos, ypos] + 1:
                neighbor[xpos, ypos] += 1
        except IndexError:
            continue

    return block_height, neighbor

# Function to calculate average growth based on block height
def average_growth(block_height, block_height_ref, time, d_int):
    numpy.sum(block_height)
    numpy.sum(block_height_ref)
    growth = (d_int * 10 ** -10) * (numpy.sum(block_height) - numpy.sum(block_height_ref)) / (time * numpy.sum(block_height_ref))
    print(growth)
    return growth