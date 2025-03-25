# Import the required functions from Definitions_sq
from Definitions_sq import *

# Define parameters
a = 7.82  # Lattice constant (Å)
b = 7.82  # Lattice constant (Å)
c = 13.54460  # Lattice constant (Å)

zeta = 0.657  # Containing the ratio of the slab to the crystal energy [-]

T_m = 273.15  # Melting temperature (K)
T = 272.15  # Current temperature (K)

# Calculate the difference in temperature (T_m - T)
dT = T_m - T

# Calculate the chemical potential (H_m is the enthalpy of fusion in kJ/mol)
H_m = 6.01 * 1000 / (6.022e23)  # Enthalpy of fusion in J (converted from kJ/mol)

# KMC parameters
dmu = H_m * dT / T_m  # Chemical potential in J
beta = dmu / (k_B * T)
alpha = zeta * H_m / (k_B * T)

# Interplanar spacing (Å)
d_int = 7.36
# Initialize iteration and time variables
iterations = 0
time = 0

# Set the end time for the simulation
time_end = 10000

# Create arrays for block heights, reference block heights, neighbors, rate constants, etc.
block_height, block_height_ref, neighbor, k_des, k_ads, k_diff, k_x_y, k_x = create_arrays()

# Initialize the neighbor array (for periodic boundary conditions)
neighbor = initialise_neighbors(neighbor)

# Update rate constants and the k_x_y array based on the current block heights and neighbors
k_des, k_ads, k_diff, k_x_y, k_tot, k_x = update_arrays(neighbor, block_height, alpha, beta, k_des, k_ads, k_diff, k_x_y, k_x)

# Main simulation loop, runs until time reaches time_end
while time < time_end:

    # Generate random uniformly distributed numbers for the simulation steps
    rnd = numpy.random.uniform(0, 1)  # Random number for event selection
    rnd_time = numpy.random.uniform(0, 1)  # Random number for time selection
    rnd_pos = numpy.random.uniform(0, k_tot)  # Random number for position selection

    # Determine the position (x, y) where the event is most likely to occur
    xpos, ypos = check_position(rnd_pos, k_x_y, k_x)

    # Check which event (adsorption, desorption, or diffusion) occurs at the selected position
    event, rnd = check_event(rnd, k_ads, k_des, k_x_y, xpos, ypos)

    # Update block heights and neighbor information based on the chosen event
    block_height, neighbor, diff_pos = update_height_nb(rnd, xpos, ypos, k_diff, block_height, neighbor, event)

    # Print the current state of the system (time, position, block height, diffusion position, and event)
    print(f"{time:<10.2e} ({xpos},{ypos})  {block_height[xpos, ypos]:<10} {diff_pos} {block_height[diff_pos[0], diff_pos[1]]} {event}\n")

    # If 10,000 iterations have passed, calculate the growth and stop the simulation
    if iterations == 10000:
        growth = average_growth(block_height, block_height_ref, time, d_int)
        break

    # Every 100 iterations, decrease the temperature slightly

    # Optionally, uncomment the following lines to stop the simulation when a layer has formed
    # if numpy.all(block_height >= 11):
    #     print('A layer has been formed after %e' %time)
    #     break  # Exit the while loop when a layer forms

    # Update rate constants and other parameters based on the current system state
    k_des, k_ads, k_diff, k_x_y, k_tot, k_x = update_arrays(neighbor, block_height, alpha, beta, k_des, k_ads, k_diff, k_x_y, k_x)

    # Increment the number of iterations and update the simulation time
    iterations += 1
    time += -numpy.log(rnd_time) / k_tot  # Time step update using the inverse of the total rate
