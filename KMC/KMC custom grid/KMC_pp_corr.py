from Definitions_pp_corr import *

zeta = 0.929 # [-]

T_m = 273.15  # melting temp [K] or 0 degr Celsius (pure water)
T = 272.15

while T >=263.15:
    dT = T_m - T  # difference in temperature [K]
    H_m = 6.01 * 1000 / (6.022e23)  # Enthalpy of fusion in J (converted from kJ/mol)
    dmu = H_m * dT / T_m  # chemical potential [J]

    # KMC parameters
    beta = dmu / (k_B * T)
    alpha = zeta * H_m / (k_B * T)

    # Interplanar spacing (A)
    d_int = 7.82
    iterations = 0
    time = 0  # counter is initially set to 0
    # Load lateral neighbors
    with open("neighborspp.json", "r") as f:
            lateral_neighbors = json.load(f)
    lateral_neighbors = {int(k): v for k, v in lateral_neighbors.items()}

    print("Neighbors loaded successfully!")
    time_end = 10000
    grid_points = create_grid()
    block_height, block_height_ref,neighbor, k_des, k_ads, k_diff, k_diff_hor, k_diff_ab, k_diff_below, k_x_y = create_arrays(grid_points)
    neighbor = initialise_neighbors(neighbor)
    k_des, k_ads, k_diff, k_diff_hor, k_diff_ab, k_diff_below, k_x_y, k_tot = create_rates(grid_points,neighbor,block_height, alpha, beta, k_des, k_ads, k_diff, k_diff_hor, k_diff_ab, k_diff_below, k_x_y,lateral_neighbors)

    while time < time_end:

        # Generation random uniformly distributed numbers
        rnd = numpy.random.uniform(0, 1)
        rnd_time = numpy.random.uniform(0, 1)
        rnd_pos = numpy.random.uniform(0, k_tot)
        # Choose the position where the event is most likely to occur
        xpos = check_position(rnd_pos, k_x_y)
        # Check which event (ads,des) occurs here
        event,rnd = check_event(rnd, k_ads, k_des, k_x_y, xpos)
        block_height, neighbor, diff_pos = update_height_nb(rnd, xpos,k_diff_hor, k_diff_ab, k_diff_below, block_height, neighbor, event, lateral_neighbors)
        # save_grid_csv(block_height,time)

        print(f"{time:<10e} {xpos}   {block_height[xpos]} {diff_pos} {block_height[diff_pos]} {event}")
        if iterations == 10000:
            growth = average_growth(block_height, block_height_ref, time, d_int)
            break


        # if block gets added, check if more blocks are already present next to it on the same level and update if so
        k_des, k_ads, k_diff, k_diff_hor, k_diff_ab, k_x_y, k_tot = update_arrays(grid_points,neighbor, block_height, alpha, beta, k_des, k_ads,k_diff, k_diff_hor, k_diff_ab, k_diff_below, k_x_y, lateral_neighbors)
        iterations +=1
        time += -numpy.log(rnd_time) / k_tot

    T = T-1