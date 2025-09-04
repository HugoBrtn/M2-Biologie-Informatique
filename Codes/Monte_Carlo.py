import random
from math import exp
from Neighbourhoods import *
from Others_function import *
from Grid import *



def MCsearch(phi, c, hp, nu, T=160):
    """
    Perform a Monte Carlo search to find a low-energy conformation of an HP sequence.

    Args:
        phi (int): Number of iterations/moves to perform.
        c (list of tuples): Current conformation as a list of (x, y) coordinates.
        hp (str): HP sequence (Ewample : "HPPHHPH").
        T (float, optional): Temperature parameter for Metropolis criterion. Defaults to 220.

    Returns:
        tuple: (final_conformation, final_energy, best_conformation, best_energy)
    """
    n = len(c)
    c_mini = c.copy()  # Best conformation found
    cp = c.copy() 
    c_courant = c.copy()
    Ep = E(cp, hp)  # Current energy
    E_mini = Ep  # Calculate initial energy

    for i in range(phi):
        c_courant = cp.copy()
        k = randint(0, n-1)  # Choose a random residue (1-based index)
        bool, c_courant = M(c_courant, k, nu)  # Apply arandom move, nu is the probability of a pull move (instead vhsd move)

        # Calculate energy differences
        E_c_courant = E(c_courant, hp)
        delta_E = E_c_courant - Ep  # Energy difference with best conformation found

        # Always accept if energy decreases or stays the same
        if delta_E <= 0:
            cp = c_courant
            Ep = E_c_courant

            # Update best conformation if this one is better
            if E_c_courant - E_mini < 0:
                c_mini = c_courant
                E_mini = E_c_courant

        else:
            q = random.random()  # Generate a random number between 0 and 1
            # Metropolis criterion : accept with certain probability if energy increases
            if q > (1 / (exp(1) ** (delta_E / T))):
                cp = c_courant
                Ep = E_c_courant
    # Return final conformation, its energy, best conformation found, and its energy
    return c_mini, E_mini




def REMCSimulation(phi, c, hp, nu, E_star, T_init=160, T_final=220, chi=5):
    """
    Perform a Replica Exchange Monte Carlo (REMC) simulation to find a low-energy conformation of an HP sequence.

    Args:
        phi (int): Number of iterations/moves to perform for each replica.
        c (list of tuples): Initial conformation as a list of (x, y) coordinates.
        hp (str): HP sequence (Example: "HPPHHPH").
        E_star (int): Target energy level for the simulation.
        T_init (float, optional): Minimum temperature. Defaults to 160.
        T_final (float, optional): Maximum temperature. Defaults to 220.
        chi (int, optional): Number of replicas to simulate. Defaults to 5.

    Returns:
        tuple: (best_conformation, best_energy)
    """
    # Create linear temperature schedule
    temperatures = [T_init + i * (T_final - T_init) / (chi - 1) for i in range(chi)]

    # Initialize replicas with the same initial conformation and energy
    replicas = [(c.copy(), E(c, hp)) for _ in range(chi)]

    # Track best conformation and energy
    best_conformation = c.copy()
    best_energy = E(c, hp)
    offset = 0

    # Maximum number of iterations to prevent infinite loops
    max_iterations = 500
    iteration = 0
    

    while best_energy > E_star and iteration < max_iterations:
        iteration += 1
        print("itération : " +str(iteration) + " /  best_energy : " + str(best_energy))
        

        # Perform MC search for each replica
        for k in range(chi):
            # Perform MC search
            new_conformation, new_energy = MCsearch(phi, replicas[k][0], hp, nu, temperatures[k])
            replicas[k] = (new_conformation, new_energy)

            # Update best conformation if needed
            if new_energy < best_energy:
                best_conformation = new_conformation.copy()
                best_energy = new_energy

        # Attempt replica exchanges between neighboring temperatures
        i = offset + 1
        while i + 1 < chi:  # Fixed condition to avoid index errors
            j = i + 1

            # Calculate exchange probability
            delta = (1/temperatures[j] - 1/temperatures[i]) * (replicas[i][1] - replicas[j][1])

            # Accept exchange with Metropolis criterion
            if delta <= 0:
                # Always accept if energy difference is favorable
                replicas[i], replicas[j] = replicas[j], replicas[i]
            else:
                # Accept with probability if energy difference is unfavorable
                if random.random() < exp(-delta):
                    replicas[i], replicas[j] = replicas[j], replicas[i]

            i += 2  # Move to next pair of replicas

        # Toggle offset for next iteration
        offset = 1 - offset

    return best_conformation, best_energy




# ----- Tests Montecarlo -----
if __name__ == "__main__" :

    test = "test_REMC"

    if test == "test_MC_search" :
        hp = "HPPHHPPHHP"
        hp = "HHHHHHHHHH"
        hp = "HPHPHPHPHP"
        hp = "HPPHPPHPPH"
        c = [(2, -2), (2, -1), (2, 0), (2, 1), (1, 1), 
            (1, 2), (0, 2), (0, 1), (-1, 1), (-1, 0)]
        #cp = MCsearch(100000,c,hp)
        #plot_molecules_side_by_side(c, cp[2], hp)
        #print(cp[3])

        hp = "HPHPPHHPHPPHPHHPPHPH"   # L'article indique E*=-9
        c = generate_random_conformation(hp)
        cp = MCsearch(100000,c,hp, 0.5)
        plot_molecules_side_by_side(c, cp[0], hp)
        print('Energie:'+ str(cp[1]))
        #print('Total changes:'+ str(cp[4]))

    elif test == "test_REMC" :
        hp = "HPHPPHHPHPPHPHHPPHPH"   # L'article indique E*=-9
        c = generate_random_conformation(hp)
        cp = REMCSimulation(1000,c,hp, 0.5, -9)
        plot_molecules_side_by_side(c, cp[0], hp)
        print('Energie:'+ str(cp[1]))
        #print('Total changes:'+ str(cp[4]))




