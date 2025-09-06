import random
from math import exp
from Neighbourhoods import *
from Others_function import *
from Grid import *
import multiprocessing

def MCsearch(hp, c=[], phi=500, nu=0.5, T=160):
    """
    Perform a Monte Carlo search to find a low-energy conformation of an HP sequence.
    Args:
        hp (str): HP sequence (Example: "HPPHHPH").
        c (list of tuples, optional): Current conformation as a list of (x, y) coordinates. If empty, a random conformation is generated.
        phi (int, optional): Number of iterations/moves to perform. Defaults to 500.
        nu (float, optional): Probability of a pull move (vs. other moves). Defaults to 0.5.
        T (float, optional): Temperature parameter for Metropolis criterion. Defaults to 160.
    Returns:
        tuple: (best_conformation, best_energy)
    """
    if c == []:
        c = generate_random_conformation(hp)
    n = len(c)
    c_mini = c.copy()  # Best conformation found
    cp = c.copy()
    c_courant = c.copy()
    Ep = E(cp, hp)  # Current energy
    E_mini = Ep  # Calculate initial energy
    for i in range(phi):
        c_courant = cp.copy()
        k = randint(0, n-1)  # Choose a random residue (1-based index)
        bool, c_courant = M(c_courant, k, nu)  # Apply a random move, nu is the probability of a pull move (instead of other moves)
        # Calculate energy differences
        E_c_courant = E(c_courant, hp)
        delta_E = E_c_courant - Ep  # Energy difference with current conformation
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
            # Metropolis criterion: accept with certain probability if energy increases
            if q > (1 / (exp(1) ** (delta_E / T))):
                cp = c_courant
                Ep = E_c_courant
    # Return best conformation found and its energy
    return c_mini, E_mini

def REMCSimulation(hp, E_star, c=[], phi=500, nu=0.5, T_init=160, T_final=220, chi=5, progress_callback=None):
    """
    Perform a Replica Exchange Monte Carlo (REMC) simulation to find a low-energy conformation of an HP sequence.
    Args:
        hp (str): HP sequence (Example: "HPPHHPH").
        E_star (int): Target energy level for the simulation.
        c (list of tuples, optional): Initial conformation as a list of (x, y) coordinates. If empty, a random conformation is generated.
        phi (int, optional): Number of iterations/moves to perform for each replica. Defaults to 500.
        nu (float, optional): Probability of a pull move (vs. other moves). Defaults to 0.5.
        T_init (float, optional): Minimum temperature. Defaults to 160.
        T_final (float, optional): Maximum temperature. Defaults to 220.
        chi (int, optional): Number of replicas to simulate. Defaults to 5.
        progress_callback (function, optional): Callback function to update progress. Defaults to None.
    Returns:
        tuple: (best_conformation, best_energy)
    """
    if c == []:
        c = generate_random_conformation(hp)

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
        print("Iteration: " + str(iteration) + " / Best energy: " + str(best_energy))

        # Perform MC search for each replica
        for k in range(chi):
            # Perform MC search
            new_conformation, new_energy = MCsearch(hp=hp, c=replicas[k][0], phi=phi, nu=nu, T=temperatures[k])
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
        # Update progress
        if progress_callback:
            progress = (iteration / max_iterations) * 100
            progress_callback(progress)
    return best_conformation, best_energy



def worker(hp, E_star, phi, nu, T_init, T_final, chi, resultat_partage, lock, index):
    """
    Worker function for multiprocessing: runs REMCSimulation with a random initial conformation.
    Args:
        hp (str): HP sequence.
        E_star (int): Target energy level.
        phi (int): Number of iterations/moves to perform for each replica.
        nu (float): Probability of a pull move.
        T_init (float): Minimum temperature.
        T_final (float): Maximum temperature.
        chi (int): Number of replicas.
        resultat_partage (list): Shared list to store results.
        lock (multiprocessing.Lock): Lock to prevent race conditions.
        index (int): Index of the process in the shared list.
    Returns:
        tuple: (best_conformation, best_energy)
    """
    c = generate_random_conformation(hp)
    best_conformation, best_energy = REMCSimulation(
        hp, E_star, c, phi, nu, T_init, T_final, chi
    )
    with lock:
        if best_energy <= E_star:
            resultat_partage[index] = (best_conformation, best_energy)
    return (best_conformation, best_energy)



def REMC_multiprocessing(hp, E_star, phi=500, nu=0.5, T_init=160, T_final=220, chi=5, nb_processus=4):
    """
    Run REMCSimulation in parallel using multiprocessing.
    Stops as soon as one process finds a conformation with energy <= E_star.
    Args:
        hp (str): HP sequence.
        E_star (int): Target energy level.
        phi (int, optional): Number of iterations/moves to perform for each replica. Defaults to 500.
        nu (float, optional): Probability of a pull move. Defaults to 0.5.
        T_init (float, optional): Minimum temperature. Defaults to 160.
        T_final (float, optional): Maximum temperature. Defaults to 220.
        chi (int, optional): Number of replicas. Defaults to 5.
        nb_processus (int, optional): Number of parallel processes. Defaults to 4.
    Returns:
        tuple: (best_conformation, best_energy)
    """
    manager = multiprocessing.Manager()
    resultat_partage = manager.list([None] * nb_processus)
    lock = manager.Lock()
    processus = []
    for i in range(nb_processus):
        p = multiprocessing.Process(
            target=worker,
            args=(hp, E_star, phi, nu, T_init, T_final, chi, resultat_partage, lock, i)
        )
        processus.append(p)
        p.start()

    # Wait for a process to find a satisfactory result
    resultat_final = None
    while True:
        with lock:
            for i in range(nb_processus):
                if resultat_partage[i] is not None:
                    resultat_final = resultat_partage[i]
                    for p in processus:
                        p.terminate()
                    return resultat_final
        # Avoid CPU overload
        multiprocessing.active_children()[0].join(timeout=0.1)

# ----- Monte Carlo Tests -----
if __name__ == "__main__":


    test = "test_REMC_multiprocessing"  # "test_REMC"   "test_MC_search"


    if test == "test_MC_search":
        hp = "HPPHHPPHHP"
        hp = "HHHHHHHHHH"
        hp = "HPHPHPHPHP"
        hp = "HPPHPPHPPH"
        c = [(2, -2), (2, -1), (2, 0), (2, 1), (1, 1),
             (1, 2), (0, 2), (0, 1), (-1, 1), (-1, 0)]
        # cp = MCsearch(hp=hp, c=c, phi=100000)
        # plot_molecules_side_by_side(c, cp[2], hp)
        # print(cp[3])
        hp = "HPHPPHHPHPPHPHHPPHPH"  # The paper indicates E*=-9
        c = generate_random_conformation(hp)
        cp = MCsearch(hp=hp, c=c, phi=100000, nu=0.5)
        plot_molecules_side_by_side(c, cp[0], hp)
        print('Energy:' + str(cp[1]))
        # print('Total changes:'+ str(cp[4]))


    elif test == "test_REMC":
        hp = "HPHPPHHPHPPHPHHPPHPH"  # The paper indicates E*=-9
        c = generate_random_conformation(hp)
        cp = REMCSimulation(hp=hp, E_star=-9, c=c, phi=1000, nu=0.5)
        plot_molecules_side_by_side(c, cp[0], hp)
        print('Energy:' + str(cp[1]))
        # print('Total changes:'+ str(cp[4]))


    elif test == "test_REMC_multiprocessing":
        hp = "HPHPPHHPHPPHPHHPPHPH"
        E_star = -9
        best_conformation, best_energy = REMC_multiprocessing(hp, E_star)
        print("Best conformation found:", best_conformation)
        print("Associated energy:", best_energy)
