import random
from math import exp
from Neighbourhoods import *
from Others_function import *
from Grid import *
import multiprocessing
import time



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
    # c_mini = c.copy()  # Best conformation found
    cp = c.copy()
    c_courant = c.copy()
    Ep = E(cp, hp)  # Current energy
    # E_mini = Ep  # Calculate initial energy

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
            # if E_c_courant - E_mini < 0:
            #     c_mini = c_courant
            #     E_mini = E_c_courant
        else:
            q = random.random()  # Generate a random number between 0 and 1

            # Metropolis criterion: accept with certain probability if energy increases
            if q > (1 / (exp(1) ** (delta_E / T))):
                cp = c_courant
                Ep = E_c_courant

    # Return best conformation found and its energy
    return c_courant, E_c_courant #c_mini, E_mini




def REMCSimulation(hp, E_star, c=[], phi=500, nu=0.5, T_init=160, T_final=220, chi=5, max_iterations=300, timeout =300):
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
    iteration = 0

    # Timeout initialization
    time_start = time.time()

    while best_energy > E_star and iteration < max_iterations and time.time() - time_start < timeout :
        iteration += 1
        print(f"Iteration {iteration}, Best Energy: {best_energy}")

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
        while i + 1 < chi:
            j = i + 1

            # Calculate exchange probability
            delta = (1/temperatures[j] - 1/temperatures[i]) * (replicas[i][1] - replicas[j][1])

            # Accept exchange with Metropolis criterion
            if delta <= 0:
                replicas[i], replicas[j] = replicas[j], replicas[i]
                temperatures[i], temperatures[j] = temperatures[j], temperatures[i]
            else:
                if random.random() <= exp(-delta):
                    replicas[i], replicas[j] = replicas[j], replicas[i]
                    temperatures[i], temperatures[j] = temperatures[j], temperatures[i]
            i += 2

        # Toggle offset for next iteration
        offset = 1 - offset

    return best_conformation, best_energy



def worker_REMC_multi(hp, E_star, phi, nu, T_init, T_final, chi, max_iteration, timeout, resultat_partage, lock, index):
    """
    Worker function for multiprocessing: runs REMC Simulation with a random initial conformation.
    Always stores the best conformation found.
    """
    c = generate_random_conformation(hp)
    best_conformation, best_energy = REMCSimulation(hp=hp, E_star=E_star, c=c, phi=phi, nu=nu, T_init=T_init, 
                                                    T_final=T_final, chi=chi, max_iterations=max_iteration, 
                                                    timeout=timeout)
    with lock:
        resultat_partage[index] = (best_conformation, best_energy)
    return (best_conformation, best_energy)



def REMC_multi(hp, E_star, phi=500, nu=0.5, T_init=160, T_final=220, chi=5, max_iteration = 300,  nb_processus=4, timeout = 300):
    """
    Run REMC Simulation in parallel using multiprocessing for calculating REMC for different initial configurations.
    Returns the best conformation found, even if no conformation satisfies E_star.
    """
    manager = multiprocessing.Manager()
    resultat_partage = manager.list([None] * nb_processus)  # Shared list for results
    lock = manager.Lock()  # Lock to avoid race conditions
    processus = []  # List to store processes

    # Start all processes
    for i in range(nb_processus):
        p = multiprocessing.Process(
            target=worker_REMC_multi,
            args=(hp, E_star, phi, nu, T_init, T_final, chi, max_iteration, timeout, resultat_partage, lock, i)
        )
        processus.append(p)
        p.start()

    # Wait for all processes to finish or a solution to be found
    best_conformation, best_energy = generate_linear_conformation(hp), 0

    while True:
        with lock:
            for i in range(nb_processus):
                if resultat_partage[i] is not None:
                    conformation, energy = resultat_partage[i]
                    if energy < best_energy:
                        best_conformation, best_energy = conformation, energy
                    if energy <= E_star:
                        # A solution was found: terminate all processes and return
                        for p in processus:
                            p.terminate()
                        return best_conformation, best_energy

        # Check if all processes are done
        all_done = all(not p.is_alive() for p in processus)
        if all_done:
            break

        time.sleep(0.1)

    # Terminate any remaining processes
    for p in processus:
        if p.is_alive():
            p.terminate()

    # Return the best conformation found
    return best_conformation, best_energy



def worker_MCsearch(args):
    """Wrapper to call MCsearch with correct arguments."""
    hp, c, phi, nu, T = args
    return MCsearch(hp=hp, c=c, phi=phi, nu=nu, T=T)


def REMC_paral(hp, E_star, c=[], phi=500, nu=0.5, T_init=160, T_final=220, chi=5, max_iterations=300, timeout=300):
    """
    Perform a Replica Exchange Monte Carlo (REMC) simulation to find a low-energy conformation of an HP sequence.
    This version uses multiprocessing to parallelize the MC search for each replica.
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
    iteration = 0
    # Timeout calculation
    start_time = time.time() 

    while best_energy > E_star and iteration < max_iterations and time.time() - start_time < timeout:
        iteration += 1
        print(f"Iteration {iteration}, Best Energy: {best_energy}")

        # Parallelize the MC search for each replica using multiprocessing
        with multiprocessing.Pool() as pool:
            # Prepare arguments for each worker
            args = [(hp, replicas[k][0], phi, nu, temperatures[k]) for k in range(chi)]
            # Map the work to the pool
            results = pool.map(worker_MCsearch, args)
            # Update replicas with results
            for k, (new_conformation, new_energy) in enumerate(results):
                replicas[k] = (new_conformation, new_energy)
                # Update best conformation if needed
                if new_energy < best_energy:
                    best_conformation = new_conformation.copy()
                    best_energy = new_energy

        # Attempt replica exchanges between neighboring temperatures
        i = offset + 1
        while i + 1 < chi:
            j = i + 1
            # Calculate exchange probability
            delta = (1/temperatures[j] - 1/temperatures[i]) * (replicas[i][1] - replicas[j][1])
            # Accept exchange with Metropolis criterion
            if delta <= 0:
                replicas[i], replicas[j] = replicas[j], replicas[i]
                temperatures[i], temperatures[j] = temperatures[j], temperatures[i]
            else:
                if random.random() < exp(-delta):
                    replicas[i], replicas[j] = replicas[j], replicas[i]
                    temperatures[i], temperatures[j] = temperatures[j], temperatures[i]
            i += 2
        # Toggle offset for next iteration
        offset = 1 - offset
    return best_conformation, best_energy




# ----- Monte Carlo / REMC Tests -----
if __name__ == "__main__":

    test = "test_REMC_paral"  # "test_REMC_multiprocessing"   "test_MC_search"   ""test_REMC_paral""

    # -- Test MCsearch -----
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

    # -- Test REMCSimulation -----
    elif test == "test_REMC":
        hp = "PPHPPHHPPPPHHPPPPHHPPPPHH"  
        c = generate_linear_conformation(hp)
        cp = REMCSimulation(hp=hp, E_star=-7, c=c, phi=500, nu=0.5)
        plot_molecules_side_by_side(c, cp[0], hp)
        print('Energy:' + str(cp[1]))
        # print('Total changes:'+ str(cp[4]))

    # -- Test REMC_multi -----
    elif test == "test_REMC_multi":
        # S1-1
        #hp = "HPHPPHHPHPPHPHHPPHPH"
        #E_star = -9  

        # S1-3
        #hp = "PPHPPHHPPPPHHPPPPHHPPPPHH"
        #E_star = -8
        
        # S1-7
        #hp = "PPHHHPHHHHHHHHPPPHHHHHHHHHHPHPPPHHHHHHHHHHHHPPPPHHHHHHPHHPHP"
        #E_star = -36 

        # S1-4
        hp = "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP"
        E_star = -14

        best_conformation, best_energy = REMC_multi(hp, E_star, max_iteration=100, nb_processus = 4)
        print("Best conformation found:", best_conformation)
        print("Associated energy:", best_energy)
        plot_molecule(best_conformation, hp)
    
    # -- Test REMC_multi -----
    elif test == "test_REMC_paral":
        # S1-1
        #hp = "HPHPPHHPHPPHPHHPPHPH"
        #E_star = -9  

        # S1-3
        #hp = "PPHPPHHPPPPHHPPPPHHPPPPHH"
        #E_star = -8
        
        # S1-7
        #hp = "PPHHHPHHHHHHHHPPPHHHHHHHHHHPHPPPHHHHHHHHHHHHPPPPHHHHHHPHHPHP"
        #E_star = -36 

        # S1-4
        hp = "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP"
        E_star = -14

        best_conformation, best_energy = REMC_paral(hp, E_star, max_iterations=00)
        print("Best conformation found:", best_conformation)
        print("Associated energy:", best_energy)
        plot_molecule(best_conformation, hp)
