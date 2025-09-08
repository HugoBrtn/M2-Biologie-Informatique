from Monte_Carlo import *
from Grid import *


# ---------- Monte Carlo / REMC ----------
if __name__ == "__main__":

    method = "MC_search"  # "REMC_multiprocessing"   "MC_search"
    plot = True # Chose True to plot the best configuration, false otherwise

    # Molecule Parameters
    hp = "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP"
    E_star = -14

    # -- REMC_multiprocessing -----
    if method == "REMC_multiprocessing":

        # REMC Method Parameters
        phi = 500                           # Iterations in Monte Carlo search
        nu = 0.4                            # Porbability of a pull move
        T_init = 160                        # Initial temperature
        T_final = 220                       # Final Temperature
        chi = 5                             # Number of replicas
        max_iteration = 500                 # Number of maximum iteration
        nb_processus = 8                    # Number of simulations
        timeout = 300                       # Timeout in seconds

        # Execution time calculation
        time_init = time.time()
        
        # Function
        best_conformation, best_energy = REMC_multiprocessing(hp=hp, E_star=E_star, phi=phi,
                                                              nu=nu, T_init=T_init, T_final=T_final,
                                                              chi=chi, max_iteration=max_iteration, 
                                                              nb_processus=nb_processus, 
                                                              timeout=timeout)
        execution_time = time.time() - time_init
        
        # Results
        print("execution time: " + str(execution_time))
        print("Best conformation found:", best_conformation)
        print("Associated energy:", best_energy)

        if plot :
            plot_molecule(best_conformation, hp)
    

    # -- REMC_multiprocessing -----
    elif method == "MC_search":

        # Monte Carlo Method Parameters
        phi = 10000                         # Iterations in Monte Carlo search
        nu = 0.4                            # Probability of a pull move
        T = 160                             # Temperature

        # Execution time calculation
        time_init = time.time()
        
        # Function
        best_conformation, best_energy = MCsearch(hp=hp, phi=phi, nu=nu, T=T)
        execution_time = time.time() - time_init
        
        # Results
        print("execution time: " + str(execution_time))
        print("Best conformation found:", best_conformation)
        print("Associated energy:", best_energy)

        if plot :
            plot_molecule(best_conformation, hp)

