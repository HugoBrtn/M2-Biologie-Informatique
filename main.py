from Monte_Carlo import *
from Grid import *


# ---------- Monte Carlo / REMC ----------
if __name__ == "__main__":

    method = "REMC_multiprocessing"  # "REMC_multiprocessing"   "MC_search"
    plot = True

    # Molecule Parameters
    hp = "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP"
    E_star = -14

    # -- REMC_multiprocessing -----
    if method == "REMC_multiprocessing":

        # REMC Method Parameters
        phi = 500                           # Iterations in Monte Carlo search
        nu = 0.4                            # Porbability of a pull move
        T_init = 160
        T_final = 220
        chi = 5
        max_iteration = 500
        nb_processus = 8
        timeout = 300                        # seconds

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

        # REMC Method Parameters
        phi = 500                           # Iterations in Monte Carlo search
        nu = 0.4                            # Porbability of a pull move
        T = 160

        # Execution time calculation
        time_init = time.time()
        
        # Function
        best_conformation, best_energy = MCsearch(hp=hp, c=c, phi=phi, nu=nu, T=T)
        execution_time = time.time() - time_init
        
        # Results
        print("execution time: " + str(execution_time))
        print("Best conformation found:", best_conformation)
        print("Associated energy:", best_energy)

        if plot :
            plot_molecule(best_conformation, hp)

