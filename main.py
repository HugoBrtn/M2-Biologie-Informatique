from Monte_Carlo import *
from Grid import *

############################################################################################
###################################### PARAMETERS ##########################################
############################################################################################

#--- Method & Plotting ---------------------------------------------------------------------
method = "REMC_multi_processes"  # "REMC_multi_processes"  "MC_search"  "REMC_parallelized"
plot = True # Chose True to plot the best configuration, false otherwise

#--- Molecule Parameters -------------------------------------------------------------------
hp = "PPHPPHHPPPPHHPPPPHHPPPPHH"                # HP Sequence
E_star = -8                                               # Target Energy

#--- REMC Parallelized Method Parameters ---------------------------------------------------
phi_paral = 500                     # Iterations in Monte Carlo search
nu_paral = 0.4                      # Porbability of a pull move
T_init_paral = 160                  # Initial temperature
T_final_paral = 220                 # Final Temperature
chi_paral = 5                       # Number of replicas
max_iteration_paral = 500           # Number of maximum iteration
timeout_paral = 300                 # Timeout (in seconds)
random_initial_config = True        # If true, the initial c is random, otherwise linear

#--- REMC Multi Method Parameters ----------------------------------------------------------
phi_multi = 500                     # Iterations in Monte Carlo search
nu_multi = 0.4                      # Porbability of a pull move
T_init_multi = 160                  # Initial temperature
T_final_multi = 220                 # Final Temperature
chi_multi = 5                       # Number of replicas
max_iteration_multi = 2000           # Number of maximum iteration
nb_processus_multi = 8              # Number of simulations (differents initial conformations)
timeout_multi = 500                 # Timeout (in seconds)

#--- Monte Carlo Method Parameters ---------------------------------------------------------
phi_mc = 10000                      # Iterations in Monte Carlo search
nu_mc = 0.4                         # Probability of a pull move
T_mc = 160                          # Temperature



############################################################################################
#################################### CODE EXECUTION ########################################
############################################################################################


#--- REMC Parallelized ---------------------------------------------------------------------
if method == "REMC_parallelized":

    # Execution time calculation
    time_init = time.time()
    
    # Function
    if random_initial_config :
        best_conformation, best_energy = REMC_paral(hp=hp, E_star=E_star, phi=phi_paral,
                                                    nu=nu_paral, T_init=T_init_paral, 
                                                    T_final=T_final_paral, chi=chi_paral, 
                                                    max_iterations=max_iteration_paral, 
                                                    timeout=timeout_paral)
    else :
        best_conformation, best_energy = REMC_paral(hp=hp, c= generate_linear_conformation(hp),
                                                    E_star=E_star, phi=phi_paral,
                                                    nu=nu_paral, T_init=T_init_paral, 
                                                    T_final=T_final_paral, chi=chi_paral, 
                                                    max_iterations=max_iteration_paral, 
                                                    timeout=timeout_paral)

    execution_time = time.time() - time_init
        
    # Results
    print("execution time: " + str(execution_time))
    print("Best conformation found:", best_conformation)
    print("Associated energy:", best_energy)

    if plot :
        plot_molecule(best_conformation, hp)


#--- REMC Multi -------------------------------------------------------------------------
elif method == "REMC_multi_processes":

    # Execution time calculation
    time_init = time.time()
    
    # Function
    best_conformation, best_energy = REMC_multi(hp=hp, E_star=E_star, phi=phi_multi,
                                                nu=nu_multi, T_init=T_init_multi, 
                                                T_final=T_final_multi, chi=chi_multi, 
                                                max_iteration=max_iteration_multi, 
                                                nb_processus=nb_processus_multi, 
                                                timeout=timeout_multi)

    execution_time = time.time() - time_init
        
    # Results
    print("execution time: " + str(execution_time))
    print("Best conformation found:", best_conformation)
    print("Associated energy:", best_energy)

    if plot :
        plot_molecule(best_conformation, hp)
    

#--- Monte Carlo Search ----------------------------------------------------------------
elif method == "MC_search":        

    # Execution time calculation
    time_init = time.time()
        
    # Function
    best_conformation, best_energy = MCsearch(hp=hp, phi=phi_mc, nu=nu_mc, T=T_mc)
    execution_time = time.time() - time_init
        
    # Results
    print("execution time: " + str(execution_time))
    print("Best conformation found:", best_conformation)
    print("Associated energy:", best_energy)

    if plot :
        plot_molecule(best_conformation, hp)

