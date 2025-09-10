# Monte Carlo lattice Project
### Description
This project aims to find the lowest-energy structure of a molecule based on its HP sequence modeled on a 2D lattice. The HP sequence simplifies the molecule by classifying each atom as either hydrophobic (H) or polar (P). The energy is calculated by counting the number of adjacent hydrophobic atom pairs on the lattice that are not consecutive in the sequence: each such pair contributes -1 to the total energy.

### Installation
To use this project, you must first install uv on your machine. Then, you can clone this repository and use uv commands to manage the environment and dependencies.

**Downloading the Project**
```bash
git clone https://github.com/HugoBrtn/Monte-Carlo-lattice.git
cd Monte-Carlo-lattice
uv sync
```
You can then use either the interface or the script to perform the functions.

**Interface**\
To run the interface, use the following command:
```bash
uv run interface.py
```

**Script**\
If you prefer to use a script directly, use this command:
```bash
uv run main.py
```

### Functions
**Monte Carlo Search (MC Search)**\
This function uses the Monte Carlo method to estimate the lowest-energy configuration.

Parameters :
- phi : Number of Monte Carlo iterations.
- nu : Probability of performing a pull move (otherwise VSHD move).
- T : Temperature.

\
**REMC on Multi-Initial Configuration (MC Search)**\
This function uses the REMC (Replica Exchange Monte Carlo) method to estimate the lowest-energy configuration.
It is parallelized to run REMC on multiple different random initial configurations of the molecule.

Parameters :
- phi : Number of Monte Carlo iterations.
- nu : Probability of performing a pull move (otherwise VSHD move).
- T_init : Temperature of the first replica.
- T_final : Temperature of the last replica.
- chi : Number of replicas.
- max_iteration : Maximum number of iterations for REMC.
- timeout : Maximum runtime before the program terminates.

\
**REMC with Parallelization for Replicas**\
This function uses the REMC method to estimate the lowest-energy configuration.
It is parallelized to execute the Monte Carlo simulations for each replica separately.

Parameters :
- phi : Number of Monte Carlo iterations.
- nu : Probability of performing a pull move (otherwise VSHD move).
- T_init : Temperature of the first replica.
- T_final : Temperature of the last replica.
- chi : Number of replicas.
- max_iteration : Maximum number of iterations for REMC.
- timeout : Maximum runtime before the program terminates.
