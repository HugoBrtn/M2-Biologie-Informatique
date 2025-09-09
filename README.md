# Monte Carlo lattice Project
## Description
This project aims to find the lowest-energy structure of a molecule based on its HP sequence modeled on a 2D grid. The HP sequence simplifies the molecule by classifying each atom as either hydrophobic (H) or polar (P). The energy is calculated by counting the number of adjacent hydrophobic atom pairs on the grid that are not consecutive in the sequence; each such pair contributes -1 to the total energy.

## Functions
### Monte Carlo Search (MC Search)
This function uses the Monte Carlo method to estimate the lowest-energy configuration.

Parameters :
- 

### REMC on Multi-Initial Configuration (MC Search)
This function uses the REMC (Replica Exchange Monte Carlo) method to estimate the lowest-energy configuration.
It is parallelized to run REMC on multiple different random initial configurations of the molecule.

### REMC with Parallelization for Replicas
This function uses the REMC method to estimate the lowest-energy configuration.
It is parallelized to execute the Monte Carlo simulations for each replica separately.

## Installation
To use this project, you must first install uv on your machine. Then, you can clone this repository and use uv commands to manage the environment and dependencies.

### Downloading the Project

```bash
git clone https://github.com/HugoBrtn/Monte-Carlo-lattice.git
cd Monte-Carlo-lattice
uv sync
uv run interface.py
```

### Interface

```bash
uv run interface.py
```


### Code
```bash
uv run main.py
```
