import pandas as pd
from random import shuffle, seed, randint



def generate_linear_conformation(hp_sequence):
    """
    Creates a linear initial configuration for an HP molecule.
    Args:
        hp_sequence (str): HP sequence (Example: "HPH")
    Returns:
        pd.DataFrame: Linear configuration as a DataFrame with columns 'x' and 'y'
    """
    # Create a DataFrame with x and y coordinates
    c_initial = pd.DataFrame({
        'x': range(len(hp_sequence)),
        'y': [0] * len(hp_sequence)
    })
    return c_initial



def generate_random_conformation(hp_sequence):
    """
    Generates a random valid conformation for an HP sequence without overlaps.
    Uses backtracking to ensure a self-avoiding walk on a 2D lattice.
    Args:
        hp_sequence (str): HP sequence (e.g., "HPPHHPHPPH")
    Returns:
        pd.DataFrame: DataFrame with columns 'x' and 'y' representing a valid self-avoiding conformation
    """
    def backtrack(current_conformation, visited, remaining_length):
        """Recursive helper function to build the conformation using backtracking."""
        # Base case: all residues placed
        if remaining_length == 0:
            return current_conformation
        last_x, last_y = current_conformation.iloc[-1]
        directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]
        shuffle(directions)  # Try directions in random order
        # Try each direction
        for dx, dy in directions:
            new_pos = (last_x + dx, last_y + dy)
            # Check if the new position is not already visited
            if new_pos not in visited:
                new_conformation = pd.concat([
                    current_conformation,
                    pd.DataFrame({'x': [new_pos[0]], 'y': [new_pos[1]]}),
                ], ignore_index=True)
                new_visited = visited | {new_pos}
                # Recursively try to place the next residue
                result = backtrack(new_conformation, new_visited, remaining_length - 1)
                if result is not None:
                    return result
        # If no valid move found, return None to trigger backtracking
        return None

    initial_position = (0, 0)
    visited = {initial_position}
    initial_conformation = pd.DataFrame({'x': [0], 'y': [0]})
    # Try to generate a valid conformation
    result = backtrack(initial_conformation, visited, len(hp_sequence) - 1)
    # If no valid conformation found, try again with a different random seed
    while result is None:
        seed(randint(0, 10000))
        result = backtrack(initial_conformation, {initial_position}, len(hp_sequence) - 1)
    return result



def is_valid_conformation(cp):
    """Vérifie si une conformation est valide (auto-évitante et connectée)."""
    n = len(cp)
    # Vérifier l'unicité des positions (auto-évitant)
    if len(cp.drop_duplicates()) != n:
        return False
    # Vérifier la connectivité des résidus adjacents
    for i in range(1, n):
        if not is_adjacent(cp.iloc[i-1], cp.iloc[i]):
            return False
    return True



def E(c, hp_sequence):
    """
    Args:
        c (pd.DataFrame): DataFrame with columns 'x' and 'y' representing the conformation.
        hp_sequence (str): String representing the HP sequence (Example: "HPPH").
    Returns:
        int: Energy of the conformation.
    """
    # List of indices of H residues
    h_indices = [i for i, residue in enumerate(hp_sequence) if residue == 'H']
    energy = 0
    # Iterate through all pairs of non-consecutive H residues
    for i in range(len(h_indices)):
        for j in range(i + 1, len(h_indices)):
            idx_i = h_indices[i]
            idx_j = h_indices[j]
            # Check if residues are not consecutive in the sequence
            if abs(idx_i - idx_j) > 1:
                # Check if residues are adjacent on the lattice
                if abs(c.iloc[idx_i]['x'] - c.iloc[idx_j]['x']) + abs(c.iloc[idx_i]['y'] - c.iloc[idx_j]['y']) == 1:
                    energy -= 1  # Each H-H contact contributes -1 to the energy
    return energy

def is_adjacent(pos1, pos2):
    """Checks if two positions are adjacent on a 2D lattice."""
    return abs(pos1['x'] - pos2['x']) + abs(pos1['y'] - pos2['y']) == 1
