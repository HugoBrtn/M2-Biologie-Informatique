from random import shuffle, seed, randint

def generate_linear_conformation(hp_sequence):
    """
    Creates a linear initial configuration for an HP molecule.
    Args:
        hp_sequence (str): HP sequence (Example: "HPH")
    Returns:
        list: Linear configuration as a list of (x, y) coordinates.
              Example: generate_linear_conformation('HPH') returns [(0, 0), (1, 0), (2, 0)]
    """
    # Create a linear configuration along the x-axis starting from (0, 0)
    c_initial = [(i, 0) for i in range(len(hp_sequence))]
    return c_initial

def generate_random_conformation(hp_sequence):
    """
    Generates a random valid conformation for an HP sequence without overlaps.
    Uses backtracking to ensure a self-avoiding walk on a 2D lattice.
    Args:
        hp_sequence (str): HP sequence (e.g., "HPPHHPHPPH")
    Returns:
        list: List of (x, y) coordinates representing a valid self-avoiding conformation.
    """
    def backtrack(current_conformation, visited, remaining_length):
        """Recursive helper function to build the conformation using backtracking."""
        # Base case: all residues placed
        if remaining_length == 0:
            return current_conformation
        last_x, last_y = current_conformation[-1]
        directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]
        shuffle(directions)  # Try directions in random order
        # Try each direction
        for dx, dy in directions:
            new_pos = (last_x + dx, last_y + dy)
            # Check if the new position is not already visited
            if new_pos not in visited:
                new_conformation = current_conformation + [new_pos]
                new_visited = visited | {new_pos}
                # Recursively try to place the next residue
                result = backtrack(new_conformation, new_visited, remaining_length - 1)
                if result is not None:
                    return result
        # If no valid move found, return None to trigger backtracking
        return None

    initial_position = (0, 0)
    visited = {initial_position}
    # Try to generate a valid conformation
    result = backtrack([initial_position], visited, len(hp_sequence) - 1)
    # If no valid conformation found, try again with a different random seed
    while result is None:
        seed(randint(0, 10000))
        result = backtrack([initial_position], {initial_position}, len(hp_sequence) - 1)
    return result

def is_valid_conformation(cp):
    """
    Checks if a conformation is valid (self-avoiding and connected).
    Args:
        cp (list of tuples): List of (x, y) coordinates representing a conformation.
    Returns:
        bool: True if the conformation is valid, False otherwise.
    """
    n = len(cp)
    # Check uniqueness of positions (self-avoiding)
    if len(set(cp)) != n:
        return False
    # Check connectivity of adjacent residues
    for i in range(1, n):
        if not is_adjacent(cp[i-1], cp[i]):
            return False
    return True

def E(c, hp_sequence):
    """
    Calculates the energy of a conformation based on H-H contacts.
    Args:
        c (list of tuples): List of (x, y) coordinates of residues.
        hp_sequence (str): String representing the HP sequence (Example: "HPPH").
    Returns:
        int: Energy of the conformation (lower is better).
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
                if abs(c[idx_i][0] - c[idx_j][0]) + abs(c[idx_i][1] - c[idx_j][1]) == 1:
                    energy -= 1  # Each H-H contact contributes -1 to the energy
    return energy

def is_adjacent(pos1, pos2):
    """
    Checks if two positions are adjacent on a 2D lattice.
    Args:
        pos1 (tuple): First position as (x, y).
        pos2 (tuple): Second position as (x, y).
    Returns:
        bool: True if the positions are adjacent, False otherwise.
    """
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1]) == 1
