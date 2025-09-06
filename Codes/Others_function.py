from random import shuffle, seed, randint
from Grid import *



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
    Uses an iterative approach with a stack to simulate recursion.
    Args:
        hp_sequence (str): HP sequence (e.g., "HPPHHPHPPH")
    Returns:
        list: List of (x, y) coordinates representing a valid self-avoiding conformation, or None if no solution is found.
    """
    if not hp_sequence:
        return []

    # Initialization: Start at the origin (0, 0)
    start = (0, 0)
    current_conformation = [start]
    visited = {start}
    # Stack stores tuples of (current_conformation, visited_positions, remaining_length)
    stack = []
    stack.append((current_conformation, visited, len(hp_sequence) - 1))

    while stack:
        # Pop the most recent state from the stack
        current_conformation, visited, remaining_length = stack.pop()

        # Base case: all residues are placed
        if remaining_length == 0:
            return current_conformation

        # Get the last position in the current conformation
        last_x, last_y = current_conformation[-1]
        # Possible directions: right, left, down, up
        directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]
        # Shuffle directions to explore them in random order
        shuffle(directions)

        # Try each direction
        for dx, dy in directions:
            new_pos = (last_x + dx, last_y + dy)
            # Check if the new position is not already visited
            if new_pos not in visited:
                # Add the new position to the conformation and mark it as visited
                new_conformation = current_conformation + [new_pos]
                new_visited = visited | {new_pos}
                # Push the new state onto the stack
                stack.append((new_conformation, new_visited, remaining_length - 1))

    # If the stack is empty and no solution was found, return None
    return None




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
                    energy -= 1  # Each H-H contact contributes -1 to the 
                    
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


import re

def expand_hp_sequence(compact_notation):
    """
    Expand a compact HP sequence notation into its full form.

    Args:
        compact_notation (str): Compact notation of HP sequence (e.g., "P3H2P2H2P5H7P2H2P4H2P2HP2")

    Returns:
        str: Expanded HP sequence
    """
    # First expand all (pattern)number notations
    while True:
        match = re.search(r'\(([^)]+)\)(\d+)', compact_notation)
        if not match:
            break
        pattern = match.group(1)
        repeat = int(match.group(2))
        expanded_pattern = pattern * repeat
        compact_notation = compact_notation.replace(match.group(0), expanded_pattern, 1)

    # Then expand all Hnumber/Pnumber notations
    expanded = []
    i = 0
    n = len(compact_notation)
    while i < n:
        if compact_notation[i] in ('H', 'P'):
            # Check if there's a number following the H or P
            j = i + 1
            num_str = ''
            while j < n and compact_notation[j].isdigit():
                num_str += compact_notation[j]
                j += 1

            if num_str:
                # Expand the notation (e.g., H3 becomes HHH)
                count = int(num_str)
                expanded.append(compact_notation[i] * count)
                i = j
            else:
                # Single character, no number
                expanded.append(compact_notation[i])
                i += 1
        else:
            # Skip any non-H/P characters (shouldn't happen in valid input)
            i += 1

    return ''.join(expanded)




# ----- Others functions Tests -----
if __name__ == "__main__":

    test = "random_conformation"    #  "linear_conformation"  # "energy"    # "linear_conformation"    #  "expanded"

    # ----- Test Energy -----
    if test == "energy":
        hp_sequence = "HPPHHPHPPH"
        c = [(0, 0), (1, 0), (1, 1), (0, 1), (0, 2),
            (1, 2), (1, 3), (0, 3), (0, 4), (1, 4)]
        energy = E(c, hp_sequence)
        print(f"Énergie de la conformation : {energy}")
        #plot_molecule(c, hp_sequence)

    # -- Test Random Conformation -----
    elif test == "random_conformation" :
        hp_sequence = "HPPHHPHPPHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPH"
        c = generate_random_conformation(hp_sequence)
        print(f"Random conformation : {c}")
        plot_molecule(c, hp_sequence)

    # -- Test Linear Conformation -----
    elif test == "linear_conformation" :
        hp_sequence = "HPPHHPHPPH"
        c = generate_linear_conformation(hp_sequence)
        print(f"Linear conformation : {c}")
        plot_molecule(c, hp_sequence)

    elif test == "expanded" :
        # Test the function with the provided sequences
        sequences = [
            "P3H2P2H2P5H7P2H2P4H2P2HP2",
            "P2H(P2H2)2P5H10P6(H2P2)2HP2H5",
            "H2(PH)3PH4PH(P3H)2P4H(P3H)2PH4(PH)4H",
            "P2H3PH8P3H10PHP3H12P4H6PH2PHP",
            "H12(PH)2(P2H2)2P2HP2H2PPH2P2HP2(H2P2)2(HP)2H12",
            "H4P4H12P6(H12P3)3HP2(H2P2)2HPH",
            "P3H2P2H4P2H3(PH2)2PH4P8H6P2H6P9HPH2PH11P2H3PH2PHP2HPH3P6H3",
            "P6HPH2P5H3PH5PH2P4H2P2H2PH5PH10PH2PH7P11H7P2HPH3P6HPH2"
        ]

        for seq in sequences:
            expanded = expand_hp_sequence(seq)

            print(f"{seq}\n-> {expanded}\n")
            print("Size : " + str(len(expanded)))