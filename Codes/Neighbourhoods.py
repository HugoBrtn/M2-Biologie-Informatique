import random
from Others_function import *
from Grid import *



def M(c, k, nu):
    """
    Applies a random move (either pull or VSHD) to residue k based on probability nu.
    Args:
        c (list of tuples): Current conformation as a list of (x, y) coordinates.
        k (int): Index of the residue to move.
        nu (float): Probability of applying a pull move (vs. VSHD move).
    Returns:
        tuple: (bool, new_conformation)
            bool: True if the move was successful, False otherwise.
            new_conformation: The new conformation after the move.
    """
    rand = random.random()
    if rand < nu:
        return pull_move(c, k)
    return M_vshd(c, k)



def M_vshd(c, k):
    """
    Applies a VSHD move (end, corner, or crankshaft) to residue k.
    Args:
        c (list of tuples): Current conformation as a list of (x, y) coordinates.
        k (int): Index of the residue to move (must be between 1 and n-3).
    Returns:
        tuple: (bool, new_conformation)
            bool: True if the move was successful, False otherwise.
            new_conformation: The new conformation after the move (or the original if no move was possible).
    """
    c_prime = c.copy()  # Create a copy of the conformation
    n = len(c_prime)

    # Case 1: End move (if k is the first or last residue)
    if k == 0 or k == n-1:
        end_move_possible, new_c = end_move(c_prime, k)
        if end_move_possible:
            return (True, new_c)
        
    # Case 2: Corner move (if k is the second-to-last residue)
    elif k == n-2:
        corner_possible, new_c = corner_move(c_prime, k)
        if corner_possible:
            return (True, new_c)

    # Case 3: For internal residues, try corner or crankshaft move
    else:
        rand = randint(1, 2)  # Randomly choose between corner and crankshaft
        corner_possible, new_c_corner = corner_move(c_prime, k)
        crankshaft_possible, new_c_crankshaft = crankshaft_move(c_prime, k)

        # Try corner move first if randomly selected
        if rand == 1 and corner_possible:
            return (True, new_c_corner)
        
        # Try crankshaft move if possible
        if crankshaft_possible:
            return (True, new_c_crankshaft)

        # If crankshaft not possible but corner is, use corner move
        if corner_possible:
            return (True, new_c_corner)
        
    # If no move is possible, return the unchanged conformation
    return (False, c_prime)



def end_move(c, k):
    """
    Applies an end move to residue k, where k must be first (0) or last residue (n-1).
    Args:
        c (list of tuples): Current conformation.
        k (int): Index of the residue to move (0 or n-1).
    Returns:
        tuple: (bool, new_conformation)
            bool: True if the move was successful, False otherwise.
            new_conformation: The new conformation after the move (or the original if no move was possible).
    """

    cp = c.copy()
    if k == 0:
        neighbour_residue = c[1]  # Residue 1
    else:  # k == n-1
        neighbour_residue = c[-2]  # Residue n-
        
    # Possible directions in a 2D lattice:
    x_nr = neighbour_residue[0]
    y_nr = neighbour_residue[1]
    directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]
    random.shuffle(directions)  # Shuffle to get a random position

    for x, y in directions:
        # Test all neighbor positions of neighbour_residue and return the first one which is empty
        new_x = x_nr + x
        new_y = y_nr + y
        new_k = (new_x, new_y)
        if (new_k not in c):
            if k == 0:
                cp[0] = new_k
            else:
                cp[-1] = new_k
            return (True, cp)
        
    return (False, cp)  # No move possible



def corner_move(c, k):
    """
    Applies a corner move to residue k, where k must be between 1 and n-2.
    Args:
        c (list of tuples): Current conformation.
        k (int): Index of the residue to move.
    Returns:
        tuple: (bool, new_conformation)
            bool: True if the move was successful, False otherwise.
            new_conformation: The new conformation after the move (or the original if no move was possible).
    """

    cp = c.copy()

    # Create variables containing the coordinates of k, its previous neighbors, and its next one.
    coord, coord_prev, coord_next = cp[k], cp[k-1], cp[k+1]
    x, y = coord[0], coord[1]
    x_prev, y_prev = coord_prev[0], coord_prev[1]
    x_next, y_next = coord_next[0], coord_next[1]
    if x_prev != x_next and y_prev != y_next:

        # The condition above checks if k and its neighbors form a corner.
        if x_prev == x and ((x_next, y_prev) not in c):
            cp[k] = (x_next, y_prev)
            return (True, cp)
        elif ((x_prev, y_next) not in c):
            cp[k] = (x_prev, y_next)
            return (True, cp)

    return (False, cp)



def crankshaft_move(c, k):
    """
    Applies a crankshaft move to residue k, where k must be between 1 and n-3.
    We consider k as the first corner of the U-shaped segment.
    Args:
        c (list of tuples): Current conformation as a list of (x, y) coordinates.
        k (int): Index of the residue to move (must be between 1 and n-3).
    Returns:
        tuple: (bool, new_conformation)
            bool: True if the move was successful, False otherwise.
            new_conformation: The new conformation after the move (or the original if no move was possible).
    """
    cp = c.copy()

    # Get coordinates of the four residues involved in the crankshaft move
    x_prev, y_prev = cp[k-1]   # Previous residue (first corner of U)
    x, y = cp[k]              # Current residue (second corner of U)
    x_next, y_next = cp[k+1]  # Next residue (third corner of U)
    x_next2, y_next2 = cp[k+2] # Residue after next (fourth corner of U)

    # Check for horizontal U-shape (residues k-1 and k+2 have same x-coordinate)
    if x_prev == x_next2 and x == x_next:
        
        # Case 1: U opens to the right - rotate 180° to the left
        if x == x_prev + 1 and ((x - 2, y) not in c) and ((x_next - 2, y_next) not in c):
            cp[k] = (x - 2, y)          # Move current residue left by 2 units
            cp[k+1] = (x_next - 2, y_next)  # Move next residue left by 2 units
            return (True, cp)
        
        # Case 2: U opens to the left - rotate 180° to the right
        elif x == x_prev - 1 and ((x + 2, y) not in c) and ((x_next + 2, y_next) not in c):
            cp[k] = (x + 2, y)          # Move current residue right by 2 units
            cp[k+1] = (x_next + 2, y_next)  # Move next residue right by 2 units
            return (True, cp)
        
    # Check for vertical U-shape (residues k-1 and k+2 have same y-coordinate)
    elif y_prev == y_next2 and y == y_next:

        # Case 3: U opens downward - rotate 180° upward
        if y == y_prev + 1 and ((x, y - 2) not in c) and ((x_next, y_next - 2) not in c):
            cp[k] = (x, y - 2)          # Move current residue down by 2 units
            cp[k+1] = (x_next, y_next - 2)  # Move next residue down by 2 units
            return (True, cp)
        
        # Case 4: U opens upward - rotate 180° downward
        elif y == y_prev - 1 and ((x, y + 2) not in c) and ((x_next, y_next + 2) not in c):
            cp[k] = (x, y + 2)          # Move current residue up by 2 units
            cp[k+1] = (x_next, y_next + 2)  # Move next residue up by 2 units
            return (True, cp)
        
    # If none of the above conditions are met, return False with original conformation
    return (False, cp)



def pull_move(c, k, max_try=3):
    """
    Applies a pull move (forward or backward) to residue k.
    Args:
        c (list of tuples): Current conformation as a list of (x, y) coordinates.
        k (int): Index of the residue to move (must be between 0 and n-3).
        max_try (int, optional): Maximum number of attempts to find a valid move. Defaults to 3.
    Returns:
        tuple: (bool, new_conformation)
            bool: True if the move was successful, False otherwise.
            new_conformation: The new conformation after the move (or the original if no move was possible).
    """
    if k <= len(c)-3:
        bool_forward, c_forward = pull_move_forward(c, k, max_try)
        if bool_forward:
            return bool_forward, c_forward
    return pull_move_backward(c, k, max_try)



def pull_move_backward(c, k, max_try=3):
    """
    Applies a pull move backward to residue k where k must be between 0 and n-3.
    Args:
        c (list of tuples): Current conformation as a list of (x, y) coordinates.
        k (int): Index of the residue to move (must be between 0 and n-3).
        max_try (int, optional): Maximum number of attempts to find a valid move. Defaults to 3.
    Returns:
        tuple: (bool, new_conformation)
            bool: True if the move was successful, False otherwise.
            new_conformation: The new conformation after the move (or the original if no move was possible).
    """
    cp = c.copy()
    n = len(c)
    cp.reverse()
    c0 = pull_move_forward(cp, n-k+1, max_try)
    c1 = c0[1]
    c1.reverse()
    return c0[0], c1



def pull_move_forward(c, k, max_try=3):
    """
    Applies a pull move forward to residue k where k must be between 0 and n-3.
    Args:
        c (list of tuples): Current conformation as a list of (x, y) coordinates.
        k (int): Index of the residue to move (must be between 0 and n-3).
        max_try (int, optional): Maximum number of attempts to find a valid move. Defaults to 3.
    Returns:
        tuple: (bool, new_conformation)
            bool: True if the move was successful, False otherwise.
            new_conformation: The new conformation after the move (or the original if no move was possible).
    """
    cp = c.copy()
    cp2 = c.copy()
    n = len(cp)
    if k > n-3:
        return (False, c)
    q = 0
    seuil = 0
    while k + q != n-3 or seuil > max_try:
        xi, yi = cp2[k + q]
        xi_next, yi_next = cp2[k + q + 1]
        xi_next2, yi_next2 = cp2[k + q + 2]

        # Two possible candidates to be adjacent to previous residue and in the corner of the current one (L)
        L1_x, L1_y = xi + (yi - yi_next), yi + (xi_next - xi)
        L2_x, L2_y = xi - (yi - yi_next), yi - (xi_next - xi)

        # And the two candidates to be in the corner (C)
        C1_x, C1_y = xi_next - (yi_next - yi), yi_next - (xi - xi_next)
        C2_x, C2_y = xi_next + (yi_next - yi), yi_next + (xi - xi_next)

        # Check if L1/2 is free and C1/2 is in the correct place with the next residue
        cond_L1 = (L1_x, L1_y) not in cp2
        cond_L2 = (L2_x, L2_y) not in cp2
        cond_C1_in_Cp2 = (C1_x, C1_y) not in cp2
        cond_C2_in_Cp2 = (C2_x, C2_y) not in cp2

        cp2 = cp.copy()
        rand = random.randint(1, 2)

        # If a candidate is free and the corresponding C is in the right place, make the move and return
        if cond_L1 and (C1_x, C1_y) == (xi_next2, yi_next2):
            cp[k + q + 1] = (L1_x, L1_y)
            if is_valid_conformation(cp):
                return (True, cp)
            elif seuil < max_try:
                cp = c.copy()
                cp2 = c.copy()
                q = 0
                seuil += 1
            else:
                return (False, c)

        elif cond_L2 and (C2_x, C2_y) == (xi_next2, yi_next2):
            cp[k + q + 1] = (L2_x, L2_y)
            if is_valid_conformation(cp):
                return (True, cp)
            elif seuil < max_try:
                cp = c.copy()
                cp2 = c.copy()
                q = 0
                seuil += 1
            else:
                return (False, c)
        elif rand == 1 and cond_L1 and cond_C1_in_Cp2:
            cp[k + q + 1] = (L1_x, L1_y)
        elif cond_L2 and cond_C2_in_Cp2:
            cp[k + q + 1] = (L2_x, L2_y)
        elif cond_L1 and cond_C1_in_Cp2:
            cp[k + q + 1] = (L1_x, L1_y)
        q += 1
    return (False, c)



# ----- Neighbourhoods Tests -----
if __name__ == "__main__":

    test = "test_pull_move"  # "test_crankshaft_move"  # "test_corner_move"  # "test_end_move"

    # ----- Test End Move -----
    if test == "test_end_move":
        hp = "HPHP"
        c = [(0, 0), (0, 1), (0, 2), (0, 3)]
        cp = end_move(c, 0)
        plot_molecules_side_by_side(c, cp[1], hp)

    # ----- Test Corner Move -----
    if test == "test_corner_move":
        hp = "HPHP"
        c = [(0, 0), (1, 0), (1, 1), (2, 1)]
        cp = corner_move(c, 2)
        plot_molecules_side_by_side(c, cp[1], hp)

    # ----- Test Crankshaft Move -----
    if test == "test_crankshaft_move":
        hp = "HPPHHPPHHP"
        c = [(2, -2), (2, -1), (2, 0), (2, 1), (1, 1), 
            (1, 2), (0, 2), (0, 1), (-1, 1), (-1, 0)]
        cp = crankshaft_move(c, 5)
        plot_molecules_side_by_side(c, cp[1], hp)

    # ----- Test Pull Move -----
    if test == "test_pull_move":
        hp = "HPHHPPHPPH"
        c = [(0, 0), (0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (2, 2), (2, 1), (2, 0)]
        # c = [(2, 0), (2, 1), (2, 2), (2, 3), (1, 3), (1, 2), (0, 2), (0, 1), (0, 0)]
        c = [(0,0), (0,1), (0,2), (1,2), (2,2), (3,2), (3,1), (2,1), (2,0), (2,-1)]
        cp = pull_move(c, 6)
        print("Pull move result 0:", pull_move(c, 0))
        print("Pull move result 1:", pull_move(c, 1))
        print("Pull move result 2:", pull_move(c, 2))
        print("Pull move result 3:", pull_move(c, 3))
        print("Pull move result 4:", pull_move(c, 4))
        print("Pull move result 5:", pull_move(c, 5))
        print("Pull move result 6:", pull_move(c, 6))
        print("Pull move result 7:", pull_move(c, 7))
        plot_molecules_side_by_side(c, cp[1], hp)
