import random
import pandas as pd
from Others_function import *
from Grid import *



def M_vshd(c, k):
    """
    Applies a VSHD move (end, corner, or crankshaft) to residue k.
    Args:
        c (pd.DataFrame): Current conformation as a DataFrame with columns 'x' and 'y'
        k (int): Index of the residue to move (must be between 1 and n-3)
    Returns:
        tuple: (bool, new_conformation)
            bool: True if the move was successful, False otherwise
            new_conformation: The new conformation after the move (or the original if no move was possible)
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
        rand = random.randint(1, 2)  # Randomly choose between corner and crankshaft
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
        c (pd.DataFrame): Current conformation as a DataFrame with columns 'x' and 'y'
        k (int): Index of the residue to move (0 or n-1)
    Returns:
        tuple: (bool, new_conformation)
            bool: True if the move was successful, False otherwise
            new_conformation: The new conformation after the move (or the original if no move was possible)
    """
    cp = c.copy()
    if k == 0:
        neighbour_residue = c.iloc[1]  # Residue 1
    else:  # k == n-1
        neighbour_residue = c.iloc[-2]  # Residue n-2
    # Possible directions in a 2D lattice:
    x_nr, y_nr = neighbour_residue['x'], neighbour_residue['y']
    directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]
    random.shuffle(directions)  # We shuffle to get a random position and not necessarily the first one
    for dx, dy in directions:
        # We test all the neighbour positions of neighbour_residue and return the first one which is empty
        new_x, new_y = x_nr + dx, y_nr + dy
        new_k = (new_x, new_y)
        if not ((c['x'] == new_x) & (c['y'] == new_y)).any():
            if k == 0:
                cp.at[0, 'x'], cp.at[0, 'y'] = new_x, new_y
            else:
                cp.at[len(cp)-1, 'x'], cp.at[len(cp)-1, 'y'] = new_x, new_y
            return (True, cp)
    return (False, cp)  # No move possible



def corner_move(c, k):
    """
    Applies a corner move to residue k, where k must be between 1 and n-2.
    Args:
        c (pd.DataFrame): Current conformation as a DataFrame with columns 'x' and 'y'
        k (int): Index of the residue to move
    Returns:
        tuple: (bool, new_conformation)
            bool: True if the move was successful, False otherwise
            new_conformation: The new conformation after the move (or the original if no move was possible)
    """
    cp = c.copy()
    # We create variables containing the coordinates of k, its previous neighbours and its next one.
    coord = cp.iloc[k]
    coord_prev = cp.iloc[k-1]
    coord_next = cp.iloc[k+1]
    x, y = coord['x'], coord['y']
    x_prev, y_prev = coord_prev['x'], coord_prev['y']
    x_next, y_next = coord_next['x'], coord_next['y']
    if x_prev != x_next and y_prev != y_next:
        # The condition above checks if k and its neighbours form a corner.
        if x_prev == x and not ((c['x'] == x_next) & (c['y'] == y_prev)).any():
            cp.at[k, 'x'], cp.at[k, 'y'] = x_next, y_prev
            return (True, cp)
        elif not ((c['x'] == x_prev) & (c['y'] == y_next)).any():
            cp.at[k, 'x'], cp.at[k, 'y'] = x_prev, y_next
            return (True, cp)
    return (False, cp)



def crankshaft_move(c, k):
    """
    Applies a crankshaft move to residue k, where k must be between 1 and n-3.
    We consider k as the first corner of the U-shaped segment.
    Args:
        c (pd.DataFrame): Current conformation as a DataFrame with columns 'x' and 'y'
        k (int): Index of the residue to move
    Returns:
        tuple: (bool, new_conformation)
            bool: True if the move was successful, False otherwise
            new_conformation: The new conformation after the move (or the original if no move was possible)
    """
    cp = c.copy()
    # Get coordinates of the four residues involved in the crankshaft move
    x_prev, y_prev = cp.iloc[k-1]['x'], cp.iloc[k-1]['y']  # Previous residue (first corner of U)
    x, y = cp.iloc[k]['x'], cp.iloc[k]['y']              # Current residue (second corner of U)
    x_next, y_next = cp.iloc[k+1]['x'], cp.iloc[k+1]['y']  # Next residue (third corner of U)
    x_next2, y_next2 = cp.iloc[k+2]['x'], cp.iloc[k+2]['y']  # Residue after next (fourth corner of U)
    # Check for horizontal U-shape (residues k-1 and k+2 have same x-coordinate)
    if x_prev == x_next2 and x == x_next:
        # Case 1: U opens to the right - rotate 180° to the left
        if x == x_prev + 1 and not ((c['x'] == x - 2) & (c['y'] == y)).any() and not ((c['x'] == x_next - 2) & (c['y'] == y_next)).any():
            cp.at[k, 'x'], cp.at[k, 'y'] = x - 2, y          # Move current residue left by 2 units
            cp.at[k+1, 'x'], cp.at[k+1, 'y'] = x_next - 2, y_next  # Move next residue left by 2 units
            return (True, cp)
        # Case 2: U opens to the left - rotate 180° to the right
        elif x == x_prev - 1 and not ((c['x'] == x + 2) & (c['y'] == y)).any() and not ((c['x'] == x_next + 2) & (c['y'] == y_next)).any():
            cp.at[k, 'x'], cp.at[k, 'y'] = x + 2, y          # Move current residue right by 2 units
            cp.at[k+1, 'x'], cp.at[k+1, 'y'] = x_next + 2, y_next  # Move next residue right by 2 units
            return (True, cp)
    # Check for vertical U-shape (residues k-1 and k+2 have same y-coordinate)
    elif y_prev == y_next2 and y == y_next:
        # Case 3: U opens downward - rotate 180° upward
        if y == y_prev + 1 and not ((c['x'] == x) & (c['y'] == y - 2)).any() and not ((c['x'] == x_next) & (c['y'] == y_next - 2)).any():
            cp.at[k, 'x'], cp.at[k, 'y'] = x, y - 2          # Move current residue down by 2 units
            cp.at[k+1, 'x'], cp.at[k+1, 'y'] = x_next, y_next - 2  # Move next residue down by 2 units
            return (True, cp)
        # Case 4: U opens upward - rotate 180° downward
        elif y == y_prev - 1 and not ((c['x'] == x) & (c['y'] == y + 2)).any() and not ((c['x'] == x_next) & (c['y'] == y_next + 2)).any():
            cp.at[k, 'x'], cp.at[k, 'y'] = x, y + 2          # Move current residue up by 2 units
            cp.at[k+1, 'x'], cp.at[k+1, 'y'] = x_next, y_next + 2  # Move next residue up by 2 units
            return (True, cp)
    # If none of the above conditions are met, return False with original conformation
    return (False, cp)



def pull_move(c, k):
    """
    Applies a pull move to residue k where k must be between 0 and n-3.
    We consider k as the first corner of the U-shaped segment.
    Args:
        c (pd.DataFrame): Current conformation as a DataFrame with columns 'x' and 'y'
        k (int): Index of the residue to move
    Returns:
        tuple: (bool, new_conformation)
            bool: True if the move was successful, False otherwise
            new_conformation: The new conformation after the move (or the original if no move was possible)
    """
    cp = c.copy()
    cp2 = c.copy()
    n = len(cp)
    q = 0
    while k + q != n-3:
        xi, yi = cp2.iloc[k + q]['x'], cp2.iloc[k + q]['y']
        xi_next, yi_next = cp2.iloc[k + q + 1]['x'], cp2.iloc[k + q + 1]['y']
        xi_next2, yi_next2 = cp2.iloc[k + q + 2]['x'], cp2.iloc[k + q + 2]['y']
        # Two candidates possibles to be adjacent to previous residue and in the corner of the current one (L)
        L1_x, L1_y = xi + (yi - yi_next), yi + (xi_next - xi)
        L2_x, L2_y = xi - (yi - yi_next), yi - (xi_next - xi)
        # And the two candidates to be in the corner (C)
        C1_x, C1_y = xi_next - (yi_next - yi), yi_next - (xi - xi_next)
        C2_x, C2_y = xi_next + (yi_next - yi), yi_next + (xi - xi_next)
        # We check if L1/2 is free and C1/2 is on the same place with next residue
        cond_L1 = not ((cp2['x'] == L1_x) & (cp2['y'] == L1_y)).any()
        cond_L2 = not ((cp2['x'] == L2_x) & (cp2['y'] == L2_y)).any()
        cond_C1_in_Cp2 = not ((cp2['x'] == C1_x) & (cp2['y'] == C1_y)).any()
        cond_C2_in_Cp2 = not ((cp2['x'] == C2_x) & (cp2['y'] == C2_y)).any()
        cp2 = cp.copy()
        rand = random.randint(1, 2)
        # If a candidate is free and the corresponding C is in the right place, we make the move and return
        if cond_L1 and (C1_x, C1_y) == (xi_next2, yi_next2):
            cp.at[k + q + 1, 'x'], cp.at[k + q + 1, 'y'] = L1_x, L1_y
            if is_valid_conformation(cp):
                return (True, cp)
            else:
                return (False, c)
        elif cond_L2 and (C2_x, C2_y) == (xi_next2, yi_next2):
            cp.at[k + q + 1, 'x'], cp.at[k + q + 1, 'y'] = L2_x, L2_y
            if is_valid_conformation(cp):
                return (True, cp)
            else:
                return (False, c)
        elif rand == 1 and cond_L1 and cond_C1_in_Cp2:
            cp.at[k + q + 1, 'x'], cp.at[k + q + 1, 'y'] = L1_x, L1_y
        elif cond_L2 and cond_C2_in_Cp2:
            cp.at[k + q + 1, 'x'], cp.at[k + q + 1, 'y'] = L2_x, L2_y
        elif cond_L1 and cond_C1_in_Cp2:
            cp.at[k + q + 1, 'x'], cp.at[k + q + 1, 'y'] = L1_x, L1_y
        q += 1
    return (False, c)


# Exemple
hp = "HPHHPPHPPH"
c = pd.DataFrame({
    'x': [0, 0, 0, 1, 2, 2, 2, 2, 2, 2],
    'y': [0, 1, 2, 2, 2, 1, 0, -1, -2, -3]
})
cp = pull_move(c, 3)
print("Pull move result 0 :", pull_move(c, 0))
print("Pull move result 1 :", pull_move(c, 1))
print("Pull move result 2 :", pull_move(c, 2))
print("Pull move result 3 :", pull_move(c, 3))
print("Pull move result 4 :", pull_move(c, 4))
print("Pull move result 5 :", pull_move(c, 5))
print("Pull move result 6 :", pull_move(c, 6))
print("Pull move result 7 :", pull_move(c, 7))
plot_molecules_side_by_side(c, cp[1], hp)
