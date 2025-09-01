import random

def is_adjacent(pos1, pos2):
    """Checks if two positions are adjacent on a 2D lattice."""
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1]) == 1


def end_move(c, k) :
    """
    Applies an end move to residue k, where k must be first (0) or last residue (n-1).
    Returns (True, new_conformation) if the move is possible, otherwise (False, unchanged_conformation).
    """
    cp = c.copy()
    if k == 0:
        neighbour_residue = c[1]  # Residue 1
    else:  # k == n-1
        neighbour_residue = c[-2]  # Residue n-2

    # Possible directions in a 2D lattice: up, down, left, right
    x_nr = neighbour_residue[0]
    y_nr = neighbour_residue[1]
    directions = [(0, 1), (0, -1), (1, 0), (-1, 0)] 
    random.shuffle(directions) # We shuffle to get a random position and not necessary the first one

    for x, y in directions:
        # We test all the neighbour position of neighbour_residue and return the first one which is empty
        new_x = x_nr + x
        new_y = y_nr + y
        new_k = (new_x, new_y)

        if  (new_k not in c) :
            if k == 0 :
                cp[0] = new_k
            else :
                cp[-1] = new_k
            return (True, cp)

    return (False, cp)  # No move possible




def corner_move(c, k):
    """
    Applies a corner move to residue k, where k must be between 1 and n-2.
    Returns (True, new_conformation) if the move is possible, otherwise (False, unchanged_conformation).
    """
    cp =  c.copy()

    # We create variable containing the coordonates of k, its previous neighbours and its next one.
    coord, coord_prev, coord_next = cp[k], cp[k-1], cp[k+1]
    x, y = coord[0], coord[1]
    x_prev, y_prev = coord_prev[0], coord_prev[1]
    x_next, y_next = coord_next[0], coord_next[1]

    if x_prev != x_next and y_prev != y_next : 
        # The condition above checks if k and its neighbours form a corner.
        if x_prev == x and ((x_next, y_prev) not in c) :
            cp[k] = (x_next, y_prev)
            return (True, cp)
        elif ((x_prev, y_next) not in c) :
            cp[k] = (x_prev, y_next)
            return (True, cp)
        
    return (False, cp)




def crankshaft_move(c, k) :
    """
    Applies a crankshaft move to residue k, where k must be between 1 and n-3. We consider k as the first corner of the U.
    Returns (True, new_conformation) if the move is possible, otherwise (False, unchanged_conformation).
    """

    cp =  c.copy()

    x_prev, y_prev = cp[k-1]
    x, y = cp[k]
    x_next, y_next = cp[k+1]
    x_next2, y_next2 = cp[k+2]

    if x_prev == x_next2 and x == x_next:
        if x == x_prev + 1 and ((x - 2, y) not in c) and ((x_next - 2, y_next) not in c) :
            cp[k] = (x - 2, y)
            cp[k+1] = (x_next - 2, y_next)
            return(True, cp)
        
        elif x == x_prev - 1 and ((x + 2, y) not in c) and ((x_next + 2, y_next) not in c) :
            cp[k] = (x + 2, y)
            cp[k+1] = (x_next + 2, y_next)
            return(True, cp)

    elif y_prev == y_next2 and y == y_next:
        if y == y_prev + 1 and ((x, y - 2) not in c) and ((x_next, y_next - 2) not in c) :
            cp[k] = (x, y - 2)
            cp[k+1] = (x_next, y_next - 2)
            return(True, cp)
        
        elif y == y_prev - 1 and ((x, y + 2) not in c) and ((x_next, y_next + 2) not in c) :
            cp[k] = (x, y + 2)
            cp[k+1] = (x_next, y_next + 2)
            return(True, cp)
        
    return(False, cp)















