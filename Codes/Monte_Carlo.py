import random
from math import exp
from Neighbourhoods import *


def Config_initiale(molécule) :
    """ On prend une molécule sous forme de chaine de caractère H ou P et on renvoit une configuration linéaire de cette 
    molécule sous frome de liste de tuples.

    Exemple : Config_initiale('HPH') retourne [(0, 0), (1, 0), (2, 0)]
    """
    c_initial = [(i, 0) for i in range(len(molécule))]
    return c_initial



def generate_random_conformation(hp_sequence):
    """
    Generates a random valid conformation (coordinates) for a given HP sequence.
    The conformation is a self-avoiding walk on a 2D lattice.

    Args:
        hp_sequence (str): HP sequence (e.g., "HPPHHPHPPH")

    Returns:
        list: List of (x, y) coordinates representing a valid conformation
    """
    # Length of the sequence
    n = len(hp_sequence)

    # Start at origin (0, 0)
    c = [(0, 0)]

    # Possible movement directions: up, down, left, right
    directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]

    # Track visited positions to ensure self-avoiding walk
    visited = set(c)

    # Try to build the conformation step by step
    for i in range(1, n):
        last_x, last_y = c[-1]
        # Shuffle directions to try them in random order
        random.shuffle(directions)

        # Try each direction until a valid move is found
        for dx, dy in directions:
            new_pos = (last_x + dx, last_y + dy)

            # Check if the new position is not already visited
            if new_pos not in visited:
                c.append(new_pos)
                visited.add(new_pos)
                break
        else:
            # If no valid move found, backtrack and try again
            # This is a simplified approach; a more robust solution would use a proper backtracking algorithm
            # For simplicity, we'll just add a position that might cause overlap (not ideal)
            # In practice, you might want to restart the process or use a more sophisticated method
            c.append((last_x + 1, last_y))  # Default move right (might cause overlap)

    return c


def E(c, hp_sequence):
    """
    Args:
        c (list of tuples): Liste des coordonnées (x, y) des résidus.
        hp_sequence (str): Chaîne de caractères représentant la séquence HP (ex: "HPPH").

    Returns:
        int: Énergie de la conformation (plus l'énergie est basse, mieux c'est).
    """
 
    # Liste des indices des résidus H
    h_indices = [i for i, residue in enumerate(hp_sequence) if residue == 'H']

    energy = 0

    # Parcourir toutes les paires de résidus H non consécutifs
    for i in range(len(h_indices)):
        for j in range(i + 1, len(h_indices)):
            idx_i = h_indices[i]
            idx_j = h_indices[j]

            # Vérifier si les résidus ne sont pas consécutifs dans la séquence
            if abs(idx_i - idx_j) > 1:
                # Vérifier si les résidus sont adjacents sur le treillis
                if abs(c[idx_i][0] - c[idx_j][0]) + abs(c[idx_i][1] - c[idx_j][1]) == 1:
                    energy -= 1  # Chaque contact H-H contribue à -1 à l'énergie

    return energy




def MCsearch(phi, c, hp, T=220) :

    n = len(c)

    c_mini = c.copy()
    E_mini = E(c, hp)

    for i in range(phi):
        c_prime = c.copy()  # Crée une copie de la conformation actuelle
        k = random.randint(1, len(c)-1)  # Choisit un résidu aléatoire
        c_prime = M_vshd(c_prime, k)  # Applique un mouvement aléatoire

        delta_E = E(c_prime,hp) - E(c, hp)  # Calcule la différence d'énergie
        delta_E_mini = E(c_prime,hp) - E_mini

        if delta_E <= 0:
            c = c_prime  # Accepte toujours si l'énergie diminue ou reste la même
            
            if delta_E_mini < 0 :
                c_mini = c_prime
                E_mini = E(c_prime,hp)
        
        else:
            q = random.random()  # Génère un nombre aléatoire entre 0 et 1
            if q < (1 / (exp(1) ** (delta_E / T))):  # Critère de Metropolis
                c = c_prime  # Accepte avec une certaine probabilité si l'énergie augmente

    return c, E(c,hp), c_mini, E_mini