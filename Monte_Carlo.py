import random
from math import exp


def Config_initiale(molécule) :
    """ On prend une molécule sous forme de chaine de caractère H ou P et on renvoit une configuration linéaire de cette 
    molécule sous frome de liste de tuples.

    Exemple : Config_initiale('HPH') retourne [(0, 0), (1, 0), (2, 0)]
    """
    c_initial = [(i, 0) for i in range(len(molécule))]
    return c_initial




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



"""
def MCsearch(phi, c, v) :

    n = len(c)

    for _ in range(phi):
        c_prime = c.copy()  # Crée une copie de la conformation actuelle
        k = random.randint(1, len(c))  # Choisit un résidu aléatoire
        c_prime = M(c_prime, k, nu)  # Applique un mouvement aléatoire

        delta_E = E(c_prime) - E(c)  # Calcule la différence d'énergie

        if delta_E <= 0:
            c = c_prime  # Accepte toujours si l'énergie diminue ou reste la même
        else:
            q = random.random()  # Génère un nombre aléatoire entre 0 et 1
            if q < (1 / (exp(1) ** (delta_E / T))):  # Critère de Metropolis
                c = c_prime  # Accepte avec une certaine probabilité si l'énergie augmente

    return c
"""