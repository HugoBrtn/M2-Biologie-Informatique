import random

def Config_initiale(molécule) :
    """ On prend une molécule sous forme de chaine de caractère H ou P et on renvoit une configuration linéaire de cette 
    molécule sous frome de liste de tuples.

    Exemple : Config_initiale('HPH') retourne [(0, 0), (1, 0), (2, 0)]
    """
    c_initial = [(i, 0) for i in range(len(molécule))]
    return c_initial



"""
def E(c)








def MCsearch(phi, c, v) :

    n = len(c)

    for i in range (phi) :
        cp = c
        k = random.randint(1, n)
        cp = cp# M()
        DeltaE = 
"""