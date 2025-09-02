from Grid import *
from Neighbourhoods import *
from Monte_Carlo import *


# ----- Test End Move -----
hp = "HPHP"
c = [(0, 0), (0, 1), (0, 2), (0, 3)]
#cp = end_move(c, 0)
#plot_molecules_side_by_side(c, cp[1], hp)



# ----- Test Corner Move -----
hp = "HPHP"
c = [(0, 0), (1, 0), (1, 1), (2, 1)]
#cp = corner_move(c, 2)
#plot_molecules_side_by_side(c, cp[1], hp)



# ----- Test Crankshaft Move -----
hp = "HPPHHPPHHP"
c = [(2, -2), (2, -1), (2, 0), (2, 1), (1, 1), 
    (1, 2), (0, 2), (0, 1), (-1, 1), (-1, 0)]
#cp = crankshaft_move(c, 5)
#plot_molecules_side_by_side(c, cp[1], hp)


# ----- Test Energie -----
hp_sequence = "HPPHHPHPPH"
c = [(0, 0), (1, 0), (1, 1), (0, 1), (0, 2),
    (1, 2), (1, 3), (0, 3), (0, 4), (1, 4)]
energy = E(c, hp_sequence)
#print(f"Énergie de la conformation : {energy}")
#plot_molecule(c, hp_sequence)


# ----- Test M_VHSD -----
hp = "HPPHHPPHHP"
c = [(2, -2), (2, -1), (2, 0), (2, 1), (1, 1), 
    (1, 2), (0, 2), (0, 1), (-1, 1), (-1, 0)]
cp = M_vshd(c, 4)
plot_molecules_side_by_side(c, cp, hp)



# Exemple
hp = "HPHHPPHPP"
c = [(0, 0), (0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (2, 2), (2, 1), (2, 0)]
c = [(2, 0), (2, 1), (2, 2), (2, 3), (1, 3), (1, 2), (0, 2), (0, 1),(0, 0) ]
#cp = pull_move(c, 1)
#print("Pull move result:", cp)
#plot_molecules_side_by_side(c, cp[1], hp)
