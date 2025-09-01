from Grid import *
from Neighbourhoods import *

# ----- Test End Move -----
hp = "HPHP"
c = [(0, 0), (0, 1), (0, 2), (0, 3)]
cp = end_move(c, 0)
#plot_molecules_side_by_side(c, cp[1], hp)



# ----- Test Corner Move-----
c = [(0, 0), (1, 0), (1, 1), (2, 1)]
cp = corner_move(c, 2)
plot_molecules_side_by_side(c, cp[1], hp)