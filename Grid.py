import matplotlib.pyplot as plt

def plot_molecule(c, hp_sequence, point_size=200, grid_color='gray', bg_color='white'):
    fig, ax = plt.subplots(figsize=(8, 8), facecolor=bg_color)
    ax.set_aspect('equal')
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, color=grid_color)

    min_x, min_y = min(x for x, y in c) - 1, min(y for x, y in c) - 1
    max_x, max_y = max(x for x, y in c) + 1, max(y for x, y in c) + 1
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)

    for i, (x, y) in enumerate(c):
        residue_type = hp_sequence[i]
        color = 'red' if residue_type == 'H' else 'blue'
        ax.scatter(x, y, color=color, s=point_size)
        ax.text(x, y, str(i), ha='center', va='center', color='white', fontsize=10)

    for i in range(len(c) - 1):
        x1, y1 = c[i]
        x2, y2 = c[i+1]
        ax.plot([x1, x2], [y1, y2], 'k-', linewidth=1.5)

    ax.legend(handles=[
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='H'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='P'),
        plt.Line2D([0], [0], color='black', linewidth=1.5, label='Bond')
    ], loc='upper right')

    plt.title("HP Molecule")
    plt.show()




def plot_molecules_side_by_side(c1, c2, hp_sequence, point_size=200, grid_color='gray', bg_color='white',
                               title1="Configuration 1", title2="Configuration 2"):

    # Verification of the size of the two sequences
    if len(c1) != len(hp_sequence) or len(c2) != len(hp_sequence):
        raise ValueError("Les configurations et la séquence HP doivent avoir la même longueur.")

    # Creation of the figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8), facecolor=bg_color)
    fig.suptitle(f"Configuration comparison: {hp_sequence}", fontsize=14, y=1.02)

    def plot_single_molecule(ax, c, title):
        ax.set_aspect('equal')
        ax.grid(True, which='both', linestyle='--', linewidth=0.5, color=grid_color)

        # Adjust the size of the grid
        min_x, min_y = min(x for x, y in c) - 1, min(y for x, y in c) - 1
        max_x, max_y = max(x for x, y in c) + 1, max(y for x, y in c) + 1
        ax.set_xlim(min_x, max_x)
        ax.set_ylim(min_y, max_y)
        ax.set_title(title)

        # Plot residues
        for i, (x, y) in enumerate(c):
            residue_type = hp_sequence[i]
            color = 'red' if residue_type == 'H' else 'blue'
            ax.scatter(x, y, color=color, s=point_size)
            ax.text(x, y, str(i), ha='center', va='center', color='white', fontsize=10)

        # Create bond between residues
        for i in range(len(c) - 1):
            x1, y1 = c[i]
            x2, y2 = c[i+1]
            ax.plot([x1, x2], [y1, y2], 'k-', linewidth=1.5)

        # Legend
        handles = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='H (Hydrophobic)'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='P (Polar)'),
            plt.Line2D([0], [0], color='black', linewidth=1.5, label='Bond')
        ]
        ax.legend(handles=handles, loc='upper right')

    # Tracer les deux configurations
    plot_single_molecule(ax1, c1, title1)
    plot_single_molecule(ax2, c2, title2)

    plt.tight_layout()
    plt.show()



# ----- Test -----
# c = [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (4, 1), (3, 1), (2, 1), (1, 1), (0, 1),  
#    (0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 3), (3, 3), (2, 3), (1, 3), (0, 3)]
# hp_sequence = "HPPHHPHPPHHPHPPHHPHP"
# plot_molecule(c, hp_sequence, point_size=200, grid_color='gray', bg_color='white')