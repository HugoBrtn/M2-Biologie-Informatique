from Monte_Carlo import *
from Neighbourhoods import *
from Others_function import *
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


def plot_molecule_interface(c, hp_sequence, point_size=200, grid_color='gray', bg_color='white'):
    """
    Plots an HP molecule configuration on a 2D grid for interface display.

    Args:
        c (list of tuples): List of (x, y) coordinates representing the molecule conformation.
        hp_sequence (str): HP sequence representing the molecule (e.g., "HPPHHPH").
        point_size (int, optional): Size of the points representing residues. Defaults to 200.
        grid_color (str, optional): Color of the grid lines. Defaults to 'gray'.
        bg_color (str, optional): Background color of the plot. Defaults to 'white'.

    Returns:
        matplotlib.figure.Figure: The figure object containing the plot.
    """

    # Create figure and axis with specified background color
    fig, ax = plt.subplots(figsize=(6, 6), facecolor=bg_color)
    ax.set_aspect('equal')

    # Set up grid with dashed lines
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, color=grid_color)

    # Automatically adjust grid limits based on molecule coordinates
    min_x, min_y = min(x for x, y in c) - 1, min(y for x, y in c) - 1
    max_x, max_y = max(x for x, y in c) + 1, max(y for x, y in c) + 1
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)

    # Plot each residue with appropriate color (red for H, blue for P)
    for i, (x, y) in enumerate(c):
        color = 'red' if hp_sequence[i] == 'H' else 'blue'
        ax.scatter(x, y, color=color, s=point_size)
        ax.text(x, y, str(i), ha='center', va='center', color='white', fontsize=10)

    # Create bonds between consecutive residues
    for i in range(len(c) - 1):
        x1, y1 = c[i]
        x2, y2 = c[i+1]
        ax.plot([x1, x2], [y1, y2], 'k-', linewidth=1.5)

    # Add legend to identify residue types and bonds
    ax.legend(handles=[
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='H (Hydrophobic)'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='P (Polar)'),
        plt.Line2D([0], [0], color='black', linewidth=1.5, label='Bond')
    ], loc='upper right')

    plt.title("HP Molecule")
    return fig


class HPModelApp:
    def __init__(self, root):
        self.root = root
        self.root.title("HP Model Simulation")
        self.root.protocol("WM_DELETE_WINDOW", self.root.quit)

        # Variables to store selected method and parameters
        self.method_var = tk.StringVar(value="Monte Carlo Search")
        self.hp_sequence = tk.StringVar(value="HPHPPHHPHPPHPHHPPHPH")
        self.E_star = tk.IntVar(value=-9)

        # Parameters for Monte Carlo Search
        self.mc_phi = tk.IntVar(value=500)
        self.mc_nu = tk.DoubleVar(value=0.5)
        self.mc_T = tk.DoubleVar(value=160.0)

        # Parameters for REMC multi Processes
        self.remc_phi = tk.IntVar(value=500)
        self.remc_nu = tk.DoubleVar(value=0.5)
        self.remc_T_init = tk.DoubleVar(value=160.0)
        self.remc_T_final = tk.DoubleVar(value=220.0)
        self.remc_chi = tk.IntVar(value=5)
        self.remc_max_iteration = tk.IntVar(value=300)
        self.remc_nb_processus = tk.IntVar(value=4)

        # Parameters for REMC parallelized
        self.remc_paral_phi = tk.IntVar(value=500)
        self.remc_paral_nu = tk.DoubleVar(value=0.5)
        self.remc_paral_T_init = tk.DoubleVar(value=160.0)
        self.remc_paral_T_final = tk.DoubleVar(value=220.0)
        self.remc_paral_chi = tk.IntVar(value=5)
        self.remc_paral_max_iterations = tk.IntVar(value=300)

        # GUI Layout
        self.setup_ui()

    def setup_ui(self):
        # Left Frame for Inputs
        left_frame = ttk.LabelFrame(self.root, text="Parameters", padding=10)
        left_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        # Method Selection
        ttk.Label(left_frame, text="Method:").grid(row=0, column=0, sticky="w")
        method_combobox = ttk.Combobox(
            left_frame,
            textvariable=self.method_var,
            values=["Monte Carlo Search", "REMC Multi Processes", "REMC Parallelized"],
            state="readonly"
        )
        method_combobox.grid(row=0, column=1, sticky="ew")
        method_combobox.bind("<<ComboboxSelected>>", self.on_method_change)

        # Common Parameters
        ttk.Label(left_frame, text="HP Sequence:").grid(row=1, column=0, sticky="w")
        ttk.Entry(left_frame, textvariable=self.hp_sequence).grid(row=1, column=1, sticky="ew")

        ttk.Label(left_frame, text="Target Energy (E*):").grid(row=2, column=0, sticky="w")
        ttk.Entry(left_frame, textvariable=self.E_star).grid(row=2, column=1, sticky="ew")

        # Monte Carlo Search Parameters Frame
        self.mc_frame = ttk.LabelFrame(left_frame, text="Monte Carlo Search Parameters", padding=5)
        self.mc_frame.grid(row=3, column=0, columnspan=2, sticky="ew", pady=5)
        ttk.Label(self.mc_frame, text="Iterations (phi):").grid(row=0, column=0, sticky="w")
        ttk.Entry(self.mc_frame, textvariable=self.mc_phi).grid(row=0, column=1, sticky="ew")
        ttk.Label(self.mc_frame, text="Pull Move Probability (nu):").grid(row=1, column=0, sticky="w")
        ttk.Entry(self.mc_frame, textvariable=self.mc_nu).grid(row=1, column=1, sticky="ew")
        ttk.Label(self.mc_frame, text="Temperature (T):").grid(row=2, column=0, sticky="w")
        ttk.Entry(self.mc_frame, textvariable=self.mc_T).grid(row=2, column=1, sticky="ew")

        # REMC multi Processes Parameters Frame
        self.remc_frame = ttk.LabelFrame(left_frame, text="REMC Simulation Parameters", padding=5)
        self.remc_frame.grid(row=4, column=0, columnspan=2, sticky="ew", pady=5)
        ttk.Label(self.remc_frame, text="MC Iterations (phi):").grid(row=0, column=0, sticky="w")
        ttk.Entry(self.remc_frame, textvariable=self.remc_phi).grid(row=0, column=1, sticky="ew")
        ttk.Label(self.remc_frame, text="Pull Move Probability (nu):").grid(row=1, column=0, sticky="w")
        ttk.Entry(self.remc_frame, textvariable=self.remc_nu).grid(row=1, column=1, sticky="ew")
        ttk.Label(self.remc_frame, text="Initial Temperature (T_init):").grid(row=2, column=0, sticky="w")
        ttk.Entry(self.remc_frame, textvariable=self.remc_T_init).grid(row=2, column=1, sticky="ew")
        ttk.Label(self.remc_frame, text="Final Temperature (T_final):").grid(row=3, column=0, sticky="w")
        ttk.Entry(self.remc_frame, textvariable=self.remc_T_final).grid(row=3, column=1, sticky="ew")
        ttk.Label(self.remc_frame, text="Number of Replicas (chi):").grid(row=4, column=0, sticky="w")
        ttk.Entry(self.remc_frame, textvariable=self.remc_chi).grid(row=4, column=1, sticky="ew")
        ttk.Label(self.remc_frame, text="Max REMC Iteration:").grid(row=5, column=0, sticky="w")
        ttk.Entry(self.remc_frame, textvariable=self.remc_max_iteration).grid(row=5, column=1, sticky="ew")
        ttk.Label(self.remc_frame, text="Number of Processes:").grid(row=6, column=0, sticky="w")
        ttk.Entry(self.remc_frame, textvariable=self.remc_nb_processus).grid(row=6, column=1, sticky="ew")

        # REMC parallelized Parameters Frame
        self.remc_paral_frame = ttk.LabelFrame(left_frame, text="REMC Simulation (Parallel) Parameters", padding=5)
        self.remc_paral_frame.grid(row=5, column=0, columnspan=2, sticky="ew", pady=5)
        ttk.Label(self.remc_paral_frame, text="MC Iterations (phi):").grid(row=0, column=0, sticky="w")
        ttk.Entry(self.remc_paral_frame, textvariable=self.remc_paral_phi).grid(row=0, column=1, sticky="ew")
        ttk.Label(self.remc_paral_frame, text="Pull Move Probability (nu):").grid(row=1, column=0, sticky="w")
        ttk.Entry(self.remc_paral_frame, textvariable=self.remc_paral_nu).grid(row=1, column=1, sticky="ew")
        ttk.Label(self.remc_paral_frame, text="Initial Temperature (T_init):").grid(row=2, column=0, sticky="w")
        ttk.Entry(self.remc_paral_frame, textvariable=self.remc_paral_T_init).grid(row=2, column=1, sticky="ew")
        ttk.Label(self.remc_paral_frame, text="Final Temperature (T_final):").grid(row=3, column=0, sticky="w")
        ttk.Entry(self.remc_paral_frame, textvariable=self.remc_paral_T_final).grid(row=3, column=1, sticky="ew")
        ttk.Label(self.remc_paral_frame, text="Number of Replicas (chi):").grid(row=4, column=0, sticky="w")
        ttk.Entry(self.remc_paral_frame, textvariable=self.remc_paral_chi).grid(row=4, column=1, sticky="ew")
        ttk.Label(self.remc_paral_frame, text="Max REMC Iterations:").grid(row=5, column=0, sticky="w")
        ttk.Entry(self.remc_paral_frame, textvariable=self.remc_paral_max_iterations).grid(row=5, column=1, sticky="ew")

        # Run Button
        ttk.Button(left_frame, text="Run Simulation", command=self.run_simulation).grid(row=8, column=0, columnspan=2, pady=10)

        # Quit Button
        ttk.Button(left_frame, text="Quit", command=self.root.quit).grid(row=9, column=0, columnspan=2, pady=10)

        # Right Frame for Results
        self.right_frame = ttk.LabelFrame(self.root, text="Results", padding=10)
        self.right_frame.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
        self.result_label = ttk.Label(self.right_frame, text="Minimum Energy: ")
        self.result_label.pack()

        # Frame for the plot
        self.plot_frame = ttk.Frame(self.right_frame)
        self.plot_frame.pack(fill=tk.BOTH, expand=True)
        self.fig, self.ax = plt.subplots(figsize=(6, 6), facecolor='white')
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack()

        # Frame for coordinates
        self.coords_frame = ttk.Frame(self.right_frame)
        self.coords_frame.pack(fill=tk.BOTH, expand=True)
        self.coords_label = ttk.Label(self.coords_frame, text="Coordinates of the best conformation:")
        self.coords_label.pack(anchor="w")
        self.coords_text = tk.Text(self.coords_frame, height=5, width=40, state='disabled')
        self.coords_text.pack(fill=tk.BOTH, expand=True)

        # Initially show only Monte Carlo parameters    
        self.on_method_change()

    def on_method_change(self, event=None):
        if self.method_var.get() == "Monte Carlo Search":
            self.mc_frame.grid()
            self.remc_frame.grid_remove()
            self.remc_paral_frame.grid_remove()
        elif self.method_var.get() == "REMC Multi Processes":
            self.mc_frame.grid_remove()
            self.remc_frame.grid()
            self.remc_paral_frame.grid_remove()
        else:  # Remc Parallelized
            self.mc_frame.grid_remove()
            self.remc_frame.grid_remove()
            self.remc_paral_frame.grid()

    def run_simulation(self):
        hp = self.hp_sequence.get()
        E_star = self.E_star.get()

        if self.method_var.get() == "Monte Carlo Search":
            phi = self.mc_phi.get()
            nu = self.mc_nu.get()
            T = self.mc_T.get()
            best_conformation, best_energy = MCsearch(hp, phi=phi, nu=nu, T=T)
        elif self.method_var.get() == "REMC Multi Processes":
            phi = self.remc_phi.get()
            nu = self.remc_nu.get()
            T_init = self.remc_T_init.get()
            T_final = self.remc_T_final.get()
            chi = self.remc_chi.get()
            max_iteration = self.remc_max_iteration.get()
            nb_processus = self.remc_nb_processus.get()
            best_conformation, best_energy = REMC_multi(
                hp, E_star, phi=phi, nu=nu, T_init=T_init, T_final=T_final, chi=chi, max_iteration=max_iteration, nb_processus=nb_processus
            )
        else:  # Remc Parallelized
            phi = self.remc_paral_phi.get()
            nu = self.remc_paral_nu.get()
            T_init = self.remc_paral_T_init.get()
            T_final = self.remc_paral_T_final.get()
            chi = self.remc_paral_chi.get()
            max_iterations = self.remc_paral_max_iterations.get()
            best_conformation, best_energy = REMC_paral(
                hp, E_star, phi=phi, nu=nu, T_init=T_init, T_final=T_final, chi=chi, max_iterations=max_iterations
            )

        # Update results
        self.result_label.config(text=f"Minimum Energy: {best_energy}")
        self.plot_conformation(best_conformation, hp)
        self.display_coordinates(best_conformation)

    def plot_conformation(self, conformation, hp_sequence):
        self.ax.clear()
        self.fig = plot_molecule_interface(conformation, hp_sequence)
        self.canvas.figure = self.fig
        self.canvas.draw()

    def display_coordinates(self, conformation):
        self.coords_text.config(state='normal')
        self.coords_text.delete(1.0, tk.END)
        self.coords_text.insert(tk.END, " , ".join([f"({x}, {y})" for x, y in conformation]))
        self.coords_text.config(state='disabled')


if __name__ == "__main__":
    root = tk.Tk()
    app = HPModelApp(root)
    root.mainloop()
