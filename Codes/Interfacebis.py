from Monte_Carlo import *
from Grid import *
from Neighbourhoods import *
from Others_function import *
import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import sys
from math import exp



def plot_molecule_interface(c, hp_sequence, point_size=200, grid_color='gray', bg_color='white'):
    fig, ax = plt.subplots(figsize=(6, 6), facecolor=bg_color)
    ax.set_aspect('equal')
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, color=grid_color)
    min_x, min_y = min(x for x, y in c) - 1, min(y for x, y in c) - 1
    max_x, max_y = max(x for x, y in c) + 1, max(y for x, y in c) + 1
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    for i, (x, y) in enumerate(c):
        color = 'red' if hp_sequence[i] == 'H' else 'blue'
        ax.scatter(x, y, color=color, s=point_size)
        ax.text(x, y, str(i), ha='center', va='center', color='white', fontsize=10)
    for i in range(len(c) - 1):
        x1, y1 = c[i]
        x2, y2 = c[i+1]
        ax.plot([x1, x2], [y1, y2], 'k-', linewidth=1.5)
    ax.legend(handles=[
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='H (Hydrophobic)'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='P (Polar)'),
        plt.Line2D([0], [0], color='black', linewidth=1.5, label='Bond')
    ], loc='upper right')
    plt.title("HP Molecule")
    return fig



# --- Interface Graphique ---
class HPModelApp:
    def __init__(self, root):
        self.root = root
        self.root.title("HP Model Simulation")
        self.root.geometry("1000x600")
        # Variables for parameters
        self.method_var = tk.StringVar(value="MCsearch")
        self.phi = tk.IntVar(value=1000)
        self.nu = tk.DoubleVar(value=0.5)
        self.T = tk.DoubleVar(value=160.0)
        self.E_star = tk.IntVar(value=-5)
        self.T_init = tk.DoubleVar(value=160.0)
        self.T_final = tk.DoubleVar(value=220.0)
        self.chi = tk.IntVar(value=5)
        self.hp_sequence = tk.StringVar(value="HPPHHPH")

        ##################
        self.progress = None
        self.progress_value = tk.DoubleVar(value=0.0)
        ##################

        # Left frame for parameters
        self.left_frame = ttk.Frame(root, padding="10")
        self.left_frame.pack(side=tk.LEFT, fill=tk.Y)

        # Right frame for results
        self.right_frame = ttk.Frame(root, padding="10")
        self.right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # Widgets for MCsearch
        self.mc_widgets = []
        # Widgets for REMCSimulation
        self.remc_widgets = []

        # Configure the window close button
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

        self.setup_left_panel()
        self.setup_right_panel()

        # Update parameter fields when method changes
        self.method_var.trace_add("write", self.update_parameter_fields)

    def on_closing(self):
        """Handle the window closing event."""
        self.root.quit()
        sys.exit(0)

    def setup_left_panel(self):
        # Method selection
        ttk.Label(self.left_frame, text="Method:", font=('Arial', 10, 'bold')).grid(row=0, column=0, sticky=tk.W)
        ttk.Radiobutton(self.left_frame, text="MCsearch", variable=self.method_var, value="MCsearch").grid(row=1, column=0, sticky=tk.W)
        ttk.Radiobutton(self.left_frame, text="REMCSimulation", variable=self.method_var, value="REMCSimulation").grid(row=2, column=0, sticky=tk.W)

        # HP sequence
        ttk.Label(self.left_frame, text="HP Sequence:", font=('Arial', 10, 'bold')).grid(row=3, column=0, sticky=tk.W)
        ttk.Entry(self.left_frame, textvariable=self.hp_sequence, width=30).grid(row=4, column=0, sticky=tk.W)

        # MCsearch parameters
        self.mc_widgets.append(ttk.Label(self.left_frame, text="Number of iterations (phi):", font=('Arial', 10, 'bold')))
        self.mc_widgets.append(ttk.Entry(self.left_frame, textvariable=self.phi, width=10))
        self.mc_widgets.append(ttk.Label(self.left_frame, text="Pull move probability (nu):", font=('Arial', 10, 'bold')))
        self.mc_widgets.append(ttk.Entry(self.left_frame, textvariable=self.nu, width=10))
        self.mc_widgets.append(ttk.Label(self.left_frame, text="Temperature (T):", font=('Arial', 10, 'bold')))
        self.mc_widgets.append(ttk.Entry(self.left_frame, textvariable=self.T, width=10))

        # REMCSimulation parameters
        self.remc_widgets.append(ttk.Label(self.left_frame, text="Target energy (E_star):", font=('Arial', 10, 'bold')))
        self.remc_widgets.append(ttk.Entry(self.left_frame, textvariable=self.E_star, width=10))
        self.remc_widgets.append(ttk.Label(self.left_frame, text="Initial temperature (T_init):", font=('Arial', 10, 'bold')))
        self.remc_widgets.append(ttk.Entry(self.left_frame, textvariable=self.T_init, width=10))
        self.remc_widgets.append(ttk.Label(self.left_frame, text="Final temperature (T_final):", font=('Arial', 10, 'bold')))
        self.remc_widgets.append(ttk.Entry(self.left_frame, textvariable=self.T_final, width=10))
        self.remc_widgets.append(ttk.Label(self.left_frame, text="Number of replicas (chi):", font=('Arial', 10, 'bold')))
        self.remc_widgets.append(ttk.Entry(self.left_frame, textvariable=self.chi, width=10))

        # Run button
        ttk.Button(self.left_frame, text="Run Simulation", command=self.run_simulation).grid(row=15, column=0, pady=10, sticky=tk.W)

        # Quit button
        ttk.Button(self.left_frame, text="Quit", command=self.on_closing).grid(row=16, column=0, pady=5, sticky=tk.W)

        #############################
        # Ajoutez la barre de progression (masquée par défaut)
        self.progress = ttk.Progressbar(self.left_frame, orient="horizontal", length=200, mode="determinate", variable=self.progress_value)
        self.progress.grid(row=14, column=0, pady=5, sticky=tk.W)
        self.progress.grid_remove()  # Masquée par défaut
        #############################

        # Show MCsearch widgets by default
        self.update_parameter_fields()
        

    def update_parameter_fields(self, *args):
        # Hide all dynamic widgets
        for widget in self.mc_widgets + self.remc_widgets:
            widget.grid_forget()

        #################
        if self.method_var.get() == "REMCSimulation":
            # Afficher la barre de progression uniquement pour REMCSimulation
            self.progress.grid()
        else:
            self.progress.grid_remove()
        #################

        # Show widgets based on selected method
        if self.method_var.get() == "MCsearch":
            self.mc_widgets[0].grid(row=5, column=0, sticky=tk.W)
            self.mc_widgets[1].grid(row=6, column=0, sticky=tk.W)
            self.mc_widgets[2].grid(row=7, column=0, sticky=tk.W)
            self.mc_widgets[3].grid(row=8, column=0, sticky=tk.W)
            self.mc_widgets[4].grid(row=9, column=0, sticky=tk.W)
            self.mc_widgets[5].grid(row=10, column=0, sticky=tk.W)
        else:
            self.remc_widgets[0].grid(row=5, column=0, sticky=tk.W)
            self.remc_widgets[1].grid(row=6, column=0, sticky=tk.W)
            self.remc_widgets[2].grid(row=7, column=0, sticky=tk.W)
            self.remc_widgets[3].grid(row=8, column=0, sticky=tk.W)
            self.remc_widgets[4].grid(row=9, column=0, sticky=tk.W)
            self.remc_widgets[5].grid(row=10, column=0, sticky=tk.W)
            self.remc_widgets[6].grid(row=11, column=0, sticky=tk.W)
            self.remc_widgets[7].grid(row=12, column=0, sticky=tk.W)


    def setup_right_panel(self):
        self.result_label = ttk.Label(self.right_frame, text="Result:", font=('Arial', 12, 'bold'))
        self.result_label.pack(anchor=tk.NW)
        self.result_text = tk.Text(self.right_frame, height=5, width=50)
        self.result_text.pack(fill=tk.X, padx=5, pady=5)
        self.canvas_frame = ttk.Frame(self.right_frame)
        self.canvas_frame.pack(fill=tk.BOTH, expand=True)

    def generate_initial_coordinates(self, hp_sequence):
        """Generate initial coordinates as a straight line."""
        n = len(hp_sequence)
        return [(i, 0) for i in range(n)]

    def run_simulation(self):
        try:
            hp = self.hp_sequence.get()
            c = self.generate_initial_coordinates(hp)
            phi = self.phi.get()
            nu = self.nu.get()

            if self.method_var.get() == "MCsearch":
                T = self.T.get()
                best_c, best_E = MCsearch(hp = hp, c = [], phi = phi, nu = nu, T=T)
            else:
                E_star = self.E_star.get()
                T_init = self.T_init.get()
                T_final = self.T_final.get()
                chi = self.chi.get()

                # Réinitialiser la barre de progression
                self.progress_value.set(0.0)

                # Définir un callback pour mettre à jour la progression
                def update_progress(value):
                    self.progress_value.set(value)
                    self.root.update_idletasks()  # Mettre à jour l'interface

                best_c, best_E = REMCSimulation(hp = hp, E_star = E_star, c = [], phi = phi, nu = nu, 
                                                T_init=T_init, T_final=T_final, chi=chi, 
                                                progress_callback=update_progress)

            # Afficher le résultat
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, f"Minimum energy found: {best_E}\n")
            self.result_text.insert(tk.END, f"Best configuration: {best_c}\n")

            # Afficher le plot
            fig = plot_molecule_interface(best_c, hp)
            self.display_plot(fig)

        except Exception as e:
            tk.messagebox.showerror("Error", f"An error occurred: {e}")


    def display_plot(self, fig):
        # Clear previous canvas
        for widget in self.canvas_frame.winfo_children():
            widget.destroy()
        # Display new plot
        canvas = FigureCanvasTkAgg(fig, master=self.canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)



# --- Lancement de l'application ---
if __name__ == "__main__":
    root = tk.Tk()
    app = HPModelApp(root)
    root.mainloop()