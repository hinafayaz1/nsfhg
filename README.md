import os
import tkinter as tk
import warnings
from tkinter import Tk, filedialog, messagebox, Menu, PhotoImage
from tkinter.ttk import Button, Label, Style, Frame, Progressbar, Separator

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Necessary for 3D plotting
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.SeqUtils import ProtParam

# Suppress specific warning
warnings.simplefilter('ignore', PDBConstructionWarning)


class ProteinStructureViewer:
    def __init__(self, master):
        self.master = master
        self.master.title("Protein Structure Viewer")
        self.master.geometry("800x600")

        self.canvas = tk.Canvas(self.master, width=800, height=600)
        self.canvas.pack(fill='both', expand=True)

        # Load background image
        bg_image_path = "path_to_your_image.png"  # Update with your image path
        if os.path.exists(bg_image_path):
            self.bg_image = PhotoImage(file=bg_image_path)
            self.canvas.create_image(0, 0, anchor='nw', image=self.bg_image)
        else:
            messagebox.showerror("Error", f"Image file not found: {bg_image_path}")

        self.style = Style()
        self.style.configure('TButton', font=('Arial', 12), foreground='black', background='#4CAF50', padx=10, pady=5)
        self.style.configure('TLabel', font=('Arial', 16), background='#f0f8ff', padx=10, pady=5)
        self.style.configure('Horizontal.TProgressbar', background='#4CAF50')
        self.style.configure('TFrame', background='#f0f8ff')

        self.create_menu()
        self.create_widgets()

    def create_menu(self):
        menubar = Menu(self.master)
        self.master.config(menu=menubar)

        file_menu = Menu(menubar, tearoff=0)
        file_menu.add_command(label="Open PDB File", command=self.browse_file)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.master.quit)
        menubar.add_cascade(label="File", menu=file_menu)

        help_menu = Menu(menubar, tearoff=0)
        help_menu.add_command(label="About", command=self.show_about)
        menubar.add_cascade(label="Help", menu=help_menu)

    def create_widgets(self):
        frame = Frame(self.master, padding="20")
        frame_window = self.canvas.create_window(400, 300, window=frame, anchor='center')

        title_label = Label(frame, text="Protein Structure Viewer", style='TLabel')
        title_label.pack(pady=10)

        browse_button = Button(frame, text="Browse PDB File", style='TButton', command=self.browse_file)
        browse_button.pack(pady=10)

        Separator(frame, orient='horizontal').pack(fill='x', padx=10, pady=10)

        self.progress_bar_load = Progressbar(frame, orient="horizontal", mode="indeterminate",
                                             style='Horizontal.TProgressbar')
        self.progress_bar_load.pack(pady=10, fill='x')

        self.status_label = Label(frame, text="", style='TLabel', foreground='green', anchor='w')
        self.status_label.pack(side='bottom', fill='x', padx=10, pady=5)

    def browse_file(self):
        filename = filedialog.askopenfilename(filetypes=[("PDB files", ".pdb"), ("All files", ".*")])
        if filename:
            self.pdb_file_path = filename
            self.progress_bar_load.start()
            self.master.after(100, self.visualize_protein)

    def visualize_protein(self):
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("protein", self.pdb_file_path)

            coords = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            coords.append(atom.coord)

            coords = list(zip(*coords))  # Transpose the coordinates
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(coords[0], coords[1], coords[2], color='blue')

            ax.set_xlabel('X Axis')
            ax.set_ylabel('Y Axis')
            ax.set_zlabel('Z Axis')
            ax.set_title('Protein Structure')
            ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio

            plt.show()
            self.predict_structure(structure)
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")
        finally:
            self.progress_bar_load.stop()

    def predict_structure(self, structure):
        try:
            seq = self.get_sequence_from_structure(structure)
            if seq:
                prot_param = ProtParam.ProteinAnalysis(seq)
                helix, turn, sheet = prot_param.secondary_structure_fraction()
                messagebox.showinfo("Secondary Structure Prediction",
                                    f"Predicted Secondary Structure:\nHelix: {helix:.2f}, Turn: {turn:.2f}, Sheet: {sheet:.2f}")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred while predicting secondary structure: {str(e)}")

    def get_sequence_from_structure(self, structure):
        seq = ""
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[0] == ' ':
                        seq += residue.get_resname().strip()
        return seq

    def show_about(self):
        messagebox.showinfo("About",
                            "Protein Structure Viewer\nVersion 1.0\n\nA simple tool to visualize protein structures from PDB files.")


if __name__ == "__main__":
    root = Tk()
    app = ProteinStructureViewer(root)
    root.mainloop()
# nsfhg
