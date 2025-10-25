"""
Flexible Sudoku Solver
"""

from gui import SudokuGUI
import tkinter as tk

if __name__ == "__main__":
    root = tk.Tk()
    app = SudokuGUI(root)
    root.mainloop()
