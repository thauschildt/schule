"""
Flexible Sudoku Solver
"""

import tkinter as tk
from tkinter import messagebox, simpledialog, filedialog
import time
from collections import defaultdict
import math
import json
from solver import SudokuSolver

class SudokuGUI:
    CELL_SIZE = 60
    OFFSET = (3, 3)

    def __init__(self, root):
        self.root = root
        root.title("Flexible Sudoku Solver")
        self.N = 9
        self.active = set()
        self.regions = [] # will be initialized in create_grid()
        self.extras = []  # unified extras
        self.selected_extra = None
        self.current_extra_build = None
        self.current_region_index = None
        self.mode = "normal"
        self.last_selected = None
        self.errors = []

        frame = tk.Frame(root)
        frame.pack(side="left", padx=10, pady=10)
        tk.Label(frame, text="Grid size N (NxN):").pack()
        self.size_var = tk.IntVar(value=self.N)
        tk.Entry(frame, textvariable=self.size_var, width=4).pack()
        tk.Button(frame, text="Create Grid", command=self.create_grid).pack(pady=4)
        tk.Frame(frame, height=2, bd=1, relief="sunken").pack(fill="x", pady=4)
        tk.Button(frame, text="Solve", command=self.solve).pack(pady=8)
        tk.Label(frame, text='Time limit (s):').pack()
        self.time_limit_var = tk.DoubleVar(value=10.0)
        tk.Entry(frame, textvariable=self.time_limit_var, width=6).pack(pady=2)
        tk.Frame(frame, height=2, bd=1, relief="sunken").pack(fill="x", pady=4)

        tk.Button(frame, text="Clear Fixed", command=self.clear_fixed).pack()
        tk.Button(frame, text="Reset Extras", command=self.reset_constraints).pack(pady=6)

        tk.Frame(frame, height=2, bd=1, relief="sunken").pack(fill="x", pady=4)
        self.region_button = tk.Button(frame, text="Edit Region", command=self.start_region)
        self.region_button.pack(pady=4)
        btn_frame = tk.Frame(frame)
        btn_frame.pack(pady=2)
        tk.Button(btn_frame, text="⟨", width=2, command=self.prev_region).pack(side="left", padx=1)
        tk.Button(btn_frame, text="⟩", width=2, command=self.next_region).pack(side="left", padx=1)
        tk.Button(btn_frame, text="+", width=2, command=self.add_region).pack(side="left", padx=1)

        tk.Frame(frame, height=2, bd=1, relief="sunken").pack(fill="x", pady=4)

        self.thermo_button = tk.Button(frame, text="Thermometer", command=self.start_thermo)
        self.thermo_button.pack(pady=6)

        self.sum_button = tk.Button(frame, text="Sums", command=self.start_sum)
        self.sum_button.pack(pady=6)

        self.str8ts_button = tk.Button(frame, text="Str8ts", command=self.start_str8ts)
        self.str8ts_button.pack(pady=6)
        self.active_button = tk.Button(frame, text="Deactivate Cells", command=self.toggle_deactivate_mode)
        self.active_button.pack(pady=6)
        tk.Frame(frame, height=2, bd=1, relief="sunken").pack(fill="x", pady=4)

        tk.Button(frame, text="Show Current Constraints", command=self.show_constraints).pack(pady=8)
        tk.Frame(frame, height=2, bd=1, relief="sunken").pack(fill="x", pady=4)
        tk.Button(frame, text="Save", command=self.save).pack(pady=8)
        tk.Button(frame, text="Load", command=self.load).pack(pady=8)
        tk.Frame(frame, height=2, bd=1, relief="sunken").pack(fill="x", pady=4)

        self.status = tk.Label(root, text="Bereit", bd=1, relief=tk.SUNKEN, anchor=tk.W)
        self.status.pack(side=tk.BOTTOM, fill=tk.X)

        self.canvas = tk.Canvas(root, width=self.N * self.CELL_SIZE + 4, height=self.N * self.CELL_SIZE + 4)
        self.canvas.pack(side="right", padx=10, pady=10)
        self.canvas.bind("<Button-1>", self.on_left_click)
        self.canvas.bind("<Button-3>", self.on_right_click)
        self.canvas.bind("<Key>", self.on_key)

        self.create_grid()
        self.canvas.focus_set()


    def create_grid(self):
        # Creating a new grid resets regions and extras (explicit user action)
        self.N = int(self.size_var.get())
        self.canvas.config(width=self.N * self.CELL_SIZE + 4, height=self.N * self.CELL_SIZE + 4)
        self.active = set((r, c) for r in range(self.N) for c in range(self.N))
        self.regions = SudokuSolver(self.N).regions
        self.extras = []
        self.fixed = {}
        self.draw_grid()

    def draw_grid(self):
        self.canvas.delete("all")
        solver_temp = SudokuSolver(self.N, active_cells=self.active, regions=self.regions or None, extras=self.extras, fixed=self.fixed)
        options = solver_temp.options
        s = math.ceil(math.sqrt(self.N)-0.0001)

        if self.mode == "region_edit" and self.current_region_index is not None:
            highlight = set(self.regions[self.current_region_index])
            for (r, c) in highlight:
                x0 = c * self.CELL_SIZE + self.OFFSET[0]
                y0 = r * self.CELL_SIZE + self.OFFSET[1]
                x1 = x0 + self.CELL_SIZE
                y1 = y0 + self.CELL_SIZE
                self.canvas.create_rectangle(x0, y0, x1, y1, fill="#fa0", outline="#f80", width=3)

        # draw thermometers
        for e in self.extras:
            if e["type"] != "thermometer":
                continue
            t = e["cells"]
            col = "#f88" if e==self.selected_extra else "#aaa"
            for i, (r, c) in enumerate(t):
                x0 = (c + 0.5) * self.CELL_SIZE + self.OFFSET[0]
                y0 = (r + 0.5) * self.CELL_SIZE + self.OFFSET[1]
                r0 = self.CELL_SIZE * 0.2
                r1 = r0*0.5
                if i == 0:
                    self.canvas.create_oval((x0 - r0, y0 - r0, x0 + r0, y0 + r0), fill=col, outline=col)
                else:
                    # draw link from previous point to this
                    prev = t[i - 1]
                    px = (prev[1] + 0.5) * self.CELL_SIZE + self.OFFSET[0]
                    py = (prev[0] + 0.5) * self.CELL_SIZE + self.OFFSET[1]
                    self.canvas.create_line(px, py, x0, y0, fill=col, width=r0*0.6)
                    self.canvas.create_oval((x0 - r1, y0 - r1, x0 + r1, y0 + r1), fill=col, outline=col)

        # draw cells + pencilmarks
        for r in range(self.N):
            for c in range(self.N):
                x0 = c * self.CELL_SIZE + self.OFFSET[0]
                y0 = r * self.CELL_SIZE + self.OFFSET[1]
                x1 = x0 + self.CELL_SIZE
                y1 = y0 + self.CELL_SIZE
                color = None if (r,c) in self.active else "#444"
                textcol = None if (r,c) in self.active else "white"
                for e in self.errors:
                    if e[0][0] == (r, c) or e[0][1] == (r, c):
                        color = "#f88"
                self.canvas.create_rectangle(x0, y0, x1, y1, fill=color, outline="#bbb", dash=(2, 2))
                if (r, c) in self.fixed:
                    self.canvas.create_text(x0 + self.CELL_SIZE / 2, y0 + self.CELL_SIZE / 2, text=str(self.fixed[(r, c)]), font=("Arial", 18), fill=textcol)
                elif (r,c) in self.active:
                    for i in range(1, self.N + 1):
                        cell_opts = options.get((r, c), set())
                        if i not in cell_opts:
                            continue
                        xx = (i - 1) % s + 0.5
                        yy = (i - 1) // s + 0.5
                        self.canvas.create_text(x0 + xx * self.CELL_SIZE / s, y0 + yy * self.CELL_SIZE / s, text=str(i), font=("Arial", 10), fill="#888")

        # heavy square outlines (if perfect square sudoku)
        if s * s == self.N:
            for br in range(s):
                for bc in range(s):
                    x0 = bc * self.CELL_SIZE * s + self.OFFSET[0]
                    y0 = br * self.CELL_SIZE * s + self.OFFSET[1]
                    x1 = x0 + self.CELL_SIZE * s
                    y1 = y0 + self.CELL_SIZE * s
                    self.canvas.create_rectangle(x0, y0, x1, y1, outline="#888", width=2)

        # draw sums (borders and clue)
        for e in self.extras:
            if e["type"] != "sum":
                continue
            col = "#f88" if e==self.selected_extra else "#080"
            cells = e["cells"]
            if not cells:
                continue
            lines = []
            minrc = (999, 999)
            for (r, c) in cells:
                if r < minrc[0] or (r == minrc[0] and c < minrc[1]):
                    minrc = (r, c)
                x0 = (c + 0.05) * self.CELL_SIZE + self.OFFSET[0]
                y0 = (r + 0.05) * self.CELL_SIZE + self.OFFSET[1]
                x1 = x0 + self.CELL_SIZE * 0.9
                y1 = y0 + self.CELL_SIZE * 0.9
                if (r - 1, c) not in cells:
                    lines.append([x0, y0, x1, y0])
                if (r + 1, c) not in cells:
                    lines.append([x0, y1, x1, y1])
                if (r, c - 1) not in cells:
                    lines.append([x0, y0, x0, y1])
                if (r, c + 1) not in cells:
                    lines.append([x1, y0, x1, y1])
            for l in lines:
                self.canvas.create_line(l, fill=col, width=2, dash=(4, 2))
            x0 = (minrc[1] + 0.25) * self.CELL_SIZE + self.OFFSET[0]
            y0 = (minrc[0] + 0.2) * self.CELL_SIZE + self.OFFSET[1]
            if e.get("value", 0) > 0:
                self.canvas.create_text(x0, y0, text=str(e["value"]), font=("Arial", 14), fill=col)

        ii=0
        for e in self.extras:
            if e["type"] != "str8ts": continue
            cells = e["cells"]
            cols=["#844","#484","#448","#884","#848","#884"]
            col=cols[ii]
            ii=(ii+1)%len(cols)
            for (r, c) in cells:
                x0 = (c + 0.05) * self.CELL_SIZE + self.OFFSET[0]
                y0 = (r + 0.05) * self.CELL_SIZE + self.OFFSET[1]
                x1 = x0 + self.CELL_SIZE * 0.9
                y1 = y0 + self.CELL_SIZE * 0.9
                self.canvas.create_rectangle(x0, y0, x1, y1, outline=col, width=3)

    def clear_fixed(self):
        self.fixed = {}
        self.draw_grid()

    def reset_constraints(self):
        self.extras = []
        self.active = [(r,c) for r in range(self.N) for c in range(self.N)]
        self.draw_grid()

    def show_constraints(self):
        region_count = len(self.regions)
        t_count = len([e for e in self.extras if e["type"] == "thermometer"])
        s_count = len([e for e in self.extras if e["type"] == "sum"])
        str_count = len([e for e in self.extras if e["type"] == "str8ts"])
        s = f"N = {self.N}\nActive cells = {len(self.active)}\nRegions = {region_count}\nThermometers = {t_count}\nSums = {s_count}\nStr8ts = {str_count}"
        messagebox.showinfo("Constraints", s)

    def save(self):
        filename = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("Sudoku files", "*.json"), ("All files", "*.*")],
            title="Save Sudoku"
        )
        if not filename:
            return
        fixed = [[list(k), self.fixed[k]] for k in self.fixed.keys()]  # store coords as lists
        # extras already serializable except tuples -> convert to lists
        extras_serial = []
        for e in self.extras:
            item = {"type": e["type"], "cells": [list(c) for c in e["cells"]]}
            if e["type"] == "sum":
                item["value"] = e.get("value", 0)
            extras_serial.append(item)
        data = {
            "N": self.N,
            "regions": [[[c[0], c[1]] for c in region] for region in self.regions],
            "fixed": fixed,
            "extras": extras_serial,
            "inactive": [[r, c] for r in range(self.N) for c in range(self.N) if (r,c) not in self.active]
        }
        with open(filename, "w") as f:
            json.dump(data, f)
        
    def load(self, filename=None):
        if not filename:
            filename = filedialog.askopenfilename(
                defaultextension=".json",
                filetypes=[("Sudoku files", "*.json"), ("All files", "*.*")],
                title="Load Sudoku"
            )
            if not filename:
                return
        with open(filename, "r") as f:
            data = json.load(f)
        self.N = data.get("N", self.N)
        self.size_var.set(self.N)
        # restore regions (convert lists to tuples)
        if "regions" in data and len(data["regions"])>0:
            self.regions = [[(c[0], c[1]) for c in region] for region in data.get("regions", [])]
        else:
            solver_temp = SudokuSolver(self.N, active_cells=self.active, regions=None)
            self.regions = solver_temp.regions
        # restore fixed
        fixed = data.get("fixed", [])
        if len(fixed)==0: fixed = data.get("givens", [])   # 
        self.fixed = {(g[0][0], g[0][1]): g[1] for g in fixed}
        # restore extras
        self.extras = []
        for e in data.get("extras", []):
            entry = {"type": e["type"], "cells": [(c[0], c[1]) for c in e["cells"]]}
            if e["type"] == "sum":
                entry["value"] = e.get("value", 0)
            self.extras.append(entry)
        # restore active/inactive cells
        inactive_cells = data.get("inactive", [])
        self.active = set((r, c) for r in range(self.N) for c in range(self.N)) - set(tuple(x) for x in inactive_cells)
        # resize canvas and redraw - DO NOT call create_grid() because it would clear extras/fixed
        self.canvas.config(width=self.N * self.CELL_SIZE + 4, height=self.N * self.CELL_SIZE + 4)
        self.draw_grid()

    def solve(self):
        solver = SudokuSolver(self.N, active_cells=self.active, regions=self.regions or None, extras=self.extras, fixed=self.fixed)
        self.errors = solver.check()
        if len(self.errors) > 0:
            self.draw_grid()
            messagebox.showinfo("Error", "\n".join(map(lambda x: x[1], self.errors)))
            return
        t0 = time.time()
        time_limit = self.time_limit_var.get()
        sol = solver.solve(max_solutions=1, time_limit=time_limit)
        t1 = time.time()
        if not sol:
            messagebox.showinfo("Result", f'Could not find a solution (time {t1-t0:.2f}s)')
            return
        self.fixed = sol
        self.draw_grid()
        self.status.configure(text=f'Found solution in {t1-t0:.2f}s')

    def on_left_click(self, event):
        self.canvas.focus_set()
        col = event.x // self.CELL_SIZE
        row = event.y // self.CELL_SIZE
        if not (0 <= row < self.N and 0 <= col < self.N):
            return

        if self.mode == "normal":
            self.last_selected = (row, col)
            self.selected_extra = None
            for e in self.extras:
                if (row, col) in e["cells"]:
                    self.selected_extra = e  # pick first matching extra
                    break
            self.draw_grid()

            self.last_selected = (row, col)
            self.draw_grid()
            x0 = col * self.CELL_SIZE + self.OFFSET[0]
            y0 = row * self.CELL_SIZE + self.OFFSET[1]
            x1 = x0 + self.CELL_SIZE
            y1 = y0 + self.CELL_SIZE
            self.canvas.create_rectangle(x0, y0, x1, y1, outline="red", width=2)
                
            if self.selected_extra:
                self.status.config(text=f"Selected {self.selected_extra['type']} for editing")
            else:
                self.status.config(text="")
            return

        if self.mode == "region_edit":
            if self.current_region_index is None or self.current_region_index >= len(self.regions):
                return
            region = self.regions[self.current_region_index]
            cell = (row, col)
            if cell in region:
                region.remove(cell)
            else:
                region.append(cell)
            self.update_region_label()
            self.draw_grid()
            return

        if self.mode == "deactivate":
            cell = (row, col)
            if cell in self.active:
                self.active.remove(cell)
            else:
                self.active.add(cell)
            self.draw_grid()
            return

        if self.mode == "thermo":
            # building a thermometer extra in last appended extras (or current builder)
            if self.current_extra_build is None or self.current_extra_build.get("type") != "thermometer":
                self.current_extra_build = {"type": "thermometer", "cells": []}
                self.extras.append(self.current_extra_build)
            t = self.current_extra_build["cells"]
            if (row, col) in t:
                # finish building if clicking an existing cell
                self.current_extra_build = None
            else:
                t.append((row, col))
                self.draw_grid()
            return

        if self.mode == "sum":
            if self.current_extra_build is None or self.current_extra_build.get("type") != "sum":
                self.current_extra_build = {"type": "sum", "cells": [], "value": 0}
                self.extras.append(self.current_extra_build)
            s = self.current_extra_build
            if (row, col) in s["cells"]:
                # prompt for sum value and finish
                n = len(s["cells"])
                N = self.N
                minval = n * (n + 1) // 2
                maxval = N * (N + 1) // 2 - (N - n) * (N - n + 1) // 2
                val = None
                while val is None:
                    ans = simpledialog.askstring("Summe", f"Gib den Summenwert ein ({minval}-{maxval}):")
                    if ans is None:
                        # user cancelled -> remove the sum and go back to normal
                        self.extras.pop()  # remove incomplete
                        self.mode = "normal"
                        self.current_extra_build = None
                        self.draw_grid()
                        return
                    try:
                        val_int = int(ans)
                    except ValueError:
                        continue
                    if minval <= val_int <= maxval:
                        val = val_int
                s["value"] = val
                self.current_extra_build = None
            else:
                s["cells"].append((row, col))
            self.draw_grid()
            return
        if self.mode == "str8ts":
            if self.current_extra_build is None or self.current_extra_build.get("type") != "str8ts":
                self.current_extra_build = {"type": "str8ts", "cells": []}
                self.extras.append(self.current_extra_build)
            s = self.current_extra_build
            if (row, col) in s["cells"]:
                self.current_extra_build = None
                self.status.config(text="")
            else:
                s["cells"].append((row, col))
            self.draw_grid()
            return

    def on_key(self, event):
        if event.keysym == "Escape":
            self.mode = "normal"
            self.current_extra_build = None
            self.current_region_index = None
            self.status.config(text="")
            self.reset_buttons()
            self.draw_grid()
            return 
        if self.mode == "normal":
            if self.last_selected is None:
                return
            cell = self.last_selected
            if event.char.isdigit():
                val = int(event.char)
                if 1 <= val <= self.N:
                    self.fixed[cell] = val
                    self.draw_grid()
            elif event.keysym in ["BackSpace", "Delete","Space"]:
                if self.selected_extra:
                    self.extras.remove(self.selected_extra)
                    self.draw_grid()
                if cell in self.fixed:
                    del self.fixed[cell]
                    self.draw_grid()
            elif event.keysym in ["Right", "Left","Up","Down"]:
                r, c = cell
                if event.keysym == "Right":
                    c = (c + 1) % self.N
                elif event.keysym == "Left":
                    c = (c - 1) % self.N
                elif event.keysym == "Up":
                    r = (r - 1) % self.N
                elif event.keysym == "Down":
                    r = (r + 1) % self.N
                self.last_selected = (r, c)
                self.draw_grid()
                x0 = c * self.CELL_SIZE + self.OFFSET[0]
                y0 = r * self.CELL_SIZE + self.OFFSET[1]
                x1 = x0 + self.CELL_SIZE
                y1 = y0 + self.CELL_SIZE
                self.canvas.create_rectangle(x0, y0, x1, y1, outline="red", width=2)
        elif self.mode == "region_edit":
            if event.keysym in ["BackSpace", "Delete"]:
                if self.current_region_index is not None and 0 <= self.current_region_index < len(self.regions):
                    del self.regions[self.current_region_index]
                    if len(self.regions) == 0:
                        self.current_region_index = None
                    else:
                        self.current_region_index = self.current_region_index % len(self.regions)
                    self.update_region_label()
                    self.draw_grid()
            elif event.keysym in ["Right", "Left","Up","Down"]:
                if self.current_region_index is None or len(self.regions) == 0:
                    return
                if event.keysym == "Right" or event.keysym == "Down":
                    self.current_region_index = (self.current_region_index + 1) % len(self.regions)
                else:
                    self.current_region_index = (self.current_region_index - 1) % len(self.regions)
                self.update_region_label()
                self.draw_grid()
            elif event.keysym == "plus":
                self.regions.append([])
                self.current_region_index = len(self.regions) - 1
                self.update_region_label()
                self.draw_grid()

    def on_right_click(self, event):

        # right click: if in normal mode and over an existing extra cell, remove that extra cell or extra
        col = event.x // self.CELL_SIZE
        row = event.y // self.CELL_SIZE
        if not (0 <= row < self.N and 0 <= col < self.N):
            return

        if self.mode=="region_edit":
            region = self.regions[self.current_region_index]
            if (row, col) in region:
                region.remove((row, col))
            # remove empty regions:
            self.regions = [r for r in self.regions if len(r) > 0]
            self.draw_grid()
            return

        e = self.selected_extra
        if not e: return
        # find extras containing this cell
        to_remove = False
        if (row, col) in e["cells"]:
            # remove the cell from the extra; if extra becomes empty, remove extra
            e["cells"].remove((row, col))
            if len(e["cells"]) == 0:
                to_remove = True
        if to_remove:
            self.extras.remove(e)
        self.draw_grid()

    def start_region(self):
        if self.mode == "region_edit":
            self.mode = "normal"
            self.current_region_index = None
            self.status.config(text="")
            self.region_button.config(relief="raised")
        else:
            if not self.regions:
                self.regions = [[]]
            if self.current_region_index is None:
                self.current_region_index = 0
            self.mode = "region_edit"
            self.update_region_label()
            self.status.config(text="Klicke Zellen, um sie zur Region hinzuzufügen oder zu entfernen.")
            self.reset_buttons()
            self.region_button.config(relief="sunken")
        self.draw_grid()

    def add_region(self):
        self.regions.append([])
        self.current_region_index = len(self.regions) - 1
        self.update_region_label()
        self.draw_grid()

    def prev_region(self):
        if self.current_region_index is None:
            self.current_region_index = 0
        else:
            self.current_region_index = (self.current_region_index - 1) % len(self.regions)
        self.update_region_label()
        self.draw_grid()

    def next_region(self):
        if self.current_region_index is None:
            self.current_region_index = 0
        else:
            self.current_region_index = (self.current_region_index + 1) % len(self.regions)
        self.update_region_label()
        self.draw_grid()
    
    def update_region_label(self):
        if self.current_region_index is None or not self.regions:
            self.status.config(text="Region: –")
        else:
            n = len(self.regions)
            k = self.current_region_index + 1
            size = len(self.regions[self.current_region_index])
            self.status.config(text=f"Region {k}/{n} ({size})")

    def start_thermo(self):
        # prepare to build a thermometer; each click adds the next cell, clicking an existing cell finishes
        self.selected_extra = None
        self.draw_grid()
        if self.mode == "thermo":
            # deactivate if already active
            self.mode = "normal"
            self.current_extra_build = None
            self.thermo_button.config(relief="raised")
            self.status.config(text="")
        else:
            self.mode = "thermo"
            self.current_extra_build = None
            self.reset_buttons()
            self.thermo_button.config(relief="sunken")
            self.status.config(text="Click thermometer cells from lowest to highest.\nTo finish: Click one of the cells again.")

    def start_sum(self):
        # prepare to build a sum group; each click adds cells, clicking an existing cell prompts for sum value and finishes
        self.selected_extra = None
        self.draw_grid()
        if self.mode == "sum":
            # deactivate if already active
            self.mode = "normal"
            self.current_extra_build = None
            self.sum_button.config(relief="raised")
            self.status.config(text="")
        else:
            self.mode = "sum"
            self.current_extra_build = None
            self.reset_buttons()
            self.sum_button.config(relief="sunken")
            self.status.config(text="Click sum cells.\nTo finish: Click one of the cells again.")

    def start_str8ts(self):
        self.selected_extra = None
        self.draw_grid()
        if self.mode == "str8ts":
            # deactivate if already active
            self.mode = "normal"
            self.current_extra_build = None           
            self.str8ts_button.config(relief="raised")
            self.status.config(text="")
        else:
            self.mode = "str8ts"
            self.current_extra_build = None
            self.status.config(text="Click cells for Str8ts area.\nTo finish: Click one of the cells again.")
            self.reset_buttons()
            self.str8ts_button.config(relief="sunken")

    def reset_buttons(self):
        self.active_button.config(relief="raised")
        self.sum_button.config(relief="raised")
        self.str8ts_button.config(relief="raised")
        self.active_button.config(relief="raised")
        self.thermo_button.config(relief="raised")
        self.region_button.config(relief="raised")
        self.current_extra_build = None
        if self.mode != "region_edit":
          self.current_region_index = None
        self.draw_grid()

    def toggle_deactivate_mode(self):
        if self.mode == "deactivate":
            self.mode = "normal"
            self.active_button.config(relief="raised")
            self.status.config(text="")
        else:
            self.mode = "deactivate"
            self.reset_buttons()
            self.active_button.config(relief="sunken")
            self.status.config(text="Click cells to toggle active/inactive.")
