"""
Flexible Sudoku Solver: GUI
"""

import tkinter as tk
from tkinter import messagebox, simpledialog, filedialog
import time
from collections import defaultdict
import math
import json
from solver import SudokuSolver

"""
This class paints the cells and extras and translates mouse clicks to cell coordinates.
A subclass has to be created for each type of Sudoku (Classic, Tridoku, Star...)
"""
class SudokuPainter:
  def __init__(self, canvas, rows, cols, cell_size=40, offset=(2,2)):
    self.CELLSIZE = cell_size
    self.OFFSET = offset
    self.rows = rows
    self.cols = cols
    self.canvas = canvas

  def xy2cell(self,x,y):
    pass
  def cell2xy(self,cell):
    pass
  def paint_cell(self, cell, value):
    pass
  def paint_extra(self, cells, fill="", outline="", text="", text_pos=(0,0)):
    pass
  def default_regions(self):
    pass


class ClassicPainter(SudokuPainter):
  def __init__(self, canvas, rows, cols, cell_size=40, offset=(2,2)):
    super().__init__(canvas, rows, cols, cell_size, offset)

  def default_regions(self):
    s = int(math.isqrt(self.cols))
    regions = []
    if self.cols == self.rows and s * s == self.cols:
      # square blocks
      for br in range(s):
        for bc in range(s):
          block = []
          for r in range(br * s, (br + 1) * s):
            for c in range(bc * s, (bc + 1) * s):
              block.append((r, c))
          regions.append(block)
    # rows, columns
    for r in range(self.rows):
      regions.append([(r, i) for i in range(self.cols)])
    for c in range(self.cols):
      regions.append([(i, c) for i in range(self.rows)])
    return regions

  def cell2xy(self, cell):
    return ((cell[1]+0.5)*self.CELLSIZE, (cell[0]+0.5)*self.CELLSIZE)

  def xy2cell(self,x,y):
    c = x // self.CELLSIZE
    r = y // self.CELLSIZE
    return (r, c)
    
  def paint_cell(self, cell, text="", fill="", outline = "", textcol="black", is_opt = False, textsize=10, width=None, dash=None, padding=0):
    x0 = cell[1] * self.CELLSIZE + self.OFFSET[0]
    y0 = cell[0] * self.CELLSIZE + self.OFFSET[1]
    x1 = x0 + self.CELLSIZE
    y1 = y0 + self.CELLSIZE
    padding *= self.CELLSIZE

    if fill != "" or outline != "":
      self.canvas.create_rectangle(x0+padding, y0+padding, x1-padding, y1-padding, fill=fill, outline=outline, width=width, dash=dash)

      # heavy square outlines (if perfect square sudoku)
      if self.rows == self.cols and padding == 0:
        s = math.isqrt(self.cols)
        if s*s == self.cols:
          if cell[0]%s == 0:
            self.canvas.create_line(x0, y0, x1,y0, fill=outline, width=2)
          if cell[1]%s == 0:
            self.canvas.create_line(x0, y0, x0,y1, fill=outline, width=2)

    if text!="":
      if not is_opt: print("draw text",text,"at",cell, textcol)
      if is_opt:
        s = math.isqrt(self.rows)
        xt = ((int(text)-1)%s + 0.5)/s
        yt = ((int(text)-1)//s +0.5)/s
        xt = x0 +xt * self.CELLSIZE
        yt = y0 +yt * self.CELLSIZE
      else:
        xt = x0 + self.CELLSIZE/2
        yt = y0 + self.CELLSIZE/2
      self.canvas.create_text(xt, yt, text=str(text), font=("Arial", textsize), fill=textcol)

  def paint_extra(self, cells, fill="", outline="", text="", text_pos=(0,0)):
    pass

class TriPainter(SudokuPainter):
  def __init__(self, canvas, rows, cols, cell_size, offset=(2,2)):
    super().__init__(canvas, rows, cols, cell_size, offset)
    self.w = self.CELLSIZE
    self.h = self.CELLSIZE*0.75**0.5

  def default_regions(self):
    regions = []
    # top triangle
    base1 = [(r, c) for r in range(3) for c in range(8-r, 9+r) ]
    # triangle below top
    base2 = [(r + 3, c) for r in range(3) for c in range(6+r, 11-r) ]

    dxy = [(0,0), (3,-3), (0,0), (3,3), (6,-6), (3,-3), (6,0), (3,3), (6,6)]
    for block in range(9):
      base = base2 if block in [2,5,7] else base1
      regions.append([(r + dxy[block][0], c+dxy[block][1]) for (r,c) in base])

    # sides:
    regions.append([(8,c) for c in range(0,18,2)]) # bottom row
    regions.append([(r,8-r) for r in range(9)]) # left side
    regions.append([(r,8+r) for r in range(9)]) # right side

    # inner lines:
    regions.append([(4,c) for c in range(4,13)]) # horizontal
    regions.append([(r,c) for r in range(4,9) for c in range(r-1,r+1) if (r,c) != (4,3)]) # left middle
    regions.append([(r,c) for r in range(4,9) for c in range(16-r,18-r) if (r,c) != (4, 13)]) # right middle
    return regions

  def xy2cell(self,x,y):
    row = int(y / self.h)
    tmpcol= int(x*2 / self.w)  # approximate column
    xy = [((c+1)*self.w/2, self.h * (row + (2-(row+c)%2)/3)) for c in range(tmpcol-1, tmpcol+2)]
    d2 = (2 * self.CELLSIZE)**2
    for c in range(-1, 2):
      dist2 = (x - xy[c][0])**2 + (y - xy[c][1])**2
      if dist2 < d2:
        d2 = dist2
        col = c+tmpcol
    return (row, col-1)    # why -1 ???

  def cell2xy(self,cell):
    r,c = cell
    upward = (r+c)%2 == 0
    x = (c+1)/2*self.w 
    y = (r + (2/3 if upward else 1/3))* self.h
    return (x,y)
    
  def paint_cell(self, cell, text="", fill="", outline = "",
      is_opt = False, textcol="black", textpos=(0.5, 0.5), textsize=10, width=1, dash=None, padding=0):
    r, c = cell
    x0 = c * self.w/2 + self.OFFSET[0]
    y0 = r * self.h + self.OFFSET[1]
    padding*=self.w
    upward = (r+c)%2 == 0
    if upward:
        points = [x0+padding, y0+self.h-padding, x0+self.w/2, y0+padding, x0+self.w-padding, y0+self.h-padding]
    else:
        points = [x0+padding, y0+padding, x0+self.w/2, y0+self.h-padding, x0+self.w-padding, y0+padding]
    self.canvas.create_polygon(points, fill=fill, outline=outline, width=width)

    if self.rows==9 and self.cols == 18 and padding ==0:
      if r%3 == 2 and upward:
        self.canvas.create_line(points[0], points[1], points[4], points[5], fill=outline, width=3*width)
      if (c-r)%3 == 2 and upward:
          self.canvas.create_line(points[2], points[3], points[4], points[5], fill=outline, width=3*width)
      if (r+c)%3 == 2 and upward:
          self.canvas.create_line(points[0], points[1], points[2], points[3], fill=outline, width=3*width)
    
    if text!="":
      if is_opt:
        s = 4
        o = int(text)
        o_r, o_c = (0,o-1)
        while o_c>(o_r if upward else 3 - o_r):
          o_r += 1
          o_c -= o_r if upward else 5 - o_r
          if o_r>3: break
        tx = o_c/s - o_r/2/s + 0.5 if upward else o_c/s + o_r/2/s
        ty = (o_r+0.5)/(s+1)
        tx = x0 + tx * self.CELLSIZE*0.8 + self.CELLSIZE*0.1 * (1 if upward else 2) 
        ty = y0 + ty * self.CELLSIZE*0.8 + self.CELLSIZE*0.2 * (1 if upward else 0.5)
      else:
        tx = sum(points[::2])/3
        ty = sum(points[1::2])/3
      self.canvas.create_text(tx, ty, text=str(text), font=("Arial", textsize), fill=textcol)

  def paint_extra(self, cells, fill="", outline="", text="", text_pos=(0,0)):
    pass

class StarPainter(TriPainter):
  def __init__(self, canvas, rows, cols, cell_size, offset=(2,2)):
    super().__init__(canvas, rows, cols, cell_size, offset)

  def default_regions(self):
    regions = []
    # upper left triangle
    base1 = [(r, c) for r in range(1,4) for c in range(1+r, 8-r) ]
    # middle left triangle
    base2 = [(r, c) for r in range(2,5) for c in range(4-r, 1+r) ]

    regions = [base1, base2]
    dxy = [(4,0), (2,6), (2,6), (-2,6)]
    for block in range(4):
      base = base2 if block in [1,3] else base1
      regions.append([(r + dxy[block][0], c+dxy[block][1]) for (r,c) in base])

    # horizontal lines:
    for r in range(1,7):
      line = set()
      for c in range(13):
        for reg in regions[0:6]:
          if (r,c) in reg: line.add((r,c))
      if r==1: line.add((0,8))
      if r==6: line.add((7,4))
      regions.append(line)

    # lines from top left to bottom right:
    for c in range(-2,9,2):
      line = set()
      for r in range(8):
        for reg in regions[0:6]:
          if (r,c+r) in reg: line.add((r,c+r))
          if (r,c+r-1) in reg: line.add((r,c+r-1))
      if c==-2: line.add((4,0))
      if c==8: line.add((3,12))
      regions.append(line)

    # lines from bottom left to top right:
    for c in range(4,15,2):
      line = set()
      for r in range(8):
        for reg in regions[0:6]:
          if (r,c-r) in reg: line.add((r,c-r))
          if (r,c-r+1) in reg: line.add((r,c-r+1))
      if c==4: line.add((1,2))
      if c==14: line.add((6,10))
      regions.append(line)

    return regions

"""
GUI class.
It contains the configuration buttons and the panel for drawing.
"""
class SudokuGUI:
    CELL_SIZE = 60
    OFFSET = (3, 3)

    def __init__(self, root):
        self.root = root
        root.title("Sudoku Solver")
        self.N = 9
        self.variant = "Classic"

        frame = tk.Frame(root)
        frame.pack(side="left", padx=10, pady=10)

        self.variant_var = tk.StringVar(value = self.variant)
        options = ["Classic", "Tridoku", "Star"]
        tk.OptionMenu(frame, self.variant_var, *options).pack()
        self.variant_var.trace_add("write", self.on_variant_change)

        tk.Label(frame, text="Grid size N (NxN):").pack()
        self.size_var = tk.IntVar(value=self.N)
        tk.Entry(frame, textvariable=self.size_var, width=4).pack()
        tk.Button(frame, text="Create Grid", command=self.create_grid).pack(pady=4)
        tk.Frame(frame, height=2, bd=1, relief="sunken").pack(fill="x", pady=4)
        tk.Button(frame, text="Solve", command=self.solve, bg="darkgreen", fg="white").pack(pady=8)
        tk.Label(frame, text='Time limit (s):').pack()
        self.time_limit_var = tk.DoubleVar(value=10.0)
        tk.Entry(frame, textvariable=self.time_limit_var, width=6).pack(pady=2)
        tk.Frame(frame, height=2, bd=1, relief="sunken").pack(fill="x", pady=4)

        tk.Button(frame, text="Clear Fixed", command=self.clear_fixed).pack()
        tk.Button(frame, text="Reset Extras", command=self.reset_constraints).pack(pady=6)

        tk.Frame(frame, height=2, bd=1, relief="sunken").pack(fill="x", pady=4)
        self.mode = "Normal"
        self.mode_var = tk.StringVar(value = self.mode)
        options = ["Normal", "Edit Regions", "Thermometer", "Sum", "Str8ts","(De)activate Cells"]
        tk.OptionMenu(frame, self.mode_var, *options).pack()
        self.mode_var.trace_add("write", self.on_mode_change)

        self.reg_btn_frame = tk.Frame(frame)
        tk.Button(self.reg_btn_frame, text="⟨", width=2, command=self.prev_region).pack(side="left", padx=1)
        tk.Button(self.reg_btn_frame, text="⟩", width=2, command=self.next_region).pack(side="left", padx=1)
        tk.Button(self.reg_btn_frame, text="+", width=2, command=self.add_region).pack(side="left", padx=1)
        
        self.reg_placeholder = tk.Frame(frame, height=1)
        self.reg_placeholder.pack(pady=2)
        self.reg_btn_frame.pack(in_=self.reg_placeholder)
        self.reg_btn_frame.lower() 

        tk.Frame(frame, height=2, bd=1, relief="sunken").pack(fill="x", pady=4)

        tk.Button(frame, text="Show Current Constraints", command=self.show_constraints).pack(pady=8)
        tk.Frame(frame, height=2, bd=1, relief="sunken").pack(fill="x", pady=4)

        btn_frame = tk.Frame(frame)
        btn_frame.pack(pady=8)
        tk.Button(btn_frame, text="Save", command=self.save).pack(side="left", padx=4)
        tk.Button(btn_frame, text="Load", command=self.load).pack(side="left", padx=4)

        self.status = tk.Label(root, text="Ready", bd=1, relief=tk.SUNKEN, anchor=tk.W)
        self.status.pack(side=tk.BOTTOM, fill=tk.X)

        self.canvas = tk.Canvas(root, width=self.N * self.CELL_SIZE + 4, height=self.N * self.CELL_SIZE + 4)
        self.canvas.pack(side="right", padx=10, pady=10)
        self.canvas.bind("<Button-1>", self.on_left_click)
        self.canvas.bind("<Button-3>", self.on_right_click)
        self.canvas.bind("<Key>", self.on_key)
        self.canvas.focus_set()

        self.active = set()
        self.assigned = []
        self.regions = [] # will be initialized in create_grid()
        self.extras = []  # unified extras
        self.selected_extra = None
        self.current_extra_build = None
        self.current_region_index = None
        self.mode = "Normal"
        self.last_selected = None
        self.errors = []

        self.create_grid()
        self.reset_constraints()

    def get_solver(self):
        fixed = getattr(self, 'fixed', None)
        active = self.active
        if len(active) == 0: active = None
        return SudokuSolver(self.N, self.regions, active_cells=active, extras=self.extras, fixed=fixed,
          check_func = self.tri_check if self.variant=="Tridoku" else None)

    def on_variant_change(self, *args):
        self.variant = self.variant_var.get()
        self.regions = self.default_regions()
        self.create_grid()
        self.reset_constraints()
        self.draw_grid()

    def default_regions(self):
      return self.painter.default_regions()

    def create_grid(self):
        # Creating a new grid resets regions and extras (explicit user action)

        self.N = int(self.size_var.get())

        if self.variant == "Classic":
            self.CELL_SIZE = 60
            self.painter = ClassicPainter(self.canvas, self.N, self.N, self.CELL_SIZE)
        elif self.variant == "Tridoku":
            self.CELL_SIZE = 80
            self.painter = TriPainter(self.canvas, self.N, 2*self.N, self.CELL_SIZE)
        elif self.variant == "Star":
            self.CELL_SIZE = 80
            self.painter = StarPainter(self.canvas, 7, 8, self.CELL_SIZE)

        self.canvas.config(width=self.N * self.CELL_SIZE + 4, height=self.N * self.CELL_SIZE + 4)

        self.regions = self.default_regions()
        self.extras = []
        self.solver = self.get_solver()
        self.active = self.solver.cells
        self.fixed = {}
        self.draw_grid()

    def draw_grid(self):
        self.canvas.delete("all")
        self.solver = self.get_solver()
        options = self.solver.options
        s = math.ceil(math.sqrt(self.N)-0.0001)

        # draw selected region
        if self.mode == "Edit Regions" and self.current_region_index is not None:
            highlight = set(self.regions[self.current_region_index])
            for (r, c) in highlight:
                if self.variant=="Classic":
                    x0 = c * self.CELL_SIZE + self.OFFSET[0]
                    y0 = r * self.CELL_SIZE + self.OFFSET[1]
                    x1 = x0 + self.CELL_SIZE
                    y1 = y0 + self.CELL_SIZE
                    self.canvas.create_rectangle(x0, y0, x1, y1, fill="#fa0", outline="#f80", width=3)
                elif self.variant=="Tridoku" or self.variant == "Star":
                    tri_h = self.CELL_SIZE*0.75**0.5
                    tri_w = self.CELL_SIZE
                    upward = (r+c)%2==0
                    x0 = c*tri_w/2 + self.OFFSET[0]
                    y0 = r*tri_h + self.OFFSET[1]
                    if upward:
                        points = [x0, y0+tri_h, x0+tri_w/2, y0, x0+tri_w, y0+tri_h]
                    else:
                        points = [x0, y0, x0+tri_w/2, y0+tri_h, x0+tri_w, y0]
                    self.canvas.create_polygon(points, fill="#fa0", outline="#f80", width=3)

        # draw thermometers
        for e in self.extras:
            if e["type"] != "thermometer":
                continue
            t = e["cells"]
            col = "#f88" if e==self.selected_extra else "#aaa"
            for i, (r, c) in enumerate(t):
                #x0 = (c + 0.5) * self.CELL_SIZE + self.OFFSET[0]
                #y0 = (r + 0.5) * self.CELL_SIZE + self.OFFSET[1]
                x0, y0 = self.painter.cell2xy((r,c))
                x0 += self.OFFSET[0]
                y0 += self.OFFSET[0]
                r0 = self.CELL_SIZE * 0.2
                r1 = r0*0.5
                if i == 0:
                    self.canvas.create_oval((x0 - r0, y0 - r0, x0 + r0, y0 + r0), fill=col, outline=col)
                else:
                    # draw link from previous point to this
                    prev = t[i - 1]
                    px, py = self.painter.cell2xy(prev)
                    px += self.OFFSET[0]
                    py += self.OFFSET[0]
                    self.canvas.create_line(px, py, x0, y0, fill=col, width=r0*0.6)
                    self.canvas.create_oval((x0 - r1, y0 - r1, x0 + r1, y0 + r1), fill=col, outline=col)

        # draw cells + pencilmarks
        for cell in self.solver.cells: # cell outline and background
          color = "#444" if cell not in self.active else ""
          if any(e[0][0]==cell or e[0][1]==cell for e in self.errors):
            color = "#600" if cell not in self.active else "#f88"
          self.painter.paint_cell(cell,
            text = "",
            fill = color,
            outline = "#bbb")
          if cell in self.fixed or cell in self.assigned: # fixed cell values
            color = "black"
            if cell in self.assigned and cell not in self.fixed: color = "#44f"
            if cell not in self.active: color = "white"
            val = self.fixed[cell] if cell in self.fixed else self.assigned[cell] if cell in self.assigned else ""
            self.painter.paint_cell(cell,
              text = str(val),
              textsize  = 18,
              textcol = color)
          elif cell in self.active: # display options
            if cell not in self.solver.options: print("CELL ",cell,"IS MISSING")
            opt = self.solver.options[cell]
            for o in opt:
              self.painter.paint_cell(cell,
                text = str(o),
                textsize  = 10,
                outline = "",
                is_opt = True,
                textcol = "#888")

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
                x0, y0 = self.painter.cell2xy((r,c))
                x0 += self.OFFSET[0]
                y0 += self.OFFSET[1]
                
                if self.variant == "Classic":
                  x0 -= 0.45 * self.CELL_SIZE
                  y0 -= 0.45 * self.CELL_SIZE
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
                elif self.variant in ["Tridoku", "Star"]:
                  w, h = (self.painter.w, self.painter.h)
                  upward = (r+c)%2 == 0
                  x0 -= 0.5 * w
                  y0 += h/3 if upward else -h/3
                  points = [x0, y0, x0+w, y0, x0+w/2, y0-h if upward else y0+h]
                  if (r, c-1) not in cells:
                      lines.append([points[0], points[1], points[4], points[5]])
                  if (r + (1 if upward else -1), c) not in cells:
                      lines.append([points[0], points[1], points[2], points[3]])
                  if (r, c + 1) not in cells:
                      lines.append([points[2], points[3], points[4], points[5]])

            for l in lines:
                self.canvas.create_line(l, fill=col, width=2, dash=(4, 2))

            x0, y0 = self.painter.cell2xy(minrc)
            x0 += self.OFFSET[0]
            y0 += self.OFFSET[1]
            upward = (minrc[0]+minrc[1])%2 == 0

            if self.variant == "Classic":
              x0 -= 0.25 * self.CELL_SIZE
              y0 -= 0.25 * self.CELL_SIZE
            else:
              x0 -= 0 if upward else 0.2 * self.painter.w
              y0 -= (0.25 if upward else 0.15) * self.painter.h
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
                self.painter.paint_cell((r,c), text="", fill="" , outline = col, width=2, padding=0.05)
                #x0 = (c + 0.05) * self.CELL_SIZE + self.OFFSET[0]
                #y0 = (r + 0.05) * self.CELL_SIZE + self.OFFSET[1]
                #x1 = x0 + self.CELL_SIZE * 0.9
                #y1 = y0 + self.CELL_SIZE * 0.9
                #self.canvas.create_rectangle(x0, y0, x1, y1, outline=col, width=3)

    def clear_fixed(self):
        self.fixed = {}
        self.draw_grid()

    def reset_constraints(self):
        self.extras = []
        self.active = {c for c in self.solver.cells}
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
            "variant": self.variant,
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
        self.variant = data.get("variant", "Classic")
        self.variant_var.set(self.variant) # create_grid gets called
        self.size_var.set(self.N)
        # restore regions (convert lists to tuples)
        if "regions" in data and len(data["regions"])>0:
            self.regions = [[(c[0], c[1]) for c in region] for region in data.get("regions", [])]
        else:
            self.regions = self.default_regions()
        # restore fixed
        fixed = data.get("fixed", [])
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

        self.solver = self.get_solver()
        self.active = set(self.solver.cells) - set(tuple(x) for x in inactive_cells)
        # resize canvas and redraw - DO NOT call create_grid() because it would clear extras/fixed
        self.canvas.config(width=self.N * self.CELL_SIZE + 4, height=self.N * self.CELL_SIZE + 4)
        self.draw_grid()

    @staticmethod
    def tri_check(cell, val, assignment):
      r,c = cell
      for dr in range (-1,2):
        if (r+c)%2==0:
          dcmax = 1 if dr==-1 else 2
        else:
          dcmax = 1 if dr==1 else 2
        for dc in range(-dcmax,dcmax+1):
          if (r+dr,c+dc) in assignment and val == assignment[(r+dr,c+dc)]: return False
      return True


    def solve(self):
        solver = self.get_solver()
        self.errors = solver.check()
        if len(self.errors) > 0:
            self.draw_grid()
            messagebox.showinfo("Error", "\n".join(map(lambda x: x[1], self.errors)))
            return
        t0 = time.time()
        time_limit = self.time_limit_var.get()
        sol = solver.solve(time_limit=time_limit)
        t1 = time.time()
        if not sol:
            messagebox.showinfo("Result", f'Could not find a solution (time {t1-t0:.2f}s)')
            return
        self.assigned = sol
        self.draw_grid()
        self.status.configure(text=f'Found solution in {t1-t0:.2f}s')

    def get_pos(self, x,y):
        x-=self.OFFSET[0]
        y-=self.OFFSET[1]
        return self.painter.xy2cell(x,y)

    def on_left_click(self, event):
        self.canvas.focus_set()
        row, col = self.get_pos(event.x, event.y)

        if (row,col) not in self.solver.cells:
            return

        if self.mode == "Normal":
            self.last_selected = (row, col)
            self.selected_extra = None
            for e in self.extras:
                if (row, col) in e["cells"]:
                    self.selected_extra = e  # pick first matching extra
                    break

            self.last_selected = (row, col)
            self.draw_grid()

            # highlight selected cell
            if self.variant == "Classic":
              x0 = col * self.CELL_SIZE + self.OFFSET[0]
              y0 = row * self.CELL_SIZE + self.OFFSET[1]
              x1 = x0 + self.CELL_SIZE
              y1 = y0 + self.CELL_SIZE
              self.canvas.create_rectangle(x0, y0, x1, y1, outline="red", width=2)
            elif self.variant == "Tridoku" or self.variant == "Star":
              tri_h = self.CELL_SIZE*0.75**0.5
              tri_w = self.CELL_SIZE
              x0 = col*tri_w/2 + self.OFFSET[0]
              y0 = row*tri_h + self.OFFSET[1]
              if (row+col)%2==0:
                  points = [x0, y0+tri_h, x0+tri_w/2, y0, x0+tri_w, y0+tri_h]
              else:
                  points = [x0, y0, x0+tri_w/2, y0+tri_h, x0+tri_w, y0]
              self.canvas.create_polygon(points, fill="white", outline="red", width=2)
                
            if self.selected_extra:
                self.status.config(text=f"Selected {self.selected_extra['type']} for editing")
            else:
                self.status.config(text="")
            return

        if self.mode == "Edit Regions":
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

        if self.mode == "(De)activate Cells":
            cell = (row, col)
            if cell in self.active:
                self.active.remove(cell)
            else:
                self.active.add(cell)
            self.draw_grid()
            return

        if self.mode == "Thermometer":
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

        if self.mode == "Sum":
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
        if self.mode == "Str8ts":
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
            self.mode_var.set("Normal")
            #self.current_extra_build = None
            #self.current_region_index = None
            #self.status.config(text="")
            #self.draw_grid()
            return 
        if self.mode == "Normal":
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
        elif self.mode == "Edit Regions":
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
        row, col = self.get_pos(event.x, event.y)

        if (row,col) not in self.solver.cells:
            return

        if self.mode=="Edit Regions":
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

    def on_mode_change(self, *args  ):
        self.mode = self.mode_var.get()
        self.selected_extra = None
        self.current_extra_build = None
        self.status.config(text="")
        if self.mode != "Edit Regions":
          self.current_region_index = None
        self.reg_btn_frame.lower()
        if self.mode == "Thermometer":
            self.status.config(text="Click thermometer cells from lowest to highest.\nTo finish: Click one of the cells again.")
        elif self.mode == "Sum":
            self.status.config(text="Click sum cells.\nTo finish: Click one of the cells again.")
        elif self.mode == "Str8ts":
            self.status.config(text="Click cells for Str8ts area.\nTo finish: Click one of the cells again.")
        elif self.mode == "(De)activate Cells":
            self.status.config(text="Click cells to toggle active/inactive.")
        elif self.mode == "Edit Regions":
            self.reg_btn_frame.lift()
            if not self.regions:
                self.regions = [[]]
            if self.current_region_index is None:
                self.current_region_index = 0
            self.update_region_label()
            self.status.config(text="Klicke Zellen, um sie zur Region hinzuzufügen oder (mit Rechtsklick) zu entfernen.")
        self.draw_grid()
