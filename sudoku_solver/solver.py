"""
Flexible Sudoku Solver
"""

import time
import math
from collections import defaultdict

class SudokuSolver:
    """
    N: number of values = number of rows and columns in classic sudoku
    regions: list of list of cells (standard sudoku: 9 rows + 9 columns + 9 blocks)
    active_cells: list of cells to be filled in
    extras: list of thermometers, of regions with given sums, and of "str8ts" e.g. {"type": "sum", "value": 15, "sum": [(0,0),(0,1),(0,2)]}
    fixed: dict of predefined cells with their values
    check_func: additional check to be performed whenever a new value should be inserted into the grid
    """
    def __init__(self, N, regions, active_cells=None, extras=None, fixed=None, check_func=None):
        self.N = N
        self.regions = regions
        if not self.regions:
          self.regions = self._default_regions()

        cells = set()
        for region in self.regions:
            for cell in region:
                cells.add(cell)
        self.cells = list(cells)
        self.cells.sort(key = lambda cell: cell[0]*100+cell[1])
        
        self.active = set(self.cells) if active_cells is None else set(active_cells)

        self.extras = extras or []  # unified list for thermometers and sums
        self.fixed = fixed or {}
        if self.regions is None or not self.regions:
            self.regions = self._default_regions()
        self.cell_to_regions = defaultdict(list)
        for idx, region in enumerate(self.regions):
            for cell in region:
                self.cell_to_regions[cell].append(idx)
        self.peers = {cell: set() for cell in self.active}
        for region in self.regions:
            for a in region:
                if a not in self.active: continue
                self.peers[a].update([b for b in region if b != a])
        self._calc_options()
        self.check_func = check_func

    """
    Solve the sudoku.
    time_limit: after this time (in seconds) the solving process will be stopped
    """
    def solve(self, time_limit=10.0):
        self.tmax = time.time() + time_limit
        self.step=0
        return self._solve_recursive(self.fixed.copy())


    def _solve_recursive(self, assignment):
        if time.time()>self.tmax: return None
        self.step+=1
        required_cells = {c for c in self.active}
        if all(c in assignment for c in required_cells):
            return assignment
        unassigned = [c for c in self.active if c not in assignment]
        # If any unassigned cell has no options => dead end
        for u in unassigned:
            if not self.options.get(u):
                return None
        cell = min(unassigned, key=lambda x: len(self.options[x]))
        for val in sorted(self.options[cell]):
            if not self._is_consistent(cell, val, assignment):
                continue
            assignment[cell] = val
            saved = []
            for p in self.peers[cell]:
                if p not in self.active: continue
                if val in self.options[p]:
                    saved.append(p)
                    self.options[p].discard(val)
            result = self._solve_recursive(assignment)
            if result:
                return result
            del assignment[cell]
            for p in saved:
                self.options[p].add(val)
        return None

    # calculate remaining number options per cell
    def _calc_options(self):
        self.options = {}
        self.cell_to_extras = { c: [] for c in self.cells }
        # lookup table for extras
        for extra in self.extras:
            for c in extra["cells"]:
                self.cell_to_extras[c].append(extra)
                
        # remove all numbers given in same group
        for c in self.active:
            if c in self.fixed:
                self.options[c] = set()
                continue
            self.options[c] = set(range(1, self.N + 1))
            for p in self.peers[c]:
                if p in self.fixed:
                    self.options[c].discard(self.fixed[p])

        # handle extras
        for c in self.active:
            for extra in self.cell_to_extras.get(c, []):
                if extra["type"] == "thermometer":
                    t = extra["cells"]
                    L = len(t)
                    for i, c in enumerate(t):
                        if c in self.fixed or c not in self.active:
                            continue
                        # each position i (0-based) must be at least (i+1) and at most N-(L-1-i)
                        min_allowed = i + 1
                        max_allowed = self.N - (L - 1 - i)
                        for val in list(self.options[c]):
                            if val < min_allowed or val > max_allowed:
                                self.options[c].discard(val)

                elif extra["type"] == "sum":
                    cells = extra["cells"]
                    L = len(cells)
                    val = extra.get("value", 0)
                    sum_fixed = 0
                    used = []
                    free_cells = []
                    for c in cells:
                        if c in self.fixed:
                            sum_fixed += self.fixed[c]
                            used.append(self.fixed[c])
                        else:
                            free_cells.append(c)
                    rest = val - sum_fixed
                    remaining = L - len(used)
                    if remaining <= 0:
                        continue
                    # compute min and max possible sums for remaining-1 cells (used for bounding single cell)
                    # build list of available numbers not used yet
                    available_numbers = [n for n in range(1, self.N + 1) if n not in used]
                    # For minsum and maxsum of (remaining-1) other cells:
                    if remaining - 1 <= 0:
                        minsum = 0
                        maxsum = 0
                    else:
                        minsum = sum(available_numbers[: remaining - 1])
                        maxsum = sum(available_numbers[-(remaining - 1) :])
                    for c in cells:
                        if c in self.fixed or c not in self.active:
                            continue
                        discard = []
                        for o in list(self.options[c]):
                            if o + maxsum < rest or o + minsum > rest:
                                discard.append(o)
                        for o in discard:
                            self.options[c].discard(o)
                elif extra["type"] == "str8ts":
                    cells = extra["cells"]
                    L = len(cells)
                    used = [self.fixed[c] for c in cells if c in self.fixed]
                    # check possible ranges
                    if c == cells[0]: # only once for each str8t
                        for _ in range(1):  # todo: repeat until no nore changes
                            start_values=[]
                            for start in range(1, self.N+1):
                                cells_ok = set()
                                for x in range(start, start+L):
                                    x_ok=False
                                    for c2 in cells:
                                        if (c2 in self.fixed and self.fixed[c2]==x) or x in self.options[c2]:
                                            cells_ok.add(c2)
                                            x_ok=True
                                    if not x_ok: break
                                if len(cells_ok)==L and x_ok: start_values.append(start)
                            valid=set()
                            for start in start_values:
                                for i in range(L): valid.add(start+i)
                            for c2 in cells:
                                self.options[c2] = {v for v in self.options[c2] if v in valid}
                    # additional check in case of some given numbers:
                    if len(used)>=1:
                        min_used = min(used)
                        max_used = max(used)
                        span = max_used - min_used + 1
                        # wenn der Bereich schon zu weit auseinander liegt, unlösbar
                        if span > L:
                            for c in cells:
                                if c in self.options:
                                    self.options[c].clear()
                    
                        # alle Zahlen müssen in einem Block von Länge L passen
                        min_possible = max(1, max_used - L + 1)
                        max_possible = min(self.N - L + 1, min_used)
                        # Zellen dürfen nur Werte aus diesem zusammenhängenden Block haben
                        allowed = set(range(min_possible, min_possible + L)) | set(range(max_possible, max_possible + L))
                        for c in cells:
                            if c not in self.fixed and c in self.active:
                                self.options[c].intersection_update(allowed)

    def _is_consistent(self, cell, val, assignment):
        # check thermometers and sums in extras
        for extra in self.cell_to_extras.get(cell, []):
            cells = extra["cells"]
            if extra["type"] == "thermometer":
                idx = cells.index(cell)
                if idx > 0 and cells[idx-1] in assignment and assignment[cells[idx-1]] >= val:
                    return False
                if idx < len(cells)-1 and cells[idx+1] in assignment and assignment[cells[idx+1]] <= val:
                    return False
            elif extra["type"] == "sum":
                dest = extra.get("value", 0)
                s_val = 0
                free = len(cells)
                for c in cells:
                    if c not in assignment:
                        continue
                    free -= 1
                    s_val += assignment[c]
                    if s_val > dest:
                        return False
                if s_val < dest - free * self.N:
                    return False
                # if this is the last free cell, check if sum ist correct
                if free == 1 and s_val != dest - val:
                    return False
            elif extra["type"] == "str8ts":
                vals = [assignment[c] for c in cells if c in assignment]
                vals.append(val)
                if len(vals) >= 2:
                    minv = min(vals)
                    maxv = max(vals)
                    span = maxv - minv + 1
                    if span > len(cells):
                        return False  # gap too large
                    # if all cells in str8t region are assigned, check if there is no gap:
                    if len(vals) == len(cells):
                        if len(set(vals))!=len(cells):
                            return False
        if self.check_func:
            if not self.check_func(cell, val, assignment):
                return False
        return True

    """
    Check for obvious inconsistencies (e.g. same number twice in the same ragion).
    Returns a list of errors (consisting of a pair of cells and a message).
    """
    def check(self):
        err = []
        for r in self.regions:
            for c1 in r:
                for c2 in r:
                    if c1[0] * 1000 + c1[1] >= c2[0] * 1000 + c2[1]:
                        continue
                    if (
                        c1 in self.fixed
                        and c2 in self.fixed
                        and self.fixed[c1] == self.fixed[c2]
                    ):
                        err.append([(c1, c2), f"{self.fixed[c1]} appears multiple times in a group."])
        return err
