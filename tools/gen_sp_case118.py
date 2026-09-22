#!/usr/bin/env python3
"""Synthetic 118-bus case in the shape of the IEEE 118-bus system.

Same cardinalities as IEEE 118 (118 buses, 186 branches of which 9 are
transformers, 54 generators, 99 loads, 3 areas, an EHV backbone of 11
buses), own topology, own parameters, own names: no IEEE data is copied.
The script validates its own output with a Newton-Raphson power flow
including generator Q-limit enforcement, an island check, an N-1 count and
a double-outage solve, then writes a MATPOWER case file.
"""
import math
import random

import numpy as np

SEED = 118
rng = random.Random(SEED)

N = 118
EHV = [8, 9, 10, 26, 30, 38, 63, 64, 65, 68, 81]   # 345 kV buses (same count as IEEE)
EHV_SET = set(EHV)
AREAS = {1: range(1, 41), 2: range(41, 81), 3: range(81, 119)}   # 40, 40, 38 buses

# ---------------------------------------------------------------- layout --
# Each area is a rough 2-D cloud; the EHV buses sit on a wide corridor that
# crosses all three areas so that the backbone ties the areas together.
pos = {}
centres = {1: (0.0, 0.0), 2: (60.0, 20.0), 3: (120.0, 0.0)}
for a, buses in AREAS.items():
  cx, cy = centres[a]
  for b in buses:
    if b in EHV_SET:
      continue
    pos[b] = (cx + rng.uniform(-28, 28), cy + rng.uniform(-22, 22))
# backbone corridor, ordered west to east
corridor_x = np.linspace(-20, 140, len(EHV))
for x, b in zip(corridor_x, EHV):
  pos[b] = (float(x), 30.0 + rng.uniform(-6, 6))

def dist(a, b):
  (x1, y1), (x2, y2) = pos[a], pos[b]
  return math.hypot(x1 - x2, y1 - y2)

# --------------------------------------------------------------- branches --
branches = []   # (f, t, kind) kind in {"L138", "L345", "T"}
edge_set = set()

def add_edge(f, t, kind):
  key = (min(f, t), max(f, t))
  if key in edge_set or f == t:
    return False
  edge_set.add(key)
  branches.append((f, t, kind))
  return True

# 1) 138 kV: per area a minimum spanning tree over distance, then chords
#    between close pairs until the target line count per area is reached
def area_lines(buses, n_target):
  buses = [b for b in buses if b not in EHV_SET]
  # Prim's MST
  in_tree = {buses[0]}
  while len(in_tree) < len(buses):
    best = None
    for u in in_tree:
      for v in buses:
        if v in in_tree:
          continue
        d = dist(u, v)
        if best is None or d < best[0]:
          best = (d, u, v)
    _, u, v = best
    in_tree.add(v)
    add_edge(u, v, "L138")
  # chords: shortest pairs not yet connected, skipping the longest ones
  pairs = sorted(((dist(u, v), u, v) for i, u in enumerate(buses) for v in buses[i + 1:]), key=lambda p: p[0])
  count = sum(1 for (f, t, k) in branches if k == "L138" and f in buses)
  for d, u, v in pairs:
    if count >= n_target:
      break
    if (min(u, v), max(u, v)) in edge_set:
      continue
    if d > 26.0:
      continue
    add_edge(u, v, "L138")
    count += 1

area_lines(list(AREAS[1]), 56)
area_lines(list(AREAS[2]), 58)
area_lines(list(AREAS[3]), 51)

# 2) 345 kV backbone: a chain along the corridor plus a few chords (12 lines)
for i in range(len(EHV) - 1):
  add_edge(EHV[i], EHV[i + 1], "L345")
for (u, v) in [(8, 10), (26, 38), (63, 65)]:
  add_edge(u, v, "L345")   # 10 + 3 = 13; one chord dropped below to reach 12
branches = [b for b in branches if not (b[2] == "L345" and {b[0], b[1]} == {26, 38})]
edge_set.discard((26, 38))

# 3) transformers: 9, each EHV bus to its nearest 138 kV bus in the same area
#    (two EHV buses share a substation with a neighbour and get none)
trafo_ehv = [8, 9, 26, 30, 38, 63, 65, 68, 81]
for e in trafo_ehv:
  area = next(a for a, r in AREAS.items() if e in r)
  cands = sorted((dist(e, b), b) for b in AREAS[area] if b not in EHV_SET)
  for _, b in cands:
    if add_edge(e, b, "T"):
      break

# 4) inter-area 138 kV ties so the areas are meshed below the backbone too
for (u, v) in [(31, 41), (36, 44), (77, 82), (79, 84), (3, 42)]:
  add_edge(u, v, "L138")

plants = {8: 200.0, 10: 380.0, 26: 300.0, 30: 180.0, 38: 200.0, 63: 180.0, 65: 350.0, 68: 200.0, 81: 300.0,
          4: 60.0, 12: 85.0, 19: 60.0, 24: 40.0, 25: 120.0, 27: 60.0, 31: 30.0, 32: 60.0, 34: 40.0, 36: 40.0,
          42: 60.0, 46: 60.0, 49: 150.0, 54: 100.0, 55: 40.0, 59: 120.0, 61: 120.0, 62: 40.0, 66: 150.0,
          70: 40.0, 72: 40.0, 74: 30.0, 76: 40.0, 77: 60.0, 80: 150.0, 85: 40.0, 87: 20.0, 89: 150.0,
          90: 40.0, 92: 60.0, 99: 40.0, 100: 150.0, 103: 40.0, 104: 40.0, 105: 40.0, 107: 30.0,
          110: 40.0, 111: 36.0, 113: 40.0, 116: 80.0, 69: 0.0}   # 69: slack
for _b in list(plants):
  if _b in AREAS[1]:
    plants[_b] = round(plants[_b] * 0.92)
  elif _b in AREAS[2] and _b != 69:
    plants[_b] = round(plants[_b] * 0.70)
  elif _b in AREAS[3]:
    plants[_b] = round(plants[_b] * 0.88)
# the north-west pocket (buses 1..7, 11..17, 29, 36, 37, 40) hangs on two
# transformers; local plants keep its N-1 cases solvable
plants.update({36: 120.0, 4: 120.0, 27: 100.0, 19: 100.0, 12: 100.0, 10: 250.0, 26: 200.0, 8: 150.0})
# small units become synchronous condensers (P = 0, reactive support only),
# their output is shifted to the remaining plants, as in the IEEE case with
# its many condenser buses
_small = [b for b, p in plants.items() if 0 < p < 45.0]
_moved = sum(plants[b] for b in _small)
for b in _small:
  plants[b] = 0.0
_big = [b for b, p in plants.items() if p > 0]
_scale = 1.0 + _moved / sum(plants[b] for b in _big)
for b in _big:
  plants[b] = round(plants[b] * _scale)
SLACK = 69
PLANT_138 = [b for b, p in plants.items() if b not in EHV_SET and (p >= 100.0 or b == SLACK)]
for pb in PLANT_138:
  deg = sum(1 for f, t, k in branches if pb in (f, t))
  cands = sorted((dist(pb, b), b) for b in range(1, N + 1) if b != pb and b not in EHV_SET)
  for _, b in cands:
    if deg >= 3:
      break
    if add_edge(pb, b, "L138"):
      deg += 1
# every bus but a handful of deliberate radial stubs gets at least two
# outlets (the IEEE case has about eight radial buses)
RADIAL = {1, 20, 35, 47, 58, 86, 97, 117}
for pb in range(1, N + 1):
  if pb in EHV_SET or pb in RADIAL:
    continue
  deg = sum(1 for f, t, k in branches if pb in (f, t))
  cands = sorted((dist(pb, b), b) for b in range(1, N + 1) if b != pb and b not in EHV_SET and b not in RADIAL)
  for _, b in cands:
    if deg >= 2:
      break
    if add_edge(pb, b, "L138"):
      deg += 1
# reinforcement of the north-west pocket: the 345/138 kV transformer bus 36
# gets a second outlet, the 17-29 corridor a parallel path, bus 11 a third
# neighbour
add_edge(36, 29, "L138")
add_edge(17, 5, "L138")
for _, b in sorted((dist(11, b), b) for b in range(1, N + 1) if b not in (11, 1, 12) and b not in EHV_SET):
  if add_edge(11, b, "L138"):
    break

def n_bridges(brs):
  """Number of bridges (edges whose outage islands a bus) by Tarjan's low-link."""
  adj = {b: [] for b in range(1, N + 1)}
  for idx, (f, t, k) in enumerate(brs):
    adj[f].append((t, idx)); adj[t].append((f, idx))
  disc = {}; low = {}; count = [0]; timer = [0]
  import sys
  sys.setrecursionlimit(10000)
  def dfs(u, pidx):
    disc[u] = low[u] = timer[0]; timer[0] += 1
    for v, idx in adj[u]:
      if idx == pidx:
        continue
      if v in disc:
        low[u] = min(low[u], disc[v])
      else:
        dfs(v, idx)
        low[u] = min(low[u], low[v])
        if low[v] > disc[u]:
          count[0] += 1
  for s in range(1, N + 1):
    if s not in disc:
      dfs(s, -1)
  return count[0]

n_lines_138 = sum(1 for b in branches if b[2] == "L138")
n_lines_345 = sum(1 for b in branches if b[2] == "L345")
n_trafo = sum(1 for b in branches if b[2] == "T")
# trim or top up 138 kV chords to exactly 186 branches
target_138 = 186 - n_lines_345 - n_trafo
all_pairs = sorted(((dist(u, v), u, v) for u in range(1, N + 1) for v in range(u + 1, N + 1)
                    if u not in EHV_SET and v not in EHV_SET), key=lambda p: p[0])
while n_lines_138 < target_138:
  for d, u, v in all_pairs:
    if (min(u, v), max(u, v)) not in edge_set and d < 30.0:
      add_edge(u, v, "L138")
      n_lines_138 += 1
      break
while n_lines_138 > target_138:
  # drop the longest 138 kV chord whose removal keeps both ends at degree >= 3
  deg = {}
  for f, t, k in branches:
    deg[f] = deg.get(f, 0) + 1
    deg[t] = deg.get(t, 0) + 1
  cands = sorted(((dist(f, t), i) for i, (f, t, k) in enumerate(branches) if k == "L138" and deg[f] >= 3 and deg[t] >= 3 and f not in PLANT_138 and t not in PLANT_138), reverse=True)
  removed = False
  for _, i in cands:
    trial = [b for j, b in enumerate(branches) if j != i]
    adj = {b: set() for b in range(1, N + 1)}
    for f, t, k in trial:
      adj[f].add(t); adj[t].add(f)
    seen = {1}; stack = [1]
    while stack:
      u = stack.pop()
      for v in adj[u]:
        if v not in seen:
          seen.add(v); stack.append(v)
    if len(seen) == N and n_bridges(trial) <= n_bridges(branches):
      f, t, _ = branches.pop(i)
      edge_set.discard((min(f, t), max(f, t)))
      n_lines_138 -= 1
      removed = True
      break
  assert removed, "cannot trim without islanding"
assert len(branches) == 186, len(branches)

# ------------------------------------------------------------ parameters --
BASE = 100.0
def line_params(f, t, kind):
  d = dist(f, t)
  if kind == "L138":
    km = 6.0 + 1.1 * d                    # 6 .. ~40 km
    zb = 138.0 ** 2 / BASE
    r = 0.10 * km / zb
    x = 0.40 * km / zb
    b = 2.8e-6 * km * zb
    rate = 300.0 if km < 30 else 250.0
    return round(r, 5), round(x, 5), round(b, 5), rate, 0.0, 0.0
  if kind == "L345":
    km = 40.0 + 3.0 * d                   # 80 .. ~150 km
    zb = 345.0 ** 2 / BASE
    r = 0.03 * km / zb
    x = 0.33 * km / zb
    b = 5.0e-6 * km * zb
    return round(r, 5), round(x, 5), round(b, 4), 900.0, 0.0, 0.0
  # transformer 345/138, 500 MVA, uk 12 %, tap on the EHV side
  x = 0.12 * BASE / 500.0
  ratio = round(rng.choice([0.96, 0.97, 0.98, 0.985, 1.0]), 3)
  return 0.0, round(x, 4), 0.0, 500.0, ratio, 0.0

branch_rows = []
for f, t, kind in branches:
  r, x, b, rate, ratio, angle = line_params(f, t, kind)
  branch_rows.append([f, t, r, x, b, rate, rate, rate, ratio, angle, 1, -360.0, 360.0])

# --------------------------------------------------------------- generators --
# 54 generator buses: 20 real plants (P > 0) and 34 synchronous condensers
# (P = 0, reactive support only), as in the IEEE case
gen_buses_all = [1, 4, 8, 10, 12, 15, 19, 24, 25, 26, 27, 30, 31, 32, 34, 36, 38, 40,
                 42, 46, 49, 54, 55, 56, 59, 61, 62, 63, 65, 66, 68, 69, 70, 72, 73, 74, 76, 77,
                 80, 81, 85, 87, 89, 90, 91, 92, 99, 100, 103, 104, 105, 107, 110, 111, 112, 113, 116]
for _drop in (1, 15, 112):
  gen_buses_all.remove(_drop)
assert len(gen_buses_all) == 54 and len(set(gen_buses_all)) == 54
SLACK = 69
gen_rows = []
for b in gen_buses_all:
  p = plants.get(b, 0.0)
  if b in plants:
    pmax = max(100.0, round(p * 1.35 / 10.0) * 10.0) if b != SLACK else 805.0
    qmax = round(pmax * 0.5)
    qmin = -round(pmax * 0.3)
    vg = rng.choice([1.0, 1.005, 1.01, 1.015, 1.02])
  else:
    pmax = 100.0
    qmax = rng.choice([24.0, 30.0, 40.0, 50.0, 60.0, 100.0])
    qmin = -rng.choice([8.0, 13.0, 20.0, 25.0, 30.0, 50.0])
    vg = rng.choice([0.985, 0.99, 0.995, 1.0, 1.005, 1.01])
  gen_rows.append([b, p, 0.0, qmax, qmin, vg, BASE, 1, pmax, 0.0])

# ------------------------------------------------------------------- loads --
load_buses = [b for b in range(1, N + 1) if b not in EHV_SET or b in (10, 30)]
load_buses = load_buses[:99] if len(load_buses) > 99 else load_buses
assert len(load_buses) == 99, len(load_buses)
raw = [min(rng.lognormvariate(3.3, 0.6), 95.0) for _ in load_buses]
scale = 4242.0 / sum(raw)
pd = {b: round(w * scale) for b, w in zip(load_buses, raw)}
qd = {b: round(pd[b] * rng.uniform(0.25, 0.45)) for b in load_buses}
# shunts: a few capacitor banks / reactors like the IEEE case
shunts = {5: (0.0, -40.0), 34: (0.0, 14.0), 37: (0.0, -25.0), 44: (0.0, 10.0), 45: (0.0, 10.0),
          46: (0.0, 10.0), 48: (0.0, 15.0), 74: (0.0, 12.0), 79: (0.0, 20.0), 82: (0.0, 20.0),
          83: (0.0, 10.0), 105: (0.0, 20.0), 107: (0.0, 6.0), 110: (0.0, 6.0)}

# ------------------------------------------------------------------- names --
prefixes = ["Alt", "Neu", "Hoch", "Nieder", "Ober", "Unter", "Wald", "Berg", "Tal", "Feld", "Bach",
            "See", "Stein", "Hain", "Lind", "Eich", "Birk", "Buch", "Wiesen", "Moor", "Sand", "Kies",
            "Rot", "Weiss", "Gruen", "Blau", "Gold", "Silber", "Nord", "Sued", "Ost", "West", "Kalt",
            "Warm", "Klein", "Gross", "Frei", "Fern", "Nah", "Lang", "Kurz", "Breit", "Tief", "Flach"]
suffixes = ["dorf", "heim", "hausen", "feld", "berg", "bach", "tal", "hof", "brueck", "furt",
            "kirchen", "stadt", "au", "roda", "walde", "hagen", "wang", "loh", "ried", "eck"]
names = []
seen = set()
while len(names) < N:
  n = rng.choice(prefixes) + rng.choice(suffixes)
  if n in seen:
    continue
  seen.add(n)
  names.append(n)
bus_name = {b: f"{names[b - 1]}_{345 if b in EHV_SET else 138}" for b in range(1, N + 1)}

# ------------------------------------------------------------------- buses --
bus_rows = []
for b in range(1, N + 1):
  area = next(a for a, r in AREAS.items() if b in r)
  gs, bs = shunts.get(b, (0.0, 0.0))
  btype = 3 if b == SLACK else (2 if b in gen_buses_all else 1)
  base_kv = 345.0 if b in EHV_SET else 138.0
  bus_rows.append([b, btype, float(pd.get(b, 0)), float(qd.get(b, 0)), gs, bs, area, 1.0, 0.0, base_kv, area, 1.06, 0.94])

# ------------------------------------------------------- validation: NR PF --
def ybus(branch_rows, bus_rows):
  Y = np.zeros((N, N), dtype=complex)
  for f, t, r, x, b, *_rest in branch_rows:
    ratio, angle, status = _rest[3], _rest[4], _rest[5]
    if not status:
      continue
    ys = 1.0 / complex(r, x)
    tap = (ratio if ratio else 1.0) * np.exp(1j * math.radians(angle))
    i, j = f - 1, t - 1
    Y[i, i] += (ys + 1j * b / 2) / (abs(tap) ** 2)
    Y[j, j] += ys + 1j * b / 2
    Y[i, j] += -ys / np.conj(tap)
    Y[j, i] += -ys / tap
  for row in bus_rows:
    i = row[0] - 1
    Y[i, i] += complex(row[4], row[5]) / BASE
  return Y

def connected(branch_rows):
  adj = {b: set() for b in range(1, N + 1)}
  for f, t, *_r in branch_rows:
    if _r[-3]:
      adj[f].add(t); adj[t].add(f)
  seen = {1}; stack = [1]
  while stack:
    u = stack.pop()
    for v in adj[u]:
      if v not in seen:
        seen.add(v); stack.append(v)
  return len(seen) == N

def newton(branch_rows, bus_rows, gen_rows, enforce_q=True, max_iter=30, tol=1e-8, warm=None):
  Y = ybus(branch_rows, bus_rows)
  G, B = Y.real, Y.imag
  Pd = np.array([r[2] for r in bus_rows]) / BASE
  Qd = np.array([r[3] for r in bus_rows]) / BASE
  Pg = np.zeros(N); Qg = np.zeros(N)
  vset = {}; qmax = {}; qmin = {}
  for g in gen_rows:
    i = g[0] - 1
    Pg[i] += g[1] / BASE
    vset[i] = g[5]; qmax[i] = g[3] / BASE; qmin[i] = g[4] / BASE
  slack = SLACK - 1
  pv = set(vset) - {slack}
  pq = set(range(N)) - pv - {slack}
  V = np.ones(N); th = np.zeros(N)
  for i, v in vset.items():
    V[i] = v
  if warm is not None:
    V = warm[0].copy(); th = warm[1].copy()
  clamped = {}
  for outer in range(12):
    for it in range(max_iter):
      Vc = V * np.exp(1j * th)
      S = Vc * np.conj(Y @ Vc)
      P, Q = S.real, S.imag
      Pspec = Pg - Pd
      Qspec = Qg - Qd
      pvl = sorted(pv); pql = sorted(pq)
      nonslack = sorted(pv | pq)
      dP = Pspec[nonslack] - P[nonslack]
      dQ = Qspec[pql] - Q[pql]
      mis = np.concatenate([dP, dQ])
      if np.max(np.abs(mis)) < tol:
        break
      # Jacobian (polar)
      n1 = len(nonslack); n2 = len(pql)
      J = np.zeros((n1 + n2, n1 + n2))
      Vm = V
      for a, i in enumerate(nonslack):
        for c, k in enumerate(nonslack):
          if i == k:
            J[a, c] = -Q[i] - B[i, i] * Vm[i] ** 2
          else:
            J[a, c] = Vm[i] * Vm[k] * (G[i, k] * math.sin(th[i] - th[k]) - B[i, k] * math.cos(th[i] - th[k]))
        for c, k in enumerate(pql):
          if i == k:
            J[a, n1 + c] = P[i] / Vm[i] + G[i, i] * Vm[i]
          else:
            J[a, n1 + c] = Vm[i] * (G[i, k] * math.cos(th[i] - th[k]) + B[i, k] * math.sin(th[i] - th[k]))
      for a, i in enumerate(pql):
        for c, k in enumerate(nonslack):
          if i == k:
            J[n1 + a, c] = P[i] - G[i, i] * Vm[i] ** 2
          else:
            J[n1 + a, c] = -Vm[i] * Vm[k] * (G[i, k] * math.cos(th[i] - th[k]) + B[i, k] * math.sin(th[i] - th[k]))
        for c, k in enumerate(pql):
          if i == k:
            J[n1 + a, n1 + c] = Q[i] / Vm[i] - B[i, i] * Vm[i]
          else:
            J[n1 + a, n1 + c] = Vm[i] * (G[i, k] * math.sin(th[i] - th[k]) - B[i, k] * math.cos(th[i] - th[k]))
      dx = np.linalg.solve(J, mis)
      for a, i in enumerate(nonslack):
        th[i] += dx[a]
      for a, i in enumerate(pql):
        V[i] += dx[n1 + a]
    else:
      return None
    if not enforce_q:
      break
    Vc = V * np.exp(1j * th)
    S = Vc * np.conj(Y @ Vc)
    Qgen = S.imag + Qd
    switched = False
    for i in sorted(pv):
      if Qgen[i] > qmax[i] + 1e-9:
        pv.discard(i); pq.add(i); Qg[i] = qmax[i]; clamped[i] = "max"; switched = True
      elif Qgen[i] < qmin[i] - 1e-9:
        pv.discard(i); pq.add(i); Qg[i] = qmin[i]; clamped[i] = "min"; switched = True
    if not switched:
      break
  Vc = V * np.exp(1j * th)
  S = Vc * np.conj(Y @ Vc)
  # branch loadings
  loads = []
  for f, t, r, x, b, rateA, *_rest in branch_rows:
    ratio, angle, status = _rest[2], _rest[3], _rest[4]
    if not status:
      loads.append(0.0); continue
    ys = 1.0 / complex(r, x)
    tap = (ratio if ratio else 1.0) * np.exp(1j * math.radians(angle))
    i, j = f - 1, t - 1
    If = (Vc[i] / tap - Vc[j]) * ys / np.conj(tap) + Vc[i] / (abs(tap) ** 2) * 1j * b / 2
    It = (Vc[j] - Vc[i] / tap) * ys + Vc[j] * 1j * b / 2
    sf = abs(Vc[i] * np.conj(If)) * BASE
    st = abs(Vc[j] * np.conj(It)) * BASE
    loads.append(100.0 * max(sf, st) / rateA)
  return dict(V=V, th=th, S=S, iters=it + 1, clamped=clamped, loads=np.array(loads), pslack=(S.real[slack] + Pd[slack]) * BASE)

assert connected(branch_rows), "base case not connected"
res = newton(branch_rows, bus_rows, gen_rows)
res0 = newton(branch_rows, bus_rows, gen_rows, enforce_q=False)
Qg0 = (res0["S"].imag + np.array([r[3] for r in bus_rows]) / BASE) * BASE
viol = sorted(((Qg0[g[0]-1] - g[3]) if Qg0[g[0]-1] > g[3] else (Qg0[g[0]-1] - g[4]) if Qg0[g[0]-1] < g[4] else 0.0, g[0], round(Qg0[g[0]-1],1), g[4], g[3]) for g in gen_rows if g[0] != SLACK)
print("unlimited pass: worst Q excursions", [v for v in viol if v[0] != 0.0][:6], "...", [v for v in viol if v[0] != 0.0][-6:])
# reactive limits derived from the unlimited solution: generous for most
# machines, deliberately BELOW the need for ten of them so that the case has
# machines at their limits (the IEEE case has several), rounded to 5 Mvar
CLAMP = {4, 19, 31, 32, 54, 72, 73, 77, 85, 87}
def r5(v):
  return 5.0 * round(v / 5.0)
for g in gen_rows:
  b = g[0]
  if b == SLACK:
    continue
  q0 = Qg0[b - 1]
  if b in CLAMP:
    if q0 >= 0:
      g[3] = max(10.0, r5(0.75 * q0)); g[4] = -max(10.0, r5(0.4 * abs(q0)) )
    else:
      g[4] = min(-10.0, r5(0.75 * q0)); g[3] = max(10.0, r5(0.4 * abs(q0)))
  else:
    g[3] = max(20.0, r5(1.3 * max(q0, 0.0) + 20.0), r5(0.3 * g[8]) if g[1] > 0 else 0.0)
    g[4] = min(-20.0, -r5(1.3 * max(-q0, 0.0) + 20.0))
res = newton(branch_rows, bus_rows, gen_rows)
assert res is not None, "base case did not converge (limits)"
# ratings sized to the base-case flow with 30 % headroom, in the usual
# classes, so that the base case is healthy and N-1 produces real overloads
CLASSES = {"L138": [200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 800.0], "L345": [900.0, 1200.0, 1500.0], "T": [500.0, 750.0, 1000.0]}
for i, (f, t, kind) in enumerate(branches):
  flow = res["loads"][i] / 100.0 * branch_rows[i][5]
  cls = CLASSES[kind]
  rate = next((c for c in cls if c >= 1.3 * flow), cls[-1])
  branch_rows[i][5] = branch_rows[i][6] = branch_rows[i][7] = rate
res = newton(branch_rows, bus_rows, gen_rows)
assert res is not None, "base case did not converge"
V = res["V"]
print(f"base case: {res['iters']} NR iterations, Vmin {V.min():.4f} at bus {V.argmin()+1}, Vmax {V.max():.4f}, "
      f"max loading {res['loads'].max():.1f} %, slack P {res['pslack']:.1f} MW, {len(res['clamped'])} generators at a Q-limit")
for a, r in AREAS.items():
  print(f"area {a}: load {sum(pd.get(b,0) for b in r):.0f} MW, plant P {sum(plants.get(b,0) for b in r):.0f} MW")
over = sorted(((l, i) for i, l in enumerate(res["loads"]) if l > 80), reverse=True)[:8]
print("overloaded:", [(branch_rows[i][0], branch_rows[i][1], branches[i][2], round(l)) for l, i in over])
low = sorted(((v, i+1) for i, v in enumerate(V)))[:8]
print("low V:", [(b, round(v,3)) for v, b in low])
assert 0.94 <= V.min() and V.max() <= 1.06, "voltage band"
assert res["loads"].max() < 90.0, "base overload"
assert 4 <= len(res["clamped"]) <= 20, "Q-limit count out of the intended range"

# N-1 count and the workshop's double outage (branches 1 and 2)
n_gen_on = sum(1 for g in gen_rows if g[7])
print(f"N-1 all: {len(branch_rows) + n_gen_on} (branches {len(branch_rows)}, generators {n_gen_on})")
double = [list(r) for r in branch_rows]
double[0][10] = 0; double[1][10] = 0
assert connected(double), "double outage of branches 1 and 2 islands the grid"
rd = newton(double, bus_rows, gen_rows)
assert rd is not None, "double outage did not converge"
print(f"double outage (branches 1, 2): converged in {rd['iters']} iterations, Vmin {rd['V'].min():.4f}")
n_island = 0; n_fail = 0
for k in range(len(branch_rows)):
  single = [list(r) for r in branch_rows]
  single[k][10] = 0
  if not connected(single):
    n_island += 1
    continue
  if newton(single, bus_rows, gen_rows, warm=(res['V'], res['th'])) is None:
    n_fail += 1
    print("   non-converging N-1:", branch_rows[k][0], branch_rows[k][1], branches[k][2], "loading", round(res["loads"][k]), "%")
print(f"N-1 branches: {n_island} islanding outages, {n_fail} non-converging")

# ------------------------------------------------------------------- write --
def fmt(v):
  if isinstance(v, float):
    s = f"{v:.6g}"
    return s if ("." in s or "e" in s) else s + ".0"
  return str(v)

out = []
out.append("function mpc = sp_case118")
out.append("% SP_CASE118  Synthetic 118-bus case in the shape of the IEEE 118-bus system.")
out.append("%")
out.append("% Same cardinalities as the IEEE 118-bus case (118 buses, 186 branches of")
n_cond = sum(1 for g in gen_rows if g[1] == 0.0 and g[0] != SLACK)
n_loads = sum(1 for r in bus_rows if r[2] > 0)
out.append(f"% which 9 transformers, 54 generators of which {n_cond} synchronous condensers,")
out.append(f"% {n_loads} loads, three areas, an 11-bus 345 kV backbone over a 138 kV grid),")
out.append("% with its OWN topology, parameters and names: no IEEE data is copied, the")
out.append("% file carries no external license. Generated by tools/gen_sp_case118.py")
out.append(f"% (seed {SEED}); validated there with a Newton-Raphson power flow including")
out.append("% generator Q-limit enforcement, an island check and an N-1 sweep.")
out.append("%")
out.append(f"% Base case: {res['iters']} NR iterations, Vmin {V.min():.4f} pu, Vmax {V.max():.4f} pu,")
out.append(f"% highest branch loading {res['loads'].max():.1f} %, {len(res['clamped'])} generators at a")
out.append(f"% reactive limit, total load {sum(pd.values())} MW / {sum(qd.values())} Mvar. N-1 over the")
out.append(f"% branches: {n_island} outages island a bus, {n_fail} do not converge. Bus 69 is the slack.")
out.append("% Ratings: sized per branch to the base-case flow with 30 percent headroom in the")
out.append("% classes 200..800 MVA (138 kV), 900..1500 MVA (345 kV) and 500..1000 MVA (transformers).")
out.append("")
out.append("%% MATPOWER Case Format : Version 2")
out.append("mpc.version = '2';")
out.append("")
out.append("%%-----  Power Flow Data  -----%%")
out.append("%% system MVA base")
out.append("mpc.baseMVA = 100;")
out.append("")
out.append("%% bus data")
out.append("%\tbus_i\ttype\tPd\tQd\tGs\tBs\tarea\tVm\tVa\tbaseKV\tzone\tVmax\tVmin")
out.append("mpc.bus = [")
for r in bus_rows:
  out.append("\t" + "\t".join(fmt(v) for v in r) + ";")
out.append("];")
out.append("")
out.append("%% generator data")
out.append("%\tbus\tPg\tQg\tQmax\tQmin\tVg\tmBase\tstatus\tPmax\tPmin")
out.append("mpc.gen = [")
for r in gen_rows:
  out.append("\t" + "\t".join(fmt(v) for v in r) + ";")
out.append("];")
out.append("")
out.append("%% branch data")
out.append("%\tfbus\ttbus\tr\tx\tb\trateA\trateB\trateC\tratio\tangle\tstatus\tangmin\tangmax")
out.append("mpc.branch = [")
for r in branch_rows:
  out.append("\t" + "\t".join(fmt(v) for v in r) + ";")
out.append("];")
out.append("")
out.append("%% bus names (imported with matpower_import.apply_bus_names)")
out.append("mpc.bus_name = {")
for b in range(1, N + 1):
  out.append(f"\t'{bus_name[b]}';")
out.append("};")
out.append("")
import os, sys
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "mpower", "sp_case118.m")
with open(OUT, "w", newline="\n") as fh:
  fh.write("\n".join(out))
print("written", OUT)